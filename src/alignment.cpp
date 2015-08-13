/**
 * @file alignment.cpp
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @brief Alignment class implementation file 
 * @details Implementation file for Alignment class used for
 * calculating local alignment between graph and sequence. Class is
 * based on https://github.com/ljdursi/poapy/blob/master/seqgraphalignment.py
 * python implementation
 */
#include <algorithm>
#include <vector>
#include <string>
#include <tuple>
#include <limits>
#include <iostream>
#include <cassert>

#include "./alignment.hpp"
#include "./graph.hpp"

using std::get;
using std::max;
using std::numeric_limits;
using std::string;
using std::tuple;
using std::vector;

int Alignment::match_score_ = 4;
int Alignment::mismatch_score_ = -2;
int Alignment::open_gap_score_ = -4;
int Alignment::extend_gap_score_ = -2;


Alignment::Alignment(const string& sequence,
                     const Graph& graph): sequence_(sequence),
                                          graph_(graph) {
    assert(sequence.size() > 0);
}

/**
 * @brief returns move score
 * @details returns move score
 *
 * @param move move
 * @return move score
 */
int move_score(const move& move) {
    return get<0>(move);
}


void Alignment::backtrack() {
    while (scores_[max_i_][max_j_] > 0 && !(max_i_ == 0 && max_j_ == 0)) {
        int next_i = backtrack_graph_idx_[max_i_][max_j_];
        int next_j = backtrack_seq_idx_[max_i_][max_j_];
        uint32_t seq_idx = max_j_ - 1;
        uint32_t node_id = index_to_nodeID_[max_i_ - 1];

        seq_idxs_.emplace_front(next_j != max_j_ ? seq_idx : -1);
        node_ids_.emplace_front(next_i != max_i_ ? node_id : -1);

        max_i_ = next_i;
        max_j_ = next_j;
    }
}


void Alignment::align() {
    uint32_t n = graph_.getNodesNum();
    uint32_t m = sequence_.length();
    max_score_ = numeric_limits<int>::min();
    max_i_ = -1;
    max_j_ = -1;

    init_dp_tables();

    vector<int> insert_cost((n+1)*(m+1), open_gap_score_);
    vector<int> delete_cost((n+1)*(m+1), open_gap_score_);

    auto match_score = [](char seq_base, char node_base) -> int {
        if (seq_base == node_base) {
            return match_score_;
        } else {
            return mismatch_score_;
        }
    };

    // modified smith-waterman(SW)
    // operation candidates
    // char values:
    // -> 'I' - insertion
    // -> 'D' - deletion
    // -> 'M' - match or mismatch
    const vector<uint32_t>& nodes_ids = const_cast<Graph&>(graph_).getNodesIds();
    for (uint32_t i = 0; i < n; ++i) {
        char base = graph_.getNode(nodes_ids[i])->base();

        for (uint32_t j = 0; j < sequence_.length(); ++j) {
            // insertion from sequence, unchanged as in SW
            move best_candidate(scores_[i + 1][j] + insert_cost[(i + 1)*m + j],
                                    i + 1, j, 'I');

            auto& node = graph_.getNode(nodes_ids[i]);
            if (node->getPredecessorsIds().size()) {
              // for every other operation I have to check for all predeccesors of
              // current node
              for (auto node_id : node->getPredecessorsIds()) {
                auto prev_index = nodeID_to_index_[node_id];

                // match/mismatch
                move mm(scores_[prev_index + 1][j] +
                    match_score(sequence_[j], base),
                    prev_index + 1, j, 'M');
                if (move_score(mm) >= move_score(best_candidate)) {
                  best_candidate = mm;
                }

                // insertion from graph to sequence / deletion
                move ins(scores_[prev_index + 1][j + 1] +
                    delete_cost[(prev_index + 1)*m + j + 1],
                    prev_index + 1, j + 1, 'D');
                if (move_score(ins) >= move_score(best_candidate)) {
                  best_candidate = ins;
                }
              }
            } else {
              int prev_index = -1;
              // match/mismatch
              move mm(scores_[prev_index + 1][j] +
                  match_score(sequence_[j], base),
                  prev_index + 1, j, 'M');
              if (move_score(mm) >= move_score(best_candidate)) {
                best_candidate = mm;
              }

              // insertion from graph to sequence / deletion
              move ins(scores_[prev_index + 1][j + 1] +
                  delete_cost[(prev_index + 1)*m + j + 1],
                  prev_index + 1, j + 1, 'D');
              if (move_score(ins) >= move_score(best_candidate)) {
                best_candidate = ins;
              }
            }

            // update dp and backtrack tables
            scores_[i + 1][j + 1] = get<0>(best_candidate);
            backtrack_graph_idx_[i + 1][j + 1] = get<1>(best_candidate);
            backtrack_seq_idx_[i + 1][j + 1] = get<2>(best_candidate);

            char move = get<3>(best_candidate);
            if (move == 'I') {
                insert_cost[(i + 1) * m + (j + 1)] = extend_gap_score_;
            } else if (move == 'D') {
                delete_cost[(i + 1) * m + (j + 1)] = extend_gap_score_;
            }

            // because of local alignment minimum score should be zero
            if (scores_[i + 1][j + 1] < 0) {
                scores_[i + 1][j + 1] = 0;
                backtrack_graph_idx_[i + 1][j + 1] = -1;
                backtrack_seq_idx_[i + 1][j + 1] = -1;
            }


            // update max score and its position
            if (scores_[i + 1][j + 1] >= max_score_) {
                max_i_ = i + 1;
                max_j_ = j + 1;
                max_score_ = scores_[i + 1][j + 1];
            }
        }
    }

    backtrack();
}


void Alignment::init_dp_tables() {
    uint32_t n = graph_.getNodesNum();
    uint32_t m = sequence_.length();

    const vector<uint32_t>& nodes_ids = const_cast<Graph&>(graph_).getNodesIds();
    int max_node_id = -1;
    for (int node_id: nodes_ids) {
      max_node_id = max(max_node_id, node_id);
    }

    nodeID_to_index_.resize(max_node_id + 1, -1);
    index_to_nodeID_.resize(nodes_ids.size());

    for (size_t i = 0; i < nodes_ids.size(); ++i) {
        nodeID_to_index_[nodes_ids[i]] = i;
        index_to_nodeID_[i] = nodes_ids[i];
    }

    // init dynamic programming smith waterman matrix
    scores_.resize(n + 1, vector<int>(m + 1, 0));


    // init backtracking matrices
    backtrack_seq_idx_.resize(n + 1, vector<int>(m + 1, 0));
    backtrack_graph_idx_.resize(n + 1, vector<int>(m + 1, 0));
}
