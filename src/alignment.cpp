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
using std::min;
using std::numeric_limits;
using std::string;
using std::tuple;
using std::vector;

namespace POA {

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


    /**
     * @brief returns move graph_idx
     * @details returns move graph_idx
     *
     * @param move move
     * @return move graph_index
     */
    int move_graph_idx(const move& move) {
        return get<1>(move);
    }


    /**
     * @brief returns move seq_idx
     * @details returns move seq_idx
     *
     * @param move move
     * @return move seq_idx
     */
    int move_seq_idx(const move& move) {
        return get<2>(move);
    }


    /**
     * @brief returns move type
     * @details returns move type
     *
     * @param move move
     * @return move type
     */
    char move_type(const move& move) {
        return get<3>(move);
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
        align_banded_starting_at(0, -1);
    }

    void Alignment::align_starting_at(const uint32_t pos) {
        align_banded_starting_at(pos, -1);
    }

    void Alignment::align_banded_starting_at(const uint32_t min_pos, const int band_width) {
        uint32_t m = sequence_.length();
        max_score_ = numeric_limits<int>::min();
        max_i_ = -1;
        max_j_ = -1;

        init_dp_tables(min_pos);

        vector<int> insert_cost((valid_nodes_num_+1)*(m+1), open_gap_score_);
        vector<int> delete_cost((valid_nodes_num_+1)*(m+1), open_gap_score_);

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
        for (uint32_t i = 0; i < valid_nodes_num_; ++i) {
            const auto node_id = index_to_nodeID_[i];
            auto& node = graph_.getNode(node_id);
            const char base = node->base();

            // calculate sequence substring that should be part of DP
            int lo = 0, hi = sequence_.length();
            if (band_width >= 0) {
                lo = max(0, (int) (graph_.min_node_distance(node_id) - band_width - min_pos));
                hi = min((int) sequence_.length(), (int) (graph_.max_node_distance(node_id) + band_width - min_pos + 1));
            }

            for (uint32_t j = lo; j < hi; ++j) {
                // insertion from sequence, unchanged as in SW
                move best_candidate(scores_[i + 1][j] + insert_cost[(i + 1)*m + j],
                                        i + 1, j, 'I');

                if (node->getPredecessorsIds().size()) {
                  // for every other operation I have to check for all predeccesors of
                  // current node
                  for (auto node_id : node->getPredecessorsIds()) {
                    auto prev_index = index_from_node_id(node_id);

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
                scores_[i + 1][j + 1] = move_score(best_candidate);
                backtrack_graph_idx_[i + 1][j + 1] = move_graph_idx(best_candidate);
                backtrack_seq_idx_[i + 1][j + 1] = move_seq_idx(best_candidate);

                char move = move_type(best_candidate);
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


    void Alignment::init_dp_tables(const int min_pos) {
        uint32_t m = sequence_.length();
        valid_nodes_num_ = 0;

        max_valid_node_id_ = numeric_limits<int>::min();
        min_valid_node_id_ = numeric_limits<int>::max();

        const vector<uint32_t>& nodes_ids = const_cast<Graph&>(graph_).getNodesIds();
        for (int node_id: nodes_ids) {
          if (graph_.max_node_distance(node_id) >= min_pos) {
            valid_nodes_num_++;
            max_valid_node_id_ = max(max_valid_node_id_, node_id);
            min_valid_node_id_ = min(min_valid_node_id_, node_id);
          }
        }

        index_to_nodeID_.resize(valid_nodes_num_);
        nodeID_to_index_.resize(max_valid_node_id_ - min_valid_node_id_ + 1, -1);

        uint32_t idx = 0;
        for (int node_id: nodes_ids) {
            if (graph_.max_node_distance(node_id) < min_pos) {
                continue;
            }

            nodeID_to_index_[node_id - min_valid_node_id_] = idx;
            index_to_nodeID_[idx] = node_id;
            idx++;
        }

        // init dynamic programming smith waterman matrix
        scores_.resize(valid_nodes_num_ + 1, vector<int>(m + 1, 0));


        // init backtracking matrices
        backtrack_seq_idx_.resize(valid_nodes_num_ + 1, vector<int>(m + 1, -1));
        backtrack_graph_idx_.resize(valid_nodes_num_ + 1, vector<int>(m + 1, -1));
    }

    int Alignment::index_from_node_id(const uint32_t node_id) const {
        if (node_id < min_valid_node_id_ || node_id > max_valid_node_id_) {
            return -1;
        }
        return nodeID_to_index_[(int) node_id - min_valid_node_id_];
    }
}
