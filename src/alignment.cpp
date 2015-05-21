/**
 * @file alignment.cpp
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @brief Alignment class implementation file 
 * @details Implementation file for Alignment class used for
 * calculating local alignment between graph and sequence. Class is
 * based on https://github.com/ljdursi/poapy/blob/master/seqgraphalignment.py
 * python implementation
 */
#include <vector>
#include <string>
#include <tuple>
#include <unordered_map>
#include <limits>

#include "./alignment.hpp"
#include "./graph.hpp"

using std::vector;
using std::string;
using std::unordered_map;
using std::tuple;
using std::get;
using std::numeric_limits;


int Alignment::match_score_ = 4;
int Alignment::mismatch_score_ = -2;
int Alignment::open_gap_score_ = -4;
int Alignment::extend_gap_score_ = -2;


Alignment::Alignment(const string& sequence,
                     const Graph& graph): sequence_(sequence),
                                          graph_(graph) {}


move Alignment::get_best_move(const vector<move>& candidates) {
    int max = numeric_limits<int>::min();
    move best_candidate;
    for (auto& candidate : candidates) {
        if (get<0>(candidate) >= max) {
            max = get<0>(candidate);
            best_candidate = candidate;
        }
    }
    return best_candidate;
}


vector<int> Alignment::node_predecessors(const shared_ptr<Node>& node) {
    vector<int> prev_indices;
    for (auto id : node->getPredecessorsIds()) {
        prev_indices.emplace_back(nodeID_to_index_[id]);
    }
    // if no predecessors, point to just before the graph
    if (prev_indices.empty()) {
        prev_indices.emplace_back(-1);
    }
    return prev_indices;
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

    // python defaultdict implemented with unordered_map and lambda function
    unordered_map<uint32_t, int> insert_map;
    auto insert_cost = [m, &insert_map](uint32_t i, uint32_t j) mutable -> int {
        if (insert_map.find(i * m + j) == insert_map.end()) {
            return open_gap_score_;
        } else {
            return insert_map[i * m + j];
        }
    };

    unordered_map<uint32_t, int> delete_map;
    auto delete_cost = [m, &delete_map](uint32_t i, uint32_t j) mutable -> int {
        if (delete_map.find(i * m + j) == delete_map.end()) {
            return open_gap_score_;
        } else {
            return delete_map[i * m + j];
        }
    };

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
            // candidates are presented by tuples, see comment for move in .hpp
            vector<move> candidates;

            // insertion from sequence, unchanged as in SW
            candidates.emplace_back(scores_[i + 1][j] + insert_cost(i + 1, j),
                                    i + 1, j, 'I');

            // for every other operation I have to check for all predeccesors of
            // current node
            auto& node = graph_.getNode(nodes_ids[i]);
            for (auto prev_id : node_predecessors(node)) {
                // match/mismatch
                candidates.emplace_back(scores_[prev_id + 1][j] +
                                        match_score(sequence_[j], base),
                                        prev_id + 1, j, 'M');
                // insertion from graph to sequence / deletion
                candidates.emplace_back(scores_[prev_id + 1][j + 1] +
                                        delete_cost(prev_id + 1, j + 1),
                                        prev_id + 1, j + 1, 'D');
            }

            auto best_candidate = get_best_move(candidates);

            // update dp and backtrack tables
            scores_[i + 1][j + 1] = get<0>(best_candidate);
            backtrack_graph_idx_[i + 1][j + 1] = get<1>(best_candidate);
            backtrack_seq_idx_[i + 1][j + 1] = get<2>(best_candidate);

            char move = get<3>(best_candidate);
            if (move == 'I') {
                insert_map[(i + 1) * m + (j + 1)] = extend_gap_score_;
            } else if (move == 'D') {
                delete_map[(i + 1) * m + (j + 1)] = extend_gap_score_;
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
