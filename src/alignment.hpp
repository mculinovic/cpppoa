/**
 * @file alignment.hpp
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @brief Alignment class header file 
 * @details Header file with declaration of Alignment class used for
 * calculating local alignment between graph and sequence. Class is
 * based on https://github.com/ljdursi/poapy/blob/master/seqgraphalignment.py
 * python implementation
 */
#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <vector>
#include <string>
#include <unordered_map>
#include <tuple>
#include <deque>

#include "./graph.hpp"

using std::vector;
using std::string;
using std::unordered_map;
using std::tuple;
using std::deque;


typedef tuple<int, int, int, char> move;

/**
 * @brief Alignment class
 * @details Calculates local alignment between graph nodes and sequence
 * using modified smith-waterman algorithm.
 * 
 */
class Alignment {
 public:
    /**
     * @brief Alignment constructor
     * @details Constructor needs two parameters - sequence and graph
     * for calculating alignment between those two.
     * 
     * @param sequence reference to string sequence
     * @param graph reference to graph
     */
    Alignment(const string& sequence, const Graph& graph);


    /**
     * @brief Aligns sequence to graph
     * @details Method performs local alignment between
     * string sequence and directed acyclic graph (DAG).
     * Algorithm is variation of smith-waterman algorithm 
     * (http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm).
     */
    void align();


    /**
     * @brief Getter for sequence
     * @return sequence
     */
    const string& sequence() const { return sequence_; }


    /**
     * @brief Getter for sequence indices in alignment
     * @details Method returns sequence indices representing
     * alignment. Index value is -1 if node in graph isn't aligned
     * to any position/character in sequence, otherwise it is position
     * in sequence string.
     * @return sequence indices
     */
    const deque<int>& seq_idxs() const { return seq_idxs_; }


    /**
     * @brief Getter for node ids in alignment
     * @details Methor returns node ids representing
     * alignment. Id value is -1 if some character in sequence
     * isn't aligned to any node in graph, otherwise it is id of node
     * in graph.
     * @return node ids
     */
    const deque<int>& node_ids() const { return node_ids_; }


 private:
    // value of match in smith-waterman score matrix
    static int match_score_;
    // value of mismatch in smith-waterman score matrix
    static int mismatch_score_;
    // value of gap opening in smith-waterman score matrix
    static int open_gap_score_;
    // value of extending gap in smith-waterman score matrix
    static int extend_gap_score_;

    // sequence to be aligned to graph
    const string& sequence_;
    // poa graph
    const Graph& graph_;

    // smith-waterman dynamic programming score matrix
    vector<vector<int>> scores_;
    // backtracking matrix for sequence
    vector<vector<int>> backtrack_seq_idx_;
    // backtracking matrix for graph
    vector<vector<int>> backtrack_graph_idx_;

    // sequence indices representing alignment
    deque<int> seq_idxs_;
    // node ids representing alignment
    deque<int> node_ids_;

    // map node ID to index in dp matrix
    unordered_map<uint32_t, uint32_t> nodeID_to_index_;
    // map index in dp matrix to node ID
    unordered_map<uint32_t, uint32_t> index_to_nodeID_;

    // max score in matrix
    int max_score_;
    // row index of maximum score
    int max_i_;
    // column index of maximum score
    int max_j_;


    /**
     * @brief Initializes matrices for smith-waterman
     * @details Initializes dynamic programming score
     * matrix, and backtracking matrices for keeping
     * track of best alignment path
     */
    void init_dp_tables();


    /**
     * @brief Backtracks best alignment path
     */
    void backtrack();


    /**
     * @brief Finds predecessors of node
     * @details Node predeccesors are nodes which outgoing
     * edges are ending in node given as function parameter
     * 
     * @param node find predeccsors of this node
     * @return predecessors node ids
     */
    vector<int> node_predecessors(const shared_ptr<Node>& node);


    /**
     * @brief calculates best move for dp matrix
     * @details calculates maximum score move between
     * candidates
     * 
     * @param candidates candidates for max score
     * @return max score move; type tuple (@see typedef for move)
     */
    move get_best_move(const vector<move>& candidates);
};

#endif  // ALIGNMENT_H
