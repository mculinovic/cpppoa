// Copyright 2015 @mculinovic
#include <string>

#include "./graph.hpp"
#include "./alignment.hpp"

using std::string;

/**
 * @brief poa test program
 * @details
 */
int main() {
    string seq("PESLLYGRFTIESDVW");
    string label("test1");

    Graph graph(seq, label);

    string seq2("PEAALYGRFTIKSDVW");
    string label2("test2");

    Alignment aln(seq2, graph);
    aln.align();
    graph.insertSequenceAlignment(aln, seq2, label2);

    string seq3("PESLAYNKFSIKSDVW");
    string label3("test3");

    Alignment aln2(seq3, graph);
    aln2.align();
    graph.insertSequenceAlignment(aln2, seq3, label3);

    string seq4("PEALNYGRYSSESDVW");
    string label4("test4");

    Alignment aln3(seq4, graph);
    aln3.align();
    graph.insertSequenceAlignment(aln3, seq4, label4);

    string seq5("PEVIRMQDDNPFSFQSDVY");
    string label5("test5");

    Alignment aln4(seq5, graph);
    aln4.align();
    graph.insertSequenceAlignment(aln4, seq5, label5);
    graph.alignment_strings();
    return 0;
}
