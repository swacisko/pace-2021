//
// Created by sylwester on 3/12/21.
//

#ifndef ALGORITHMSPROJECT_TEST_GRAPHS_H
#define ALGORITHMSPROJECT_TEST_GRAPHS_H

#include "Makros.h"

class CE_test_graphs{
public:
    /**
     * Given a list of edges, creates the VVI structure
     * @param edges
     * @return
     */
    static VVI createVVIfromVPII(VPII & edges);

    static VVI flowcutter_test2;
    static VVI flowcutter_test3;

    static VPII swpcndeo_edges;
    static VI swpcndeo_partition;

    static VVI cluster_graph_test;
    static VI cluster_graph_test_partition;

    static VPII cluster_graph_test2_edges;
    static VI cluster_graph_test2_partition;

    static VPII neg_interchange_edges;
    static VI neg_interchange_partition;

    static VPII comp_exp_rep_edges;
    static VI comp_exp_rep_partition;

    static VVI pathP9;
    static VVI pathP3;

    static VVI big_test_50_104;

    /**
     * Three disjoint edges, one triangle and one P3 (all 5 structures disjoint).
     */
    static VVI kern1;
    static VVI kern_ed;
    static VVI kern2;
    static VVI kern3;

    static VVI kern4positive;
    static VVI kern4negative;

    static VVI kern6positive;
    static VVI kern6negative;

    static VVI kern7positive;
    static VVI kern7negative;

    static VVI kern8positive;
    static VVI kern8positive2;
    static VVI kern8negative;

    static VVI kern9positive;
    static VVI kern9negative;

    static VVI kern15;

    static VVI kern_heur_1;
    static VVI kern_heur_2;
    static VVI kern_heur_2_neg;

    static VVI kern_heur_3;
    static VVI kern_heur_4;
    static VVI kern_heur_5;
    static VVI kern_heur_6;
    static VPII kern_heur_7_edges;


    static vector<tuple<int,int,int>> EOgraphEdges; // edges (v,w,weight) in the graph - easier to write
    static VI EOnw; // node weights of the EOgraph

    static VPII swpcndedge_test1_edges;
    static VI swpcndedge_test1_partition;

    static VPII swpcndedge_test2_edges;
    static VI swpcndedge_test2_partition;

};

#endif //ALGORITHMSPROJECT_TEST_GRAPHS_H
