//
// Created by sylwester on 3/16/21.
//

#ifndef ALGORITHMSPROJECT_CLUSTERGRAPH_H
#define ALGORITHMSPROJECT_CLUSTERGRAPH_H

#include <graphs/GraphInducer.h>
#include "Makros.h"


class ClusterGraph{
public:

    ClusterGraph() : origV(nullptr) {}

    /**
     * Creates cluster graph for given partition
     * If partition[i] == -1, then the node i will not be considered for graph construction.
     * @param V
     * @param partition
     */
    ClusterGraph( VVI * V, VI partition );

    /**
     * Just an abbreviation
     * @param x
     * @return reference to node_weights[x]
     */
    int& nw(int x){ return node_weights[x]; }

    static void test();

    friend ostream& operator<<( ostream& str, ClusterGraph& g );

    /**
     * Creates a vector of pairs {a,b} such that there are exactly b nodes with nw(b) == a
     * @param cnt
     */
    VPII getClusterSizesCount();

//private:

    /**
     * Pointer to the original graph
     */
    VVI * origV;
//    VI origPartition;

    /**
     * Structure of the graph. V[i][j] = ( id, weight )
     */
    VVPII V;
    int N; // V.size()

    VI node_weights;

    /**
     * clusterNodes[i] is the list of nodes from the original graph that are in a cluster represented by node i
     */
    VVI clusterNodes;

    /**
     * partition[i] is the reverse of clusterNodes
     */
    VI partition;

    /**
     * Creates cluster graph structure
     */
    void createGraph();

};

class InducedClusterGraph : public InducedGraphPI{
public:
    /**
     * Default constructor. Setting par to nullptr.
     */
    InducedClusterGraph() : clPar(nullptr) {}

    /**
     * Creates structure of the induced graph - swaps with those of [ig].
     * @param ig
     */
    InducedClusterGraph( ClusterGraph& clg, InducedGraphPI ig );

    friend ostream& operator<<( ostream& str, InducedClusterGraph& cl );

    VI node_weights;

    /**
     * Pointer to the cluster graph parent
     */
    ClusterGraph* clPar;
};

/**
 * Induces the graph
 * @param clg
 * @param nodes
 * @return
 */
InducedClusterGraph induce( ClusterGraph & clg, VI & nodes );

/**
 * Induces the graph, WITHOUT creating [perm] map. This is faster, but does not admit mapping nodes from original
 * graph to present graph.
 * @param clg
 * @param nodes
 * @param helper
 * @return
 */
InducedClusterGraph induce( ClusterGraph & clg, VI & nodes, VI & helper );

#endif //ALGORITHMSPROJECT_CLUSTERGRAPH_H
