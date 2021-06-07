//
// Created by sylwester on 3/8/21.
//

#ifndef CESWAT_STATE_H
#define CESWAT_STATE_H

#include "Makros.h"
#include "Cluster.h"

enum StateInitializationType{
    SINGLE_NODES, // all nodes are in separate clusters
    RANDOM_MATCHING,  // we take random matching, then make cluster of size 2 for each edge in the matching
    ONE_CLUSTER, // one giant cluster
    SQRT_RANDOM, // sqrt(N) random subsets, each of size roughly sqrt(N)
    MAXIMUM_MATCHING, // take a largest matching in the graph
    RANDOM_STATE_PERM, // takes a random permutation, then for each nodes considers nodes that are before it, and joins it to the best cluster
    LEAF_TRIMMING, // takes a leaf and makes it in cluster with its neighbor, then 'removes both' nodes from graph and repeats
    EXPANSION_ORDER
};

/**
 * Represents 'global' data for current state - mapping arrays, etc.
 * CAUTION! In clusters, last position cluster.back() is ALWAYS AN EMPTY CLUSTER, for implementation reasons.
 */
class State{
public:

    /**
     * By default, RANDOM_MATCHING initialization method is preferred.
     * @param clg
     * @param init_type
     */
    State( ClusterGraph& clg, StateInitializationType init_type = RANDOM_MATCHING );

    /**
     * Resets all arrays - fill with 'initializer' values.
     * Does not clear
     */
    void clearState();

    /**
     * Creates state for given partition [part] of clg->V.
     * @param part
     */
    void applyPartition( VI & part );

    /**
     * Applies swap for given nodes - moves to_swap[i].first from its cluster to cluster to_swap[i].second
     *
     * CAUTION!
     * Changes the whole structure of state, clusters and all arrays. It works by creating completely new state, then
     * merges clusters properly.
     *
     * If there are some sets of nodes that should be moved to separate new clusters, then corresponding nodes should
     * have marked to_swap[i].second value at least getIdOfEmptyCluster()
     */
    void applySwap( VPII & to_swap );

    /**
     * Creates all necessary arrays. Initially all nodes of cluster graph clg are in separate clusters.
     * Initial setup may be a bit memory consuming, as each cluster keeps induced subgraph.
     */
    void initializeStateData( StateInitializationType init_type );

    /**
     * Initializes state using sparseGraphTrimmer.
     */
    void sparseGraphTrimming();

    /**
     * For each pair (a,b) in [part], merges clusters a and b.
     * If there are pairs (a,b) and (b,c) then a,b,c are merged together.
     * Merging clusters means creating a cluster that contains nodes of all merged clusters.
     *
     * An empty cluster is added to the end of [clusters].
     *
     * ************ !!!!!!!!!!!!!1 *************
     * CAUTION!!! [cl_neigh_graph] structure is invalidated !!!
     * ************ !!!!!!!!!!!!!1 *************
     *
     * Works in time O( N +  \sum deg(v) ), where sum runs over all nodes in merged clusters
     * @param c1
     * @param c2
     */
    void mergeClusters( VPII & part );
    void mergeClusters( VVI to_merge );

    /**
     * Creates cl_neigh_graph structure (see [cl_neigh_graph] for description).
     */
    void createClNeighGraph();

    int getIdOfEmptyCluster(){ assert( !clusters.empty() ); return clusters.size()-1; }

    friend ostream& operator<<(ostream& str, State& st);

    /**
     * Finds the number of edge modifications needed to transform cluster graph [clg] to the state represented by
     * [clusters]
     */
    int calculateResultForState();


    //*********************//*********************//*********************//********************* fields

    /**
     * size of the graph
     */
    int N;

    /**
     * inCl[i] is the id of the cluster where node i belongs
     */
    VI inCl;

    /**
     * degInCl[i] is the degree of node i in cluster inCl[i]. By degree we mean sum of weights of edges.
     */
    VI degInCl;

    /**
     * idInCl[i] is the id of node i in graph induced by cluster i.
     * In other words indInCl[i] == clusters[ inCl[i] ]->g.perm[i].
     * This is much faster, as we do not use unordered_map, but just an array
     */
    VI idInCl;

    /**
     * Vector containing clusters of the current state
     */
    vector<Cluster> clusters;

    /**
     * Cluster graph created for
     */
    ClusterGraph *clg;

    /**
     * cl_neigh_graph is a structure, that for node i keeps pairs (cl,w) where cl is a cluster neighbor of node i, and
     * w is the sum of weights of edges with one end in i and the other in cluster cl.
     *
     * Cluster inCl[i] is NOT included in cl_neigh_graph[i] list.
     * For each i pairs (cl,w) should be sorted by non-increasing w
     *
     * This is used to quickly find cluster neighbors of given node.
     */
    VVPII cl_neigh_graph;

    /**
     * Some hashes
     */
    VLL hashes;


    VB was, was2, was3; // 'Global' helper arrays of size N
    VI marker; // helper int array
//    VVI marker2; // helper VVI array
};

#endif //CESWAT_STATE_H
