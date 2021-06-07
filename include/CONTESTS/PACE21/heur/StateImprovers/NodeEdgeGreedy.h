//
// Created by sylwester on 4/22/21.
//

#ifndef ALGORITHMSPROJECT_NODEEDGEGREEDY_H
#define ALGORITHMSPROJECT_NODEEDGEGREEDY_H

#include <CONTESTS/PACE21/heur/Config.h>
#include "../State.h"
#include "../SwapCandidates/SwapCandidate.h"
#include "NEG.h"

/**
 * Algorithm works in iterations.
 * In each iteration a random permutation of nodes is chosen. Then, in that order, nodes are tried to be moved to some
 * other cluster. If, at any point, a node with NONPOSITIVE swap value is found, then it is instantaneously moved.
 *
 * The same goes for checking edges.
 *
 * If during whole iteration no improvement was found, algorithm stops.
 */
class NodeEdgeGreedy : public NEG{
public:

    NodeEdgeGreedy( State & st );

    virtual void initializeForState(State & st) override;
    virtual void improve() override;

//private:
    template<class _T> void localShuffle(_T & v);

    int countNonemptyClusters() override;
    virtual void resizeStructuresForEmptyCluster( int empty_cl ) override;

    vector< unordered_set<int> > cluster_nodes;


    int countClusterNumerationGaps() override;


    //****************************************  PERTURBATIONS   ************************

    /**
     * edges_to_cluster[v] is a map such that edges_to_cluster[v][cl] is the number of edges between node v and cluster
     * cl.
     */
    vector< unordered_map<int,int,quick_hash> > edges_to_cluster;
//    vector< unordered_map<int,int> > edges_to_cluster;
//    vector< unordered_map<int,int,quick_hash> > edges_to_cluster;
//    vector< map<int,int> > edges_to_cluster;

    /**
     * Checks all nodes from perm[a], perm[a+1], ..., perm[b], then returns best move in the form
     * (node, cluster, swap_value)
     */
    virtual tuple<int,int,int> getBestNodeMoveForRange( VI & perm, int a, int b ) override;


    virtual SwapCandidateAdapter getBestTriangleAll(int v) override;

    virtual void moveNodeTo(int v, int to) override;

    virtual tuple<PII, PII, int> getBestEdgeMoveForRange( VI & perm, int a, int b ) override;

    void createEdgesToCluster(int v, bool use_sort) override;

    int findEdgesToCluster(int v, int cl) override;

    VPII getEdgesToCluster(int v) override;

    void createClusterNodes() override;

    VI getClusterNodes(int c) override;

    int maxClusterId() override;
};

#endif //ALGORITHMSPROJECT_NODEEDGEGREEDY_H
