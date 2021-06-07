//
// Created by sylwester on 4/22/21.
//

#ifndef ALGORITHMSPROJECT_NODEEDGEGREEDYNOMAP_H
#define ALGORITHMSPROJECT_NODEEDGEGREEDYNOMAP_H

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
class NodeEdgeGreedyNomap : public NEG{
public:

    NodeEdgeGreedyNomap( State & st );

    virtual void initializeForState(State & st) override;
    virtual void improve() override;

//private:

    template<class _T> void localShuffle(_T & v);


    //****************************************  PERTURBATIONS   ************************

    /**
     * edges_to_cluster[v] is a map such that edges_to_cluster[v][cl] is the number of edges between node v and cluster
     * cl.
     */

    VVPII edges_to_cluster;

    /**
     * CAUTION! If use_sort is set to false, then [findEdgesToCluster] may not work properly if it uses binary search!!
     * use sort can be set to false to quicker create edges_to_cluster[v], if we are only interested in iterating over it
     */
    virtual void createEdgesToCluster(int v, bool use_sort = true);
    virtual int findEdgesToCluster( int v, int cl );
    VI helper_etocl;
    VB helper_was_etocl;
//    VI degInCl;
    bool use_binary_search_etocl = true;
    const int MAGIC_SORT_THR = 15;


    virtual tuple<int,int,int> getBestNodeMoveForRange( VI & perm, int a, int b ) override;
    vector<tuple<int,int,int>> best_node_move_results;

    virtual  SwapCandidateAdapter getBestTriangleAll(int v) override;
    VB helper_was4; // used in getBestTriangleAll

    virtual void moveNodeTo(int v, int to) override;

    virtual tuple<PII, PII, int> getBestEdgeMoveForRange(VI & perm, int a, int b ) override;
//    vector<tuple< PII, PII,int >> best_edge_move_results;

    virtual int swapValueForNode(int v, int trg_cl, int edges_clv, int tot_clv_possible_edges, int edges_trg) override;
    virtual VPII getEdgesToCluster(int v) override;
    virtual void createClusterNodes() override;
    virtual VI getClusterNodes(int c) override;
    virtual int maxClusterId() override;
};

#endif //ALGORITHMSPROJECT_NODEEDGEGREEDYNOMAP_H
