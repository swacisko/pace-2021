//
// Created by sylwester on 3/20/21.
//

#ifndef ALGORITHMSPROJECT_SWPCNDEDGE_H
#define ALGORITHMSPROJECT_SWPCNDEDGE_H

#include "SwapCandidate.h"

/**
 * Swap Candidate - moving edges.
 */
class SwpCndEdge : public SwapCandidate{
public:
    SwpCndEdge() = default; // default
    /**
     * Creates a swap candidate for moving edge (u,v) to cluster [trg_cl].
     */
    SwpCndEdge( int swpval, int u, int v, int trg_cl );


    VI getAffectedClusters(State & st) override{
        unordered_set<int> zb;
        zb.insert( st.inCl[u] );
        zb.insert( st.inCl[v] );
        zb.insert( move_node_to );
        return VI(ALL(zb));
    }

    VPII getNodesToSwap() override{ return {{u,move_node_to}, {v,move_node_to}}; }
    LL swpVal() override{ return swap_value; }

    virtual VI getNodes(){ return {u,v}; }
    virtual VI getMoveNodesTo(){ return {move_node_to, move_node_to}; }

    int u,v, move_node_to, swap_value;

    friend ostream& operator<<(ostream& str, SwpCndEdge& cnd);

};

class SwpCndEdgeCreator : public SwapCandidateCreatorAdapter{
public:
    SwpCndEdgeCreator(State &s);

    vector<SwapCandidate *> createSwapCandidatesRaw() override;

    vector<SwpCndEdge> createSwapCandidates();

    /**
     * Considers all edges as swap candidates. For each edge e in [cl] it checks, whether that edge (both nodes) can
     * be moved to other cluster.
     *
     * If [only_common_neighbors] is set, then both ends of considered edge must have an endpoint in a cluster to which
     * they are supposed to be moved. Otherwise it suffices that only one node has. If [only_common_neighbors] is true,
     * then running time is O( E * sqrt(E) ), otherwise it might get pessimistically O(E * C), where C is the total number of
     * clusters.
     *
     * If [keep_only_nonpositive_candidates] is set, then only found candidates with swap_vale <=0 will be returned.
     *
     * If [keep_only_best_cluster_to_move_to] i set, then for each edge in cluster [cl] only the best option will
     * be returned. If there is a tie, the first one found will be returned. If that option has positive swap value and
     * [keep_only_nonpositive_candidates] is set, then for this edge nothing will be returned.
     */
    vector<SwpCndEdge> create_MoveTo_SwapCandidatesForCluster( Cluster& cl, const bool only_common_neighbors = true );

    /**
     * Checks whether there exists an edge (u,v) with u and v in different clusters, such that that both its ends can
     * be moved to other cluster. Ther rest is similar as in [create_MoveTo_SwapCandidatesForCluster] function.
     *
     * This function is very similar to [create_MoveTo_SwapCandidatesForCluster], changes are done in just a few places.
     */
    vector<SwpCndEdge> create_MoveTo_SwapCandidates_DifferentClusters( const bool only_common_neighbors = true );

};

#endif //ALGORITHMSPROJECT_SWPCNDEDGE_H
