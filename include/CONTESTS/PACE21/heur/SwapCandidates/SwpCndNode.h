//
// Created by sylwester on 3/20/21.
//

#ifndef ALGORITHMSPROJECT_SWPCNDNODE_H
#define ALGORITHMSPROJECT_SWPCNDNODE_H

#include "SwapCandidate.h"

/**
 * Swap Candidate - moving triangles.
 */
class SwpCndNode : public SwapCandidate{
public:
    /**
     * Creates a swap candidate for moving node [v] to cluster [trg_cl].
     */
    SwpCndNode( int swpval, int v, int trg_cl );

    VI getAffectedClusters(State & st) override{ return {st.inCl[node], move_node_to}; }

    VPII getNodesToSwap() override{ return {{node,move_node_to}}; }
    LL swpVal() override{ return swap_value; }

    int node, move_node_to, swap_value;

    VI getNodes(){ return {node}; }
    VI getMoveNodesTo(){ return {move_node_to}; }

    friend ostream& operator<<(ostream& str, SwpCndNode& cnd);
};

/**
 * Considers swap candidates by moving single nodes from/to clusters. There are (wll be) a few possibilities.
 */
class SwpCndNodeCreator : public SwapCandidateCreatorAdapter{
public:
    SwpCndNodeCreator( State & s );

    vector<SwapCandidate *> createSwapCandidatesRaw() override;

    vector<SwpCndNode> createSwapCandidates();

    /**
     * Creates swap candidates for given cluster. In the resulting vector there are only single-node candidates, that
     * can be moved from one cluster to another (possible empty).
     *
     * If [keep_only_nonpositive_candidates] is set, then only found candidates with swap_vale <=0 will be returned.
     *
     * If [keep_only_best_cluster_to_move_to] i set, then for each node v in cluster [cl] only the best option will
     * be returned. If there is a tie, the first one found will be returned. If that option has positive swap value and
     * [keep_only_nonpositive_candidates] is set, then for this node v nothing will be returned.
     */
    vector<SwpCndNode> create_MoveTo_SwapCandidatesForCluster( Cluster& cl );


};

#endif //ALGORITHMSPROJECT_SWPCNDNODE_H
