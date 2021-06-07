//
// Created by sylwester on 3/20/21.
//

#ifndef ALGORITHMSPROJECT_SWPCNDEO_H
#define ALGORITHMSPROJECT_SWPCNDEO_H

#include <CONTESTS/PACE21/heur/ExpansionOrder.h>
#include <CONTESTS/PACE21/heur/ConvexHullTrickDynamic.h>
#include "SwapCandidate.h"

/**
 * Swap Candidate for Expansion Order
 */
class SwpCndEO : public SwapCandidate{
public:
    SwpCndEO( ExpansionOrder& eo, int k, int move_to, int swap_value, LL hash );
    VPII getNodesToSwap() override;
    LL swpVal() override;

    VI getAffectedClusters(State& st) override;

    friend ostream& operator<<(ostream& str, SwpCndEO& cnd);

    ExpansionOrder* eo;
    int k; // number such that in this swap candidate there are nodes eo.ord[i] for i <= k
    int move_to; // id of the cluster this swap candidate should be moved to
    int swap_value;

    /**
     * set_hash is the hash of the set represented by this object. It is used to keep only one copy, if many same sets
     * would be created for different expansion orders.
     */
    LL set_hash;
};

class SwpCndEOCreator : public SwapCandidateCreatorAdapter{
public:
    /**
     * Constructing the object works in time O(N), where N is the size of the state cluster graph.
     * @param s
     */
    SwpCndEOCreator(State & s);

    void createHashes();

    vector<SwapCandidate *> createSwapCandidatesRaw() override;

    /**
     * Creates swap candidates for given expansion order [eo].
     * Expansion order is created for cluster eo.cl. We need also structure of the cluster graph.
     *
     * For each of eo.ord.size() sets X_i = {v_0,..., v_i} of expansion order eo, we find the best cluster to move the
     * set X_i to.
     *
     * If [keep_only_nonpositive_candidates] is set, then all setx X_i for which swpVal() is negative will be discarded.
     *
     * If [keep_only_best_cluster_to_move_to] is true, then we will use ConvexHullOptimization trick to calculate only
     * the best cluster-candidate-to-move-to. In that case the complexity for given expansion order is pessimistically
     * O( E*logC ), where E is the number of edges incident to at least one node in [eo] and C is the number of
     * 'edges; in state->cluster_neigh_graph with an end in a node from [eo]
     *
     * If [keep_only_best_cluster_to_move_to] is false, then all neighboring clusters will be checked. In this case the
     * complexity might grow to O( E*C ). Depending on test characteristics, however, it may work ok or not.
     */
    vector<SwpCndEO> createSwapCandidates(ExpansionOrder & eo);

    /**
     * For each expansion order creates SwpCndEO using createSwapCandidates(ExpansionOrder & eo). Returns all found
     * candidates, without repetitions (because there may be many same candidates originating from different orders).
     *
     * It is necessary to store expansion orders as long as SwpCndEO are used, since they store a pointer to the
     * order from which they originate.
     *
     * Considers only those expansion orders that are in subarray exp_orders[beg:end), excluding [end].
     */
    vector<SwpCndEO> createSwapCandidates( vector<ExpansionOrder*> & exp_orders, int beg = 0, int end = -1 );


    /**
     * For each cluster cl creates an expansion order for each of fun(cl) sets. Then, for each
     * expansion order, creates candidates using createSwapCandidates(ExpansionOrder & eo) function. Returns the union
     * of all swap candidates.
     *
     * Creating expansion orders is done separately, since they need to be in scope to create swap candidates.
     */
    vector<ExpansionOrder> createExpansionOrders( function<VVI(Cluster*)> fun );

    /**
     * Number of edges between set X and Y (Y is just the rest of the cluster C \ X ).
     */
    int cut_value;

    /**
     * sum of weights of nodes in the set X
     */
    int sum_nw;

    /**
     * ConvexHull for quick acces to best clusters
     */
    ConvexHullTrickDynamic hull;

    /**
     * min_b_coef_for_cluster_size[i] is the minimal values of the b coefficient among all lines added for cluster
     * with size i. This way, we will add new line to [hull] only if
     * -2 * edges_to_cluster[cl] < min_b_coef_for_cluster_size[cl]. Otherwise adding the line will have no effect and
     * can be omitted.
     */
    VI min_b_coef_for_cluster_size;

    /**
     * Vector of random variables used to create set hashes
     */
    VLL hashes;

//private:




};

#endif //ALGORITHMSPROJECT_SWPCNDEO_H
