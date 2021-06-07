//
// Created by sylwester on 4/8/21.
//

#ifndef ALGORITHMSPROJECT_COMPONENTEXPANSIONREPULSION_H
#define ALGORITHMSPROJECT_COMPONENTEXPANSIONREPULSION_H


#include <CONTESTS/PACE21/heur/ExpansionOrder.h>
#include "datastructures/Heap.h"
#include "SwapCandidate.h"

class ExpansionOrderRepulsion{
public:
    /**
     * When 'expanding - attracting' given cluster [cl], node v[i] should be moved to cluster move_to[i]
     */
    ExpansionOrderRepulsion( Cluster& cl, VI v, VI to );
    VI ord; // just the order
    VI move_to; // move_to[i] is the id of the cluster to which node ord[i] should be moved
    Cluster* cl; // pointer to the cluster for which the EO was created.

    /**
     * From this expansion order creates expansion orders for all neighboring cluster of cluster [cl].
     * Those created expansion orders are 'induced' by nodes that are moved to those clusters.
     * (e.g. if ord = {1,5,2,7,3} and move_to = {1,2,2,1,8}, then expansion order induced for cluster 1 will be {1,7}).
     * All induced expansion orders are created for cluster [eo.cl]
     *
     * Since nodes in ExpansionOrderRepulsion can be moved to many new empty clusters, it is at this point impossible to
     * create those clusters. That is why, for each new empty cluster a new expansion order is created, but the the
     * cluster (reference) is st.clusters[ st.getIdOfEmptyCluster() ]. This way for an empty clusters many expansion
     * orders may be returned.
     */
    static vector<ExpansionOrder> induceOrders(State & st, ExpansionOrderRepulsion & eo );
};

/**
 * For some more details see SwapCandidate and SwapCandidateAdapter
 */
class SwpCndEORepulsion : public SwapCandidate{
public:
    SwpCndEORepulsion( ExpansionOrderRepulsion * eo, int k, int swap_value, LL hash );

    VPII getNodesToSwap() override;

    VI getNodes(){ return StandardUtils::getSubarray(eo->ord,0,k); }
    VI getMoveNodesTo(){ return StandardUtils::getSubarray(eo->move_to,0,k); }
    virtual LL swpVal(){ return swap_value; }

    friend ostream& operator<<( ostream& str, SwpCndEORepulsion& cnd );

//private:
    ExpansionOrderRepulsion* eo;
    int k; // number such that in this swap candidate there are nodes eo.ord[i] for i <= k
    int swap_value;
    LL set_hash;
};

/**
 * Given a cluster S, it iteratively selects node in S, then moves this node to some other cluster to which it has
 * smallest swap value. Repeats until S is not empty.
 *
 * The swap value for node u and cluster c is
 * nw_u * (nw_c - sumNWinS) + nw_u^2 - 2*edges_to_cluster[u][c] + 2*eInS[u].
 * Since  -nw_u * sumNWinS + nw_u^2 + 2*eInS[u] is independent of c, for fixed u we can select cluster c that minimizes
 * nw_u * nw_c - 2*edges_to_cluster[u][c].    <--- core value
 * We keep those core values for each node u to each of the clusters, where given node u has neighbors.
 * This way, for given u, we can in O(1) access the cluster to which it is best to move it.
 * Whenever node u is moved from S to some cluster c, then core_value[v][c] for all remaining v \in S nodes are updated.
 *
 * Nodes in S will be updated at most S.size() times.
 * Works in time O( S^2 * logD + D ), where D is the number of clusters that are neighbors to cluster S.
 * However, the constant is probably not small... :(
 *
 * CAUTION!! If there is a node u that is connected to every other node in the graph, then D will always be the number
 * of clusters, st->clusters.size(). If so, then calling the procedure for every cluster in the graph may take time
 * E * logN + N^2  !!!   <--- CAUTION!!
 */
class ComponentExpansionRepulsion{
public:

    ComponentExpansionRepulsion( State & st );

    /**
     * Creates swap candidates for all clusters in state, then aggregates and return them.
     * @return
     */
    pair<vector<ExpansionOrderRepulsion*>, vector<SwpCndEORepulsion>> createSwapCandidates();

    /**
     * Creates expansion order and swap candidates for given cluster.
     * Works in tim
     *
     * **********   !!! **********
     * CAUTION!!! Created SwpCndEORepulsion may contain many nodes that are to be moved to empty clusters!! Some nodes
     * may be moved to the same newly created clusters.
     * IDS OF CLUSTERS TO WHICH NODES ARE MOVED MAY BE GREATER THAN st->getIdOfEmptyCluster()
     * **********   !!! **********
     *
     * IMPORTANT!! If the whole set S should be moved to a single empty cluster, then we should not return that as a
     * swap candidate
     *
     * @return pair containing expansion order and a vector of candidates.      CAUTION!!!!!! Order needs to be deleted
     * after using to avoid memory leak.
     */
    pair<ExpansionOrderRepulsion*, vector<SwpCndEORepulsion> > createSwapCandidates(Cluster& cl);

//protected:
    int min_cluster_size = 1;

    /**
     * Pointer to the cluster graph in [state].
     * THIS IS NOT INDUCED_CLUSTER_GRAPH, AS IT IS IN USUAL COMPONENT EXPANSION (ComponentExpansion.h)
     */
    ClusterGraph * clg;

    State* st; // pointer to the state
    Cluster* cl; // pointer to the current cluster from which we move nodes to other clusters

    VVPII * V; // V = &clg.V

    int N; // clg.V.size()

    /**
     * eInS[v] is the SUM OF WEIGHTS OF EDGES of node v in set S.
     * When a node u is removed from [cl] to cluster C, then eInS[v] is decreased for all v \in (N(u) \ (C+S))
     */
    VI eInS;

    VB inS; // boolean marker for S
    int sumNWinS; // sum of node weights of nodes in S

    /**
     * clusters that are neighbors to cluster [cl].
     */
    VI cl_neigh;

    // ************************************************************

    struct ClusterSet{
        ClusterSet( ComponentExpansionRepulsion* cer, int v, int c ){ this->cer = cer; this->v = v; this->c = c; }
        bool operator<( const ClusterSet& oth ) const{
            assert( v == oth.v );
            int core1 =  cer->swpval_core_for_node[v][c];
            int core2 = cer->swpval_core_for_node[v][oth.c];
            if( core1 != core2 ) return core1 < core2;
            else if( cer->cluster_weights[c] != cer->cluster_weights[oth.c] ){
                // we return smaller set as the better candidate
                return cer->cluster_weights[c] < cer->cluster_weights[oth.c] ;
            }
            else return c < oth.c; // we return cluster with larger id
        }
        bool operator==( const ClusterSet& oth ) const{ return cer == oth.cer && v == oth.v && c == oth.c; }
        ComponentExpansionRepulsion* cer;
        int v,c;
        friend ostream& operator<<(ostream& oth, ClusterSet& cs);
    };

    /**
     * cluster_sets[v] is a set of pairs (swpval,i) such that moving node v to cluster i results in swap value swpval.
     * This is used to quickly access best cluster to which node v should be moved
     */
    vector< set<ClusterSet> > cluster_sets;

    /**
     * swpval_for_node[v][c] is the core for swap value of moving node v to cluster c. We cannot store exact swap
     * values because each time we remove node from S we would have to update ALL nodes remaining in S for ALL
     * neighboring cluster. This would lead to O( S^2*D ) complexity, where D is the maximum number of node in S.
     * By keeping only the core, we do not know the exact value of swap value, but we have access in O(1) time to
     * cluster to which given node from S has minimal swap value. Knowing that cluster we can in O(1) calculate that
     * swap value. This lead to O(S) complexity instead of O(S*D) for finding best choice for move. This comes, however,
     * at a cost of slightly increased constant and O(logD) update time for each of up to S^2 updates.
     * Thus we get complexity O( S^2 * logD ).
     */
    vector< unordered_map<int,int> > swpval_core_for_node;

    /**
     * Updates swap value core for moving node u \in S to cluster c.
     */
    void updateSwpValCoreForNode(int u, int c );

    /**
     * Finds pair (u,c) - node u \in S and cluster c such that u should be moved to c.
     */
    PII getNodeToMove();

    // ************************************************************

    /**
     * Expanded set.
     */
    unordered_set<int> S;

    /**
     * cluster_weights[i] is the weight of cluster st->clusters[i]. This value is changing during 'repulsion
     */
    VI cluster_weights;

    /**
     * For v \in cl.g.nodes we have edges_to_cluster[v][c] is the sum of weights of edges of node v to cluster c
     */
    vector< unordered_map<int,int> > edges_to_cluster;

    VB was, was2, was3; // helper arrays

    /**
     * Moves node v from S to cluster c
     */
    void moveNodeTo(int v, int c);

    /**
     * Calculates and returns swap value for moving node u to cluster c
     */
    int calculateSwpVal(int u, int c );


    /**
     * If true, then only those swap candidates (  ord[0], ord[1], ... , ord[k] ) for which swpval <= 0 will be kept.
     */
    bool keep_only_nonpositive_candidates = true;

};

#endif //ALGORITHMSPROJECT_COMPONENTEXPANSIONREPULSION_H
