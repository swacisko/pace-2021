//
// Created by sylwester on 4/8/21.
//

#ifndef ALGORITHMSPROJECT_COMPONENTEXPANSIONATTRACTION_H
#define ALGORITHMSPROJECT_COMPONENTEXPANSIONATTRACTION_H

#include "../ExpansionOrder.h"
#include "SwapCandidate.h"

class ExpansionOrderAttraction{
public:
    ExpansionOrderAttraction( Cluster& cl, VI v ){ this->cl = &cl; ord = v; }
    VI ord; // just the order
    Cluster* cl; // pointer to the cluster to which nodes are moved

    static vector<ExpansionOrder> induceOrders(State & st, ExpansionOrderAttraction & eo );
};

class SwpCndEOAttraction : public SwapCandidate {
public:
    // we do not need set hash value, because ComponentExpansionAttraction will not be called more that once for each
    // cluster. It is still possible, however, that induced orders will 'coincide' with other ExpansionOrders created
    // using some other SwpCndCreators
    SwpCndEOAttraction( ExpansionOrderAttraction & eo, int k, int swpval ){
        this->k = k; this->eo = &eo; swap_value = swpval;
    }

    VPII getNodesToSwap() override{
        VI nodes = getNodes(); VI move_to = getMoveNodesTo();
        return StandardUtils::zip(nodes, move_to);
    }

    VI getNodes(){ return StandardUtils::getSubarray(eo->ord,0,k); }
    VI getMoveNodesTo(){ return VI(k+1, eo->cl->id); } // all nodes are moved to cluster eo->cl
    virtual LL swpVal(){ return swap_value; }

    friend ostream& operator<<( ostream& str, SwpCndEOAttraction& cnd );

//private:
    ExpansionOrderAttraction* eo;
    int k; // number such that in this swap candidate there are nodes eo.ord[i] for i <= k
    int swap_value;

};

/**
 * This component expansion works as follows:
 * Consider a cluster S and a set {C_1,C_2, .., C_p} of neighbors of cluster S. We iteratively select from  C = \cup C_i
 * a node u that minimizes swap value swp(u,S), then move that node do S and update all necessary data.
 *
 * Let also denote C = \sum |C_i|.
 * In order to select proper node u, we keep a list of different node sizes
 * (there are at most O(sqrt(C)) different sizes). For each such node size we select the best node u of that size,
 * that minimizes swap value swp(u,S). We do that in O(1) time - the 'core values' are stored in a set in order to update
 * them quickly.
 *
 *
 * swap value for node u is the following:
 * nw_u * ( sumNWinS - cluster_weights[c] ) + 2*edges_in_cl[u] - 2*eToS[u] + (nw_u)^2   <---- swap value swp(u,S)
 * -nw_u * cluster_weights[c] + 2*edges_in_cl[u] - 2*eToS[u] + (nw_u)^2    <----- core value ('free value')
 *
 *
 * Hence, in each iteration a node u can be chosen in (most pessimistically) time O( sqrt(C) ). This gives the
 * total runtime for given cluster S: O( C * sqrt(C) + E(C) * logC ), where E(C) is the number of edges in the graph
 * induced by S \cup \bigcup C_i.
 *
 * It is possible (perhaps) to keep a binary search tree (just a set) for swap values calculated for best representatives of
 * all possible node_weights, to replace sqrt(C) factor by logC factor.
 *
 */
class ComponentExpansionAttraction{
public:

    ComponentExpansionAttraction( State & st );

    /**
     * Creates swap candidates for all clusters in state, then aggregates and return them.
     * @return
     */
    pair<vector<ExpansionOrderAttraction*>, vector<SwpCndEOAttraction>> createSwapCandidates();

    /**
     * Creates swap candidates.
     *
     * For more details about how those candidates are created, see the description of the class.
     *
     * @return a pair containing pointer to the dynamically allocated ExpansionOrderAttraction and a vector of swap
     * candidates. Those swap candidates store pointers to the EOA. Given EOA need to be deleted after when no longer
     * used, in order to avoid memory leak.
     */
    pair<ExpansionOrderAttraction*, vector<SwpCndEOAttraction> > createSwapCandidates(Cluster& cl);



    void initialize();

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
     * eToS[v] is the SUM OF WEIGHTS OF EDGES from node u to set S
     */
    VI eToS;

    /**
     * edges_in_cl[u] is the number of edges between node u and the rest of its cluster
     */
    VI edges_in_cl;

    VB inS; // boolean marker for S
    int sumNWinS; // sum of node weights of nodes in S

    /**
     * clusters that are neighbors to cluster [cl].
     */
    VI cl_neigh;

    /**
     * 'set' of different cluster weights among all neighboring clusters of S
     */
    unordered_map<int,int> node_neigh_weights;

    // ************************************************************

    /**
     * cluster_cores[x] is the set that contains pairs (free_value,v) for all nodes v \notin S with nw_v = x.
     *
     * swap value for node u is the following:
     * nw_u * ( sumNWinS - cluster_weights[c] ) + 2*edges_in_cl[u] - 2*eToS[u] + (nw_u)^2
     *
     * For a fixed nw_u, by keeping for only 'free value' (core value)
     * -nw_u * cluster_weights[c] + 2*edges_in_cl[u] - 2*eToS[u] + (nw_u)^2
     *
     * in sorted order, we are able to quickly node that minimizes its swap value (among all nodes with given nw_u).
     */
    vector< set<PII> > cluster_cores;

    /**
     * Calculates and returns 'free value' for given node d.
     */
    int getFreeValue( int d );

    // ************************************************************


    /**
     * Removes/inserts all neighbors of node [node_to_move] from [cluster_cores]. It is necessary to remove all
     * neighbors of moved node, before moving the node, because otherwise we could not calculate the 'free value' used
     * as key to store nodes in sets.
     * It also removes/insert all nodes that are in the same cluster as [node_to_move]. Since the density of each
     * cluster should be at least 0.5, this is not a worry.
     *
     * @param insert if true, then entry if node u will be inserted, otherwise it will be erased.
     */
    void modifyNodeEntryInClusterCores( int node_to_move, bool insert );

    /**
     * Finds node u that should be moved to S.
     */
    int getNodeToMove();

    // ************************************************************

    /**
     * Set to which nodes are attracted.
     */
    unordered_set<int> S;

    /**
     * cluster_weights[i] is the weight of cluster st->clusters[i]. This value is changing during node attraction
     */
    VI cluster_weights;

    VB was, was2, was3; // helper arrays

    /**
     * Moves node v from its cluster to set S
     */
    void moveNodeToS(int v);

    /**
     * If true, then only those swap candidates (  ord[0], ord[1], ... , ord[k] ) for which swpval <= 0 will be kept.
     */
    bool keep_only_nonpositive_candidates = true;
};


#endif //ALGORITHMSPROJECT_COMPONENTEXPANSIONATTRACTION_H
