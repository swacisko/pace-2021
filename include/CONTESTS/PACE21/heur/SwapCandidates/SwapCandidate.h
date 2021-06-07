//
// Created by sylwester on 3/8/21.
//

#ifndef ALGORITHMSPROJECT_SWAPCANDIDATE_H
#define ALGORITHMSPROJECT_SWAPCANDIDATE_H

#include <CONTESTS/PACE21/heur/State.h>
#include <utils/StandardUtils.h>
#include <Constants.h>
#include "Makros.h"
#include "CollectionOperators.h"
#include "CONTESTS/PACE21/heur/Cluster.h"

class SwapCandidate{
public:
//protected:

    /**
     * Type for number of edges - in large grapgs, starting from a single cluster, number of edges may need to be
     * of type long long.
     */
    typedef int int_type;

    /**
     * Default virtual destructor.
     * Though it is purely virtual, it has definition in SwapCandidate.cpp file, to avoid compilation problems.
     */
    virtual ~SwapCandidate() = default;

    /**
     * res[i] is a pair(a,id) such that [a] should be moved to cluster with given [id].
     * Ids of nodes and clusters should correspond to those of the current state (of current ClusterGraph).
     * @return
     */
    virtual VPII getNodesToSwap() = 0;

    /**
     * @return swap value, that is the difference in scores (score after change minus current score)
     */
    virtual LL swpVal() = 0;


    /**
     * Finds and returns all affected clusters by the swap, given a state [st]. It may be impossible to access all
     * affected clusters without the state. By default, uses [getAffectedClusters()] function
     */
    virtual VI getAffectedClusters( State & st ){
        unordered_set<int> zb;
        auto p = StandardUtils::unzip( getNodesToSwap() );
        zb += p.second;
        for( int d : p.first ) zb.insert( st.inCl[d] );
        return VI(ALL(zb));
    }

    /**
     * Creates a vector of pointers to cluster that will be affected by the swap.
     * @param clusters
     * @param cl cl[i] should denote id of the cluster where node i belongs
     * @return
     */
    static vector<Cluster*> getAffectedClusters( VPII & toSwap, vector<Cluster> & clusters, VI& cl );


public:
    static void test();

};

class SwapCandidateAdapter : public SwapCandidate{
public:

    /*int_type ein(){ return Ein; }
    int_type eout(){ return Eout; }*/

    virtual ~SwapCandidateAdapter(){}

    /**
     * Ein - number of edges between nodes from swap candidate. It includes all edges that are 'internal' nodes of
     * ClusterGraph's nodes (there are  (nw() choose 2)  such edges for each node in the cluster ).
     * Eout - total number of edges going out of the swap candidate (that is sum of weights)
     */
//    int_type Ein, Eout;

    virtual LL swpVal(){ return swap_value; }
    virtual VPII getNodesToSwap(){ return StandardUtils::zip(getNodes(), getMoveNodesTo()); }

    /**
     * nodes - Nodes that are in the swap candidate.
     * move_node_to - nodes[i] should be moved to cluster move_node_to[i]
     */
    VI nodes, move_node_to;
    virtual VI getNodes(){ return nodes; }
    virtual VI getMoveNodesTo(){ return move_node_to; }

    /**
     * Swap value - see [swpVal] function.
     */
    int swap_value;
};

class SwapCandidateCreator{
    /**
     * Create a set of swap candidates and return them as pointers.
     * CAUTION!! They will need to be deleted after using!!
     */
    virtual vector<SwapCandidate*> createSwapCandidatesRaw() = 0;
};

class SwapCandidateCreatorAdapter : public SwapCandidateCreator{
public:
    SwapCandidateCreatorAdapter( State & s );

    /**
     * Just overriding
     */
    virtual vector<SwapCandidate*> createSwapCandidatesRaw() override { return vector<SwapCandidate*>(); }

    /**
     * clEout[i] is the number of edges between this swap candidate and cluster i. Onlt record with positive values are
     * stored.
     */
//    unordered_map<int,int> clEout;

    State * state; // pointer to the state
    ClusterGraph* clg; // pointer to the cluster graph of state [state]

    VI edges_to_cluster; // helper array of size(cluster.size()

    VB was, was2, was3; // helper boolean arrays

    /**
     * For each neighbor p of node d, if p and d are ain different clusters (if p is not in [in_cl_d]), then:
     * - if in_cl_p was not visited yet, in_cl_p is added to [neigh].
     * - number of edges edges_to_cluster[in_cl_p] is increased byt weight of edge (d,p).
     *
     * CAUTION! If state->cl_neigh_graph is empty(), then state->createClNeighGraph() is called.
     */
    void updateEdgesInNeighboringClustersForNode( int d, int in_cl_d, VI& neigh );

    /**
     * IMPORTANT TO REMEMBER WHAT swp_val1 is !!
     * * * * * * * * * * *  IMPORTANT! Number of edges is accessed by edges_to_cluster[swp_trg_cl]
     *
     * @param swp_val1 this is total number of possible edges between swap candidate and the rest of its cluster, decreased by the
     * number of present edges.
     * @param swpcnd_w weight of the swap candidate
     * @param swp_trg_cl if of the cluster to which the swap candidate shall be moved
     * @param deg_in_cl degree of moved node in its cluster (just state->degInCl[d])
     * @return final swap value of moving swap candidate to target cluster
     */
    int getSwpValForMove( int swp_val1, int swpcnd_w, int deg_in_cl, int swp_trg_cl );

    /**
     * For given swap candidate, it checks all clusters c in neigh_cl.
     * If the requirements are met (e.g. nonpositive and only_best), then a candidate is appended to [res]. A candidate
     * object is created by calling   constructor(final_swpval, trg_cl).
     *
     * @param swpval1 traditional - number (sum of weights) of all possible edges between considered swpcnd and the rest of its cluster,
     * decreased the number (sum of weights) of present edges there
     * @param swpcnd_w weight of the considered swap candidate
     * @param swpcnd_deg_in_cl degree in cluster - total weight of edges between this swpcnd and the rest of its cluster
     * @param neigh_cl clusters that are connected to the cluster of considered swap candidate
     * @param res vector where swap candidates are to be added
     * @param constructor creates swap candidate objects
     */
    template<class _T, class _cr>
    void createAllCandidatesForSwpCnd( int swpval1, int swpcnd_w, int swpcnd_deg_in_cl, VI & neigh_cl, vector<_T> & res, _cr constructor ){
        const bool debug = false;

        if(debug){
            DEBUG(swpval1); DEBUG(swpcnd_w); DEBUG(swpcnd_deg_in_cl); DEBUG(neigh_cl);
        }

        int best_swpval = Constants::INF;
        int best_trg_cl = -1;

        for( int trg_cl : neigh_cl ){
            int final_swpval = getSwpValForMove( swpval1, swpcnd_w, swpcnd_deg_in_cl,  trg_cl ); // swpval after moving d to cluster trg_cl

            if(debug) DEBUG2(trg_cl,final_swpval);

            if( final_swpval > 0 ){ // positive solution (worse than current one)
                if( !keep_only_nonpositive_candidates ){
                    if( keep_only_best_cluster_to_move_to ){
                        if( final_swpval < best_swpval ){
                            best_swpval = final_swpval;
                            best_trg_cl = trg_cl;
                        }
                    }
                    else res.push_back( constructor(final_swpval, trg_cl) );
                }
            }else{
                // swpval nonpositive
                if( keep_only_best_cluster_to_move_to ){
                    if( final_swpval < best_swpval ){
                        best_swpval = final_swpval;
                        best_trg_cl = trg_cl;
                    }
                }
                else res.push_back( constructor(final_swpval, trg_cl));
            }
        }

        // the second condition is to make sure that we do not add invalid data in case no nonpositive swap
        // candidate exists.
        if( keep_only_best_cluster_to_move_to && ( best_trg_cl != -1 ) ){
            res.push_back( constructor(best_swpval, best_trg_cl) );
        }
    }

    /**
     * If true, then the creator will discard (at any time) all swap candidate with negative swap value.
     */
    bool keep_only_nonpositive_candidates = true;

    /**
     * If a swap candidate can be moved to many clusters, then:
     * if [keep_only_best_cluster_to_move_to] is true, then only the best one cluster to move to will be kept. If there
     * is a tie, a 'random' one will be retained (probably the first met)
     *
     * If false, then all candidates will be kept (where candidates are clusters that are incident to some considered
     * at the moment structure). Bear in mind, that this can lead to O(clusters.size()) candidates for each considered
     * swap candidate. This can lead to O(N^2) time and memory complexity.
     *
     * Mark that this function will not by default consider and return all other clusters to which given swap candidate
     * can be moved. For example, if a single node can be removed from its cluster and moved to any other cluster with
     * positive swap value, it does not mean, that every such possibility will be found. In order to find all possible
     * movement of the swap candidate to other clusters, set [find_possible_movements_to_all_clusters]. Note also, that
     * find_possible_movements_to_all_clusters may not be respected by every SwapCandidateCreator.
     */
    bool keep_only_best_cluster_to_move_to = true;

    /**
     * If true, then for each swap candidate all other clusters will be checked. This can lead to O(clusters.size())
     * possible movements of a single swap candidate and will surely result in a (possibly very significant) increase
     * in running time of the algorithm. This parameter is dependent on the specific instance of SwapCandidateCreator.
     */
    bool find_possible_movements_to_all_clusters = false;
};

#endif //ALGORITHMSPROJECT_SWAPCANDIDATE_H
