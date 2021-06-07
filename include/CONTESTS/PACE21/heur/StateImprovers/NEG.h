//
// Created by sylwester on 5/24/21.
//

#ifndef ALGORITHMSPROJECT_NEG_H
#define ALGORITHMSPROJECT_NEG_H


#include <CONTESTS/PACE21/heur/Config.h>
#include "../State.h"
#include "../SwapCandidates/SwapCandidate.h"

/**
 * Algorithm works in iterations.
 * In each iteration a random permutation of nodes is chosen. Then, in that order, nodes are tried to be moved to some
 * other cluster. If, at any point, a node with NONPOSITIVE swap value is found, then it is instantaneously moved.
 *
 * The same goes for checking edges.
 *
 * If during whole iteration no improvement was found, algorithm stops.
 *
 *
 * This is a base class from which NodeEdgeGreedy, NodeEdgeGreedyNoamp and NodeEdgeGreedyW1 can inherit to
 * implement their versions of solutions
 */
class NEG{
public:

    NEG( State & st );

    virtual ~NEG(){}

    /**
     * Initializes this object for given state.
     * @param st
     */
    virtual void initializeForState(State & st) = 0;

    virtual void initializeIndependentData(State & st);

    /**
     * Sets parameters used by this state improved that are specified in [cnf]
     * @param cnf
     */
    virtual void setConfigurations( Config & cnf );

    /**
     * Improves the state using algorithm described in the description of this class.
     * @param max_nonnegative_iters maximum number of iterations in which all candidates had swap value >= 0
     */
    virtual void improve();

//private:

    int perturbations_done = 0;

    /**
     * prefer cluster mode is used to direct swaps to clusters of:
     * 0 - do not prefer, just take random
     * 1 - prefer moving to smaller clusters
     * 2 - prefer moving to larger clusters
     */
    int prefer_cluster_mode = 0;

    /**
     * Only this fraction of random permutation will be considered in a single 'queue iteration'
     */
    double perm_fraction = 1.0;

    virtual void shuffleClg();

    /**
     * Shuffles elements using predefined random permutation
     */
    template<class _T> void localShuffle(_T & v);
    int shuffle_ind = 0;
    VI shuffle_seq;

    /**
     * Needs to create some structure that will be returned later by getEdgesToCluster(v).
     */
    virtual void createEdgesToCluster(int v, bool use_sort = true) = 0;
    virtual int findEdgesToCluster( int v, int cl ) = 0;
    virtual VPII getEdgesToCluster(int v) = 0;

    VVI cluster_nodes;
    /**
     * Creates vector [cluster_nodes]
     */
    virtual void createClusterNodes() = 0;
    virtual VI getClusterNodes(int c) = 0;
    virtual int maxClusterId() = 0;

    UniformIntGenerator rnd;

    /**
     * Maximum number of iterations to do until termination.
     * This may be used to create some solution relatively quickly and do not wait until improve() terminates
     */
    int max_iterations_to_do = 1e9;

    /**
     * Array, easier to use than st.clusters[d].cluster_weight.
     */
    VI cluster_weights, inCl, degInCl;

    /**
     * in_queue[v] is true if node v is in queue, waiting for processing
     */
    VB in_queue;
    deque<int> queue;

    /**
     * If true, then instead of taking random permutations, all nodes from affected clusters will be added to queue
     * and processed in that order. Then the queue emties, new permutation will be considered.
     */
    bool use_queue_propagation = true;

    /**
     * Keeps track of the number of nonempty clusters.
     */
    int nonempty_clusters_cnt = 0;
    virtual int countNonemptyClusters();

    /**
     * If triangle swaps are used, then every [triangle_swaps_frequency] iterations since last perturbation triangles
     * will be checked.
     * The same holds for edge_swaps_frequency.
     */
    int triangle_swaps_frequency = 35;
    int edge_swaps_frequency = 19;
    int node_interchanging_frequency = 25; // #TEST
    int chain2_swaps_frequency = 29;
    int join_clusters_frequency = 47;

    /**
     * If false, then instead of doing perturbations we will terminate
     */
    bool allow_perturbations = true;
    int perturb_mode = 0; // only splitting   0 - splitting, 1 - joining

    /**
     * Sets of ids of free clusters
     */
    set<int> first_free_cluster;
    VB is_empty_cluster;

    /**
     * Counts and returns the number of numbers x between 0 and C (C is largest index such that cluster_nodes[C] is
     * not empty) with cluster_nodes[x].empty() == true
     *
     * This can be used to check if
     */
    virtual int countClusterNumerationGaps();

    ClusterGraph *clg; // st.clg;
    int N; // clg.V.size()

    //****************************************  PERTURBATIONS   ************************
    /**
     * Maximum number of perturbations made.
     * A perturbation is made based on the following:
     * 0. Split cluster randomly in two clusters.
     * 1. Join clusters into pairs.
     * 3. Apply ExpOrdRep to a maximal set of mutually non-neighboring clusters
     *
     * X. Apply all swap candidates with swap value <= perturb_swp_thr
     */
    int max_perturbations = 10;
    int perturb_swp_thr = 0;

    int max_nonnegative_iters = 10; // #TEST - original 10

    void perturb(int perturbation_mode);

    /**
     * Splits all clusters into two clusterm randomly
     */
    virtual void splitClustersIntoTwo();

    /**
     * Joins clusters into pairs.
     * Pairs are created by selecting from all pairs of neighboring clusters greedily the
     * pair that minimizes the density of edges between clusters.
     */
    void joinClustersInPairs();

    /**
     * Sorts clusters from largest to smallest, then applies ExpOrdREP to a maximal set of mutually unaffected
     * clusters.
     */
    void applyExpOrdREP();

    //****************************************  PERTURBATIONS   ************************

    /**
     * Checks all nodes from perm[a], perm[a+1], ..., perm[b], then returns best move in the form
     * (node, cluster, swap_value)
     */
    virtual tuple<int,int,int> getBestNodeMoveForRange( VI & perm, int a, int b ) = 0;
    vector<tuple<int,int,int>> best_node_move_results;

    /**
     * For given node v, checks all cluster, to which v has sum of edge weights at least [PERC] that cluster weight.
     * Then for each node u in that cluster, it checks, whether interchanging nodes v and u would be profitable.
     *
     * Additionally it check all edges (v,u) for that change, even if u is in cluster to which there is less than [PERC]
     * weight of edges.
     *
     * CAUTION! Now it is used in outside loop - craeteEdgesToCluster and createClusterNodes are called in each call
     * to getBestInterchangeNodePairForInterval
     */
    virtual SwapCandidateAdapter getBestInterchangeNodePair( int v );
    virtual SwapCandidateAdapter getBestInterchangeNodePairForInterval( VI & nodes_to_check, int a, int b );
    double PERC = 0.4;
    VB checked_for_v;
    bool use_node_interchanging = true;
    bool use_node_interchanging_in_inner_loop = false;
    /**
     * @return swap value of interchanging given two nodes, assuming that there is edge with weight [int edge_vu_weight]
     * between v and u (may be 0 if there is no edge)
     */
    virtual int getInterchangeValue(int v, int u, int w_vu);
    virtual int getInterchangeValue(int &v, int &u, int w_vu, int &cl_v, int &nw_v, int &clw_v, int &cl_u, int &clw_u,
                            int &ev_to_clv, int ev_to_clu);

    virtual SwapCandidateAdapter getBestTriangleAll(int v) = 0;
    VI edges_to_cluster_triangle, weight_ac_triangle;
    int max_best_cl_size_triangle_swaps = 2;
    bool use_triangle_swaps_to_other_clusters = true;

    /**
     * Checks triangles using [getBestTriangleDiffClToMove(v)] starting from each v from nodes_to_check[a:b]
     */
    virtual SwapCandidateAdapter getBestTriangleDiffClToMove(VI & nodes_to_check, int a, int b);

    /**
     * Moves no [v] to cluster [to].
     */
    virtual void moveNodeTo(int v, int to) = 0;

    /**
     * Adds all nodes in cluster with given id [cl_id] to [queue].
     * @param cl_id
     */
    void addClusterNodesToQueue( int cl_id );

    /**
     *  Checks all edges with one end in { perm[a], perm[a+1], ..., perm[b] }, then returns best move in the form
     * (edge, cluster, swap_value).
     *
     * For given x only edges (x,y) are considered, where y has not greater degree than x. Otherwise we could get
     * complexity O(N^2). If both x and y have the same degree, then y is checked if it occurs in perm later than x.
     *
     * This function is called only if there are no node candidates with negative swap value in range perm[a:b]
     */
    virtual tuple<PII, PII, int> getBestEdgeMoveForRange(VI & perm, int a, int b );
    vector<tuple< PII, PII,int >> best_edge_move_results;

    /**
     * If true. then the hull trick will be used, otherwise all neighboring clusters will be checked by 'brute force'.
     */
    const bool USE_HULL_TRICK_IN_EDGE_SWAPS = true;

    /**
     * Moves edge (v,w) to cluster [to]. Calls twice function [moveNodeTo].
     */
//    void moveEdgeTo( int v, int w, int to );

    /**
     * one node per each move_frequency nodes will be moved.
     */
    int move_frequency;

    /**
     * Checks whether there are two clusters that can be joined together with negative/nonpositive swap value
     * @param allow_zero_swpval if true, then also cluster with zero swap value will be joined
     * @return vector of triples with ids of clusters to join and corresponding swap value
     */
    vector<tuple<int,int,int>> getBestClustersToJoin(bool allow_zero_swpval = false);

    /**
     * Calculates swap value for moving node [v] to cluster [trg_cl]. Node v should have [edges_clv] edges in its
     * cluster and [edges_trg] edges in cluster [trg_cl]. [tot_clv_possible_edges] is total possible number of
     * edges between v and clv \ v.
     */
    virtual int swapValueForNode(int v, int trg_cl, int edges_clv, int tot_clv_possible_edges, int edges_trg);


    /**
     * Keeps track of the overall result of current state
     */
    int current_result = 0;

    /**
     * Best partition (and corresponding best result) of original graph [clg.origV]
     */
    VI best_partition;
    int best_result = 1e9;
    int initial_result = 1e9;

    VI helper, helper2; // helper array;
    VB helper_was, helper_was2, helper_was3;

    /**
     * If true, then if there should be perturbations done, they are not done if the result was improved. That is
     * if best_result (equal to current_result) < initial_result, then we do not use perturbations, and terminate.
     * This may be used in Solver::localSearch to enable quick interchanging between ExpOrdCreators and NEGs.
     */
    bool do_not_perturb_if_improved = false;

    /**
     * This is used to reduce adding edges to hull, when finding edge swap candidates.
     * If there are two clusters with the same total weight, then we need to add only the one to which the number of
     * edges is greater.
     */
    VI min_w2_to_cluster_of_given_weight;

    /**
     * If true, then node swaps will be considered, otherwise they will not be done.
     */
    bool use_node_swaps = true;

    /**
     * If true, then edge swaps will be considered.
     * Edge swaps may take time up to E * sqrt(E) for whole iteration.
     * Edge swaps for given block perm[a:b] will be done only if no negative node swaps were made.
     */
    bool use_edge_swaps = true;
    bool use_edge_repulsion = true;

    /**
     * If true, then not only edges will be checks in edge-swaps, but all pairs of nodes in the same cluster will
     * be checked.
     */
    bool use_two_node_swaps = true; // #TEST - originally was not admitted

    /**
     * If true, then some functions will return as soon as the first negative swap-valued change is found. This may
     * speed up performance, perhaps at  a slight cost of results... who knows...
     */
    bool return_on_first_negative_swap = false; // #TEST - originally false

    /**
     * If 0, then triangle moves will not be checked.
     * If 1, then all nodes of 'almost triangle' need to be in different clusters, and move is made only to empty
     * cluster - !!!! this is perhaps not supported yet
     */
    int use_triangle_swaps = 1;

    /**
     * If true, then we will check at the end of each iteration if there are two clusters that can be joined
     * with negative/non-positive swap value.
     */
    bool use_join_clusters = true;

    bool use_chain2_swaps = true;

    /**
     * Pointer to the number of iterations done. This can be used only inernally, as the object to which this points
     * is local in improve() function.
     */
    int * iterations_done;

    /**
     * Calculates result using PaceUtils::EvaluateSolution and compares it to [current_result]
     * @return true if both results are the same, false otherwise
     */
    virtual bool compareCurrentResultWithBruteResult();


    /************************************** CHAIN SWAPS
    /**
     * Functions checks, whether there exist two nodes u and v in different clusters and a cluster C such that it is
     * profitable to move u from cl_u to cl_v and simultaneously move v from cl_v to C.
     *
     * First, a list of node moves is created. Then we check each node u:
     * For each cluster D neighboring to u, to which u has sum of weights at least 0.4 weight of cluster D we consider
     * all nodes in D. For each node v in D it pays to try to move v to the best cluster (different from cl_u) to which
     * it had least swap value.
     *
     * If there exists a negative-valued swap, then it is performed straightaway and affected clusters are marked and
     * not considered for further swap (because the swap values would have to be updated). This way we find a
     * maximal set of swaps with pairwise-disjoint-affected-clusters.
     *
     * @return true if a negative swap value was achieved, false otherwise
     */
    virtual bool makeChain2Swaps();

    /**
     * Resizes some structures to be able to handle more clusters.
     * @param empty_cl
     */
    virtual void resizeStructuresForEmptyCluster( int empty_cl );

};


#endif //ALGORITHMSPROJECT_NEG_H
