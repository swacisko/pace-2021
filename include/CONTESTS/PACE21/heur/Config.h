//
// Created by sylwester on 3/8/21.
//

#ifndef CESWAT_CONFIG_H
#define CESWAT_CONFIG_H

#include "State.h"

/**
 * Represents swap candidate creators
 */
enum SwpCndCrId{
    node = 0, // SwpCndNode
    edge_same_cl = 1, // SwpCndEdge, only edges in the same cluster
    edge_diff_cl = 2, // SwpCndEdge, only edges with ends in different clusters
    edge_all = 3, // (1) + (2)
    triangle = 4, // SwpCndEdgeTriangle
    exp_ord = 5, // ExpansionOrder
    exp_ord_rep = 6, // ExpansionOrderRepulsion
    exp_ord_attr = 7 // ExpansionOrderAttraction
};

/**
 * Represents mode for updating state - applying swap
 */
enum SwpCndSwapMode{
    GREEDY_MAXIMAL_DISJOINT = 0, // candidates are sorted by non-decreasing swap value. Then, greedily, candidates are
    // applied if possible. It is possible to apply considered candidate if none of applied earlier candidates has
    // affected a any cluster affected by current candidate

    ONLY_BEST_ONE = 1 // only the best candidate (with least swap value) is applied. This may 'converge' very slowly
            // to local optimum.
};

/**
 * Speed mode used for algorithms. The slower, the more swap candidates are created and the greater (in theory) chance
 * of finding a better solution.
 */
enum SPEED_MODE{
    super_fast,
    very_fast,
    fast,
    medium_fast,
    medium,
    medium_slow,
    slow,
    very_slow
};

/**
 * Mask values must be used - powers of 2
 */
enum COARSEN_MODE{
    remove_edges = 1, // all edges that have both ends always in different clusters will be removed
    contract_all = 2, // all nodes that are always in the same cluster will be 'contracted'
    contract_matching = 4 // for all sets X of nodes that are always in the same cluster, a random maximal matching
            // will be created and that pairs of nodes will be contracted.
};

/**
 * Represents some configuration for the algorithms, e.g. values for different parameters.
 *
 * The main algorithm works in 'Large' iterations.
 *
 * Before large iterations a cluster graph is created for given partition.
 * In each 'Large' iteration:
 * 1. A local search is called for the cluster graph. The LS method consists in:
 * a) first create an initial state using some fast heuristic.
 * b) improve created state by calling each of swap candidate creators specified in [swpCndCreatorsToUse] and updating
 * the state if possible.
 * After [max_nonnegative_iterations] iterations that did not improve the state's result, the LS finishes and the best
 * solution found is returned.
 * That best solution is further improved using SwpCndNode and SwpCndEdge called to the
 * original graph.
 *
 * After [granularity_frequency] Large iterations are done (and therefore exactly granularity_frequency solutions
 * known), we update the partition and the initial graph. The partition is updated using the following rules (or only
 * some of them):
 * a) if there exists an edge (u,v) such that u and v are in different clusters in each of known solutions, then remove
 * that edge from the graph
 * b) if there exist two nodes (u,v) that are always together in a cluster in all known solutions, then 'merge' those
 * nodes by marking their partition to the same cluster (they will be considered as single node in cluster graph in all
 * further Large iterations).
 * c) for each set X of nodes that are always in the same cluster in known solutions, select a random maximal matching
 * of that 'cluster'. The matching is selected, by greedily selecting an edge (a,b,w) that minimizes
 * nw_a * nw_b - w
 */
class Config{
public:
    Config();

    /**
     * Creates configuration to run.
     */
    void createConfiguration();

//private:

    /**
     * Swap candidate creators that are to be used. They should be used in the exact order specified by the vector.
     */
    vector<SwpCndCrId> swpCndCreatorsToUse;

    /**
     * This is the maximal recursion depth used in [Solver::run_recursive].
     */
    int max_recursion_depth = 5;

    /**
     * Function that creates initial sets for expansion order creation. For each of those sets an independent expansion
     * order will be created.
     * @return a vector containing 'sets' for which expansion orders are to be created. These sets contain ids from the
     * INDUCED graph, that is from range [0, cl->size())
     */
    static VVI expansionOrderInitialSetProvider(Cluster* cl);

    /**
     * SwpCandidate creators are run in a series, one after another. If at some point some creator returns at least
     * one swap candidate with negative value, then:
     *
     * If [apply_swap_on_first_negative] is true, state is updated without running the rest of creators
     * from [swpCndCreatorsToUse].
     * If [apply_swap_on_first_negative] is false, then state will be updated after all creators
     * from [swpCndCreatorsToUse] are run.
     */
    bool apply_swap_on_first_negative = true;

    /**
     * Type of swaps made to apply current state.
     */
    SwpCndSwapMode swap_application_mode = GREEDY_MAXIMAL_DISJOINT;

    /**
     * If true, then swap candidates with swap value 0 will be applied only if there exists a cluster to which a node
     * is about to be moved, with smaller size than the cluster from which it is moved.
     */
    bool apply_neutral_swaps_only_to_smaller_clusters = false;

    /**
     * Maximum number of iterations that are to be made since the last time a state was IMPROVED (by a candidate with
     * negative swap value).
     * An iterations consists of calling swap candidate creators specified in [swpCndCreatorsToUse].
     */
    int max_nonnegative_iterations = 1;

    /**
     * Number of 'Large' iterations of the main algorithms after which partition or graph will be updated.
     */
    int granularity_frequency = 10;

    /**
     * If true, then in the constructor of Solver object, kernelization algorithms will be used to create initial
     * partition. Kernelization will be run not longer (roughly, until current function ends its work) than
     * [max_kernelization_time_in_sec] seconds.
     */
    bool use_kernelization = false;
    bool use_only_fast_exact_kernelization = true; // uses only rules 1, 15 and 16
    bool use_heuristic_kernelization = true;
    double max_kernelization_time_in_sec = 60;

//    ******************* NEG
    /**
     * Maximum number of perturbations done by NEG
     */
    int neg_max_perturb = 10;
    /**
     * Maximum number of nonnegative iterations done by NEG
     */
    int neg_max_nonneg_iters = 10;

    /**
     * if true, then edge moves will be checked (in E * sqrt(E))
     */
    bool neg_use_edge_swaps = true;
    int neg_edge_swaps_frequency = 15;
    bool neg_use_edge_repulsion = true;

    double neg_perm_fraction = 1.0;

    int neg_max_best_cl_size_triangle_swaps = 2;
    bool neg_use_triangle_swaps_to_other_clusters = true;

    /**
     * Only clusters of that size or large will be considered in Solver::localSearch in ComponentExpansionRepulsion
     * and ComponentExpanionAttraction
     */
    int min_cluster_size_for_eo_rep = 3;
    int min_cluster_size_for_eo_attr = 4;

    /**
     * If 0, then triangle moves will not be checked.
     * If 1, then triangle moves will be checked only to an empty cluster (in E * sqrt(E))
     * If 2, then triangle moves will be checked to all clusters (TIME-CONSUMING!!)
     */
    int neg_use_triangle_swaps = 1;
    int neg_triangle_swaps_frequency = 50;
    int neg_chain2_swaps_frequency = 29;
    int neg_node_interchanging_frequency = 25;

    /**
     * If true, then node interchanging between two clusters will be used (swapping nodes)
     */
    bool neg_use_node_interchange = true;

    /**
     * If true, then cluster joining will be checked in each iteration of NEG
     */
    bool neg_use_join_clusters = false;

    bool neg_use_queue_propagation = true;

    bool neg_use_chain2_swaps = true;

    bool neg_use_two_node_swaps = true;

    int neg_move_frequency = 2;

    /**
     * If true, then NEG_Nomap and NEG_W1 will not be used (as they are much slower on some graphs).
     */
    bool use_neg_map_version = false;

    int neg_max_iterations_to_do = 1e9;

    //    ******************* NEG

    bool solver_use_only_neg_to_create_known_solutions = true;
    bool solver_improve_best_known_solution_using_local_search = true;

    int solver_max_rec_depth_run_fast = 3;

    /**
     * Those values will be propagated to swap candidate creators.
     */
    bool keep_only_nonpositive_candidates = true;
    bool keep_only_best_cluster_to_move_to = true;

    /**
     * If true, then some algorithms will be run with parameters reducing complexity
     */
    SPEED_MODE speed_mode = medium;

    /**
     * Sets all parameters that rely on [speed_mode].
     * @param mode
     */
    void setSpeedMode( SPEED_MODE mode );

    /**
     * This parameter will be forwarded to    SwpCndTriangleCreator::create_MoveTo_SwapCandidates
     */
    bool use_only_empty_cluster_in_swp_cnd_triangle = false;

    bool use_only_common_neighbors_in_swp_cnd_edge = false;

    /**
     * If true, then in Solver, no candidates will be removed until the end of smallIteration(). It will enable
     * possible future algorithms such as chainAlgorithm to operate on present candidates.
     * If false, all candidates will be removed as soon as they were considered for swap and not applied.
     */
    bool keep_all_swap_candidates = false;

    /**
     * Type of state initialization. By default it is RANDOM_MATCHING
     */
    StateInitializationType state_init_type = RANDOM_STATE_PERM;

    /**
     * This is a MASK!! That is if [coarsen_mode] & option is > 0, then the option will be applied.
     */
    int coarsen_mode = contract_matching; //(contract_all | remove_edges);

    /**
     * If true, then the first solution created in run_fast will be induced from the best solution from lower
     * recursion level.
     * This ensures that at greater recursion depths we will get at least one result that is not worse than currently
     * found best one.
     */
    bool solver_run_fast_induce_first_solution_from_lower_levels = false;

    bool neg_do_not_perturb_if_improved = false;
};

#endif //CESWAT_CONFIG_H
