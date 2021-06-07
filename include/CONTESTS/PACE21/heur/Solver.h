//
// Created by sylwester on 4/14/21.
//

#ifndef ALGORITHMSPROJECT_SOLVER_H
#define ALGORITHMSPROJECT_SOLVER_H

#include "State.h"
#include "SwapCandidates/SwapCandidate.h"
#include "Config.h"
#include "StateImprovers/NEG.h"

class Solver{
public:
    /**
     * @param V structure of the original graph
     * @param initial_partition partition of the graph. For original graph it will be identity array. If a [Solver] object will be
     * created 'recursively', the partition may be needed to pass, induced from the 'upper recursion level' cluster
     * graph.
     * @param rec_depth depth of the recurrence
     */
    Solver(VVI & V, VI initial_partition, Config & cnf, int rec_depth = 0);

    ~Solver();

    /**
     * Runs algorithms to find solution.
     * Runs [iters] large iterations. After each cnf->granularity_frequency large iterations, current solution is
     * subjected to coarsening using [granulateSolution].
     */
    void run(int iters = 500);

    //**************************** RUN2 ********************
    /**
     * Other version of running the algorithms.
     * It first creates some set of known solutions using [createKnownSolutions] function.
     * Then it coarsens the graph.
     * Finally, it tries to refine found solution.
     */
    void run_recursive();

    /**
     * Runs quickly to get quickly the results.
     * It creates some very good solution using a very quick heuristic, contract matching and recursively finds
     * solution for the graphs after contraction. At the end it refines solution using NEG.
     *
     * Look function definition to see parameters used in NEG to find solutions before recursion and in refinement.
     */
    void run_fast();

    /**
     * Craetes and returns a pointer to the newly created NEG object. Depending on current state data (such as
     * recurrence depth, graph statistics, etc) the NEG can be either NodeEdgeGreedy, NodeEdgeGreedyNomap or
     * NodeEdgeGreedyW1.
     */
    NEG* createNegForState(State * st);

    /**
     * This function is called to create a set of know solutions.
     *
     */
    void createKnownSolutions();

    /**
     * Creates a set of known solutions using run_fast()
     */
    void createKnownSolutionsUsingRunFast();

    /**
     * Refines result after coarsening. Uses [partition] to refine.
     */
    void refineAfterCoarsening();

    //**************************** RUN2 ********************


    /**
     * @return set of edge modifications necessary to obtain best result found by the solver.
     */
    VPII getModifications();

    /**
     * Executes single Large iteration (see Config class for description). In each large iteration a new solution is
     * created.
     * @return first element contains partition of [origV], second contains partition of [clg]
     */
    pair<VI,VI> largeIteration(int iter_cnt = -1);

//private:

    Config * cnf;

    /**
     * Pointer to the original graph.
     */
    VVI * origV;

    VVI V; // original structure of the graph, initially equal to [origV], but subject to modifications..

    /**
     * Cluster graph created for [V] and [partition].
     */
    ClusterGraph clg;
    int N; // V.size()


    /**
     * This is the partition of [origV] found by solver that created this object. If this is not empty, then the first
     * solution created by run_fast() will induce the solution from this partition.
     * Hence, one solution of best quality so far will be given 'at hand' just to try to improve it further
     */
    VI lower_level_best_partition_to_induce;

    /**
     * Makes [partition] of [clg] from partition [lower_level_best_partition_to_induce] of origV.
     * [clg] needs to be created to do that.
     */
    bool inducePartitionFromLowerLevelPartition();

    /**
     * State created for cluster graph [clg] in each Large iteration. It is modified during local search.
     */
    State *st = nullptr;

    VLL hashes, hashes2; // hashes of integers from range [0,N-1]

    /**
     * Partition.
     * This is fishy. Sometimes it is used as a partition of V and sometimes as a partition of clg->V.
     */
    VI partition;

    /**
     * Writes [depth] tabs to clog.
     * @param depth
     */
    void logSpacing( int depth = 0 );

    /**
     * The depth of the recurrence at which the [Solver] is created. If [run] function will create instances of this
     * class that will be used to find solutions to some subgraphs, the depth of the recurrence may be used to
     */
    int recurrence_depth = 0;

    /**
     * Vector containing know solutions of [clg->origV] This vector will be cleared every cnf->granularity_frequency calls to
     * [largeIteration]. Each Large iterations should provide one solution.
     */
    VVI known_solutions;

    /**
     * The vector of partitions of [clg] rather than [origV]. These solutions should correspond to those in
     * [know_solutions].
     */
    VVI known_clg_partitions;

    /**
     * Updates partition and/or structure of [V] using methods described in class in Config.h (granularity).
     */
    void granulateSolution();

    /**
     * Function used to create cluster graph in each Large iteration. Uses [V] and [partition] to create [clg].
     * Then, for that [clg] a new state [st] is created and used in local search.
     */
    void createClusterGraph();

    /**
     * Executes smallIteration as long as there are any swaps plus additional cnf.max_nonnegative_iterations times.
     * @param iter_cnt number of large iteration executed
     * @return best partition of [origV] that was achieved in the local search procedure and best partition of [clg]
     */
    pair<VI,VI> localSearch(int itert_cnt = -1);

    /**
     * Tries to make any changes to state [st] basing on swap candidates provided by parameter [candidates].
     * Applies all candidates, according to the rule (e.g ONLY_VEST_CANDIDATE), regardless of their swap value.
     * If only nonpositive candidates are to be applied, the positive ones need to be filtered before calling this
     * function.
     *
     * If there are many candidates that are moved to a new empty cluster, then each candidate 'creates' its own new
     * cluster. Hence many candidates from different cluster can be moved simultaneously to new clusters.
     *
     * @return true if state was improved (if there was a candidate with negative swap value), false otherwise.
     */
    bool apply_swap_for_candidates( vector<SwapCandidate*> & candidates );

    /**
     * Makes some perturbations to [state], e.g. using ComponentExpansionRepulsion, or halving clusters, etc. Options
     * are specified by [cnf]. Perturbation will only be made if nonneg_iter_cnt > 0
     * @param nonneg_iter_cnt
     */
    void perturbState(int nonneg_iter_cnt);

    /**
     * Executes all algorithms in cnf.swpCndCreatorsToUse.
     * @param iter_cnt number of small iteration executed in given [localSearch] run
     * @param nonneg_iter number of nonnegative iteration - that is this number of times smallIteration was called
     * since last time state was improved. It may be used for perturbations.
     * @return true if any improvement was done, false otherwise
     */
    bool smallIteration(int iter_cnt, int nonneg_iter_cnt = 0);

    /**
     * Creates the partition of original graph [origV] for given state [st]. Reads clusters from the state, then joins
     * all nodes from the same cluster to common partition.
     *
     * Creates the partition of clg.V into clusters. It is passed as the second argument.
     *
     * If there are nodes in [origV] that are not 'represented' (contained) by any node in cluster graph [clg], then
     * their partition is marked to the value from [partition].
     */
    pair<VI,VI> createPartitionsForGivenState( State &st );

    /**
     * Used to store best result found so far and the partition of [origV] that gives that result
     */
    int best_result;
    VI best_partition;

    /**
     * Compares given partition to best solution and updates it if the new one is better.
     * @return true if [part] is better than [best_partition], false otherwise
     */
    bool compareToBestSolutionAndUpdate(VI part);

    /**
     * local_search_creator_calls[s] is a pair (called,improved) where called is the number of times given swap
     * candidate creator was called, and improved is the number of times a negative swap candidate was found
     * using given creator.
     */
    unordered_map<string,PII> local_search_creator_calls;

};

#endif //ALGORITHMSPROJECT_SOLVER_H
