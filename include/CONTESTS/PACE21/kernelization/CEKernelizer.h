//
// Created by sylwester on 3/8/21.
//

#ifndef CESWAT_CEKERNELIZER_H
#define CESWAT_CEKERNELIZER_H

#include <utils/RandomNumberGenerators.h>
#include <Constants.h>
#include "Makros.h"
#include "CriticalClique.h"
#include "datastructures/FAU.h"

/**
 * Class represents an object that performs kernelization of given graph.
 * Rules according to "An overview of different kernelization algorithms for the cluster editing problem"
 * by Joost Besseling.
 *
 *
 * CAUTION!!
 * This algorithm may not be deterministic. In order to make it deterministic comment sort in rule1().
 * However, this sort seems to magically make a difference to the results (even though ccs are
 * randomly shuffled beforehand)
 * CAUTION!!
 *
 */
class CEKernelizer{
public:

    CEKernelizer(VVI &V, unsigned seed = 179'162'867);

    static void test();

    bool rule1();

    bool rule2();

    bool rule3();

    bool rule4();

    bool rule5( int k ); // requires parameter k

    /**
     * Lemma 3 is stirnger than just rule6, so instead of rule 6, lemma 3 is implemented
     *
     */
    bool rule6_lemma3();

    /**
     * Can rule7 be executed simultaneously with rule8 ??
     */
    bool rule7();

    /**
     * Can rule8 be executed simultaneously with rule7 ??
     */
    bool rule8();

    /**
     * Pendulum algorithm
     */
    bool rule9();

    /**
     * Similar neighborhood rule in unweighted case. (Exact Algorithms for Cluster Editing: Evaluation and Experiments)
     *
     * If there exists an edge (a,b) with N(b) \subset (N(a) \ b) and |N(b)| = |N(a)|+1, then a and b belong to the
     * same cluster in some optimal solution.
     *
     * CAUTION!! Structure of V is modified - neighborhood lists are sorted! Edges may be added to the graph.
     *
     * Is the rule correct?
     *
     * @return
     */
    bool rule15();

    bool rule16(const bool debugRuleApplicationOnTheFly = false); // almost clique rule for triangles, K4 and diamonds (it is correct)

    /**
     * Consider a set A = {K_0, K_1, ... } of critical cliques that have the same OPEN neighborhood N(K_i).
     * Let G1 be a graph after removing all K_i from V.
     *
     * If now there exists an index i and a critical clique K' in G1 such that in V we have N(K_i) \subset K'
     * and |N(K_i)| <= |K'| / 2 and |K_i| < |K'| / 2,  then make K_i a single cluster and remove it from V.
     *
     * Is it correct?? If not, why?
     */
    bool rule_heur_1();

    /**
     * Consider a set A = {u_0, u_1, ...} of nodes for which N(u_i) is a clique.
     * Let G1 be a graph after removing from V all nodes from A.
     *
     * If there exists in G1 a critical clique K such that in V there exists node u \in A with property
     * N(u) \subset K and |N(u)| <= |K| / 2, then make u a single cluster and remove it from V.
     *
     * @return
     */
    bool rule_heur_2();

    /**
     * Rule similar to rule2.
     *
     * Consider a set A = {K_0, K_1, ...} of critical cliques for which N(K_i) is a clique.
     * Let G1 be a graph after removing from V all cc's from A.
     *
     * If there exists in G1 a critical clique K such that in V there exists cc K_i \in A with property
     * N(K_i) \subset K and |N(K_i)| <= |K| / 2, then make u a single cluster and remove it from V.
     *
     * If [apply_strict_size_condition] is set, then also condition |K_i| <= |K| / 2 must be met.
     * If [apply_strict_size_condition] is NOT set, then condition |K_i| <= |K| must be met.
     *
     * @param apply_strict_size_condition
     * @param min_cc_size minimum size of a critical clique that will be considered to remove
     * @param max_cc_size maximum size of a critical clique that will be considered to remove. Taking 1,
     * we should get the same as using rule2.
     * @return
     */
    bool rule_heur_3(int min_cc_size=1, int max_cc_size=Constants::INF, bool apply_strict_size_condition = true);

    /**
     * Makes iterations calling rule_heur_3 by changing parameters.
     *
     * FIXME: enable changing values of MAXIMAL_MIN_S and MAX_S_RELATIVE_RANGE in function definition.
     * @return
     */
    bool iterative_rule_heur_3();
    bool USE_HEUR_RULE_3_WITHOUT_STRICT_SIZE_CONDITION = false;

    /**
     * LOSSY KERNELIZATION RULE.
     * It is possible, that every pair created using this rule will increase the optimal result by 1.
     *
     * Let A = {v_1,...,v_a} be a set of nodes with the same OPEN neighborhood N(v_i). Then, we take the first
     * |A| * perc nodes and pair them, {v_1,v_2}, {v_3,v_4}, ...
     *
     * Each pair is marked to be in the same cluster, but the edge IS ([add_edge]) / IS NOT (![add_edge]) added to the
     * structure.
     * The idea behind is that the two nodes will either be in the same cluster, or will end as two separate nodes
     * in the final result. Adding the edge may help do other heuristic rules.
     *
     * By taking perc = 0.99 we pair all nodes, leaving 1 node (if |A| is odd) or 2 nodes (if |A| is even).
     */
    bool rule_heur_4( double perc = 0.99, const bool add_edge = false );
    /**
     * The collection of paired nodes, to be able to restore the 'lossy' kernelization at the end and try to improve
     * the result.
     */
    VPII rule_heur_4_paired_nodes;

    /**
     * LOSSY KERNELIZATION RULE.
     * This may pessimistically cause large deviations from optimal result!
     *
     * If there exist an edge(u,v) with
     *
     * N(u) \subset N(v) and |N(u)| > |N(v)| * perc
     * and
     * |N(v)| - |N(u)| <= absolute_diff
     * and | N(u) \cap N(v) | >= min_common_neighbors
     *
     * then u and v are marked to be in the same cluster. Setting perc=0 and max_absolute_diff = 1 we get rule15.
     *
     * It may make sense to add any restriction on the structure of N(u) (e.g. to contain at least
     * |N(v)| - |N(u)| edges).
     */
    bool rule_heur_5(double perc = 0.75, int max_absolute_diff = 30, int min_common_neighbors = 3);

    /**
     * If there exists a triangle (v_1,v_2,v_3) with deg(v_i) <= 2 + [max_out_deg] and
     * sum_deg(v_i) <= 6 + [max_out_deg_sum], then mark all nodes with degree <= 2 + [max_out_deg_to_merge] to the
     * same cluster.
     *
     * Default values let consider only triangles with at most 3 out-of-triangle edges and at most 2
     * out-of-triangle-edges-per-node.
     *
     * Is this rule heuristic for max_out_deg_to_merge=1  and max_out_deg_sum = 3 or 4?
     * For  max_out_deg_sum = 2, it is the almost clique rule applied to triangle, so not heuristic
     */
    bool rule_heur_6(int max_out_deg_sum = 4, int max_out_deg = 2, int max_out_deg_to_merge = 1 );

    /**
     * If there exists a cycle C4 (v_1,v_2,v_3,v4) with maximum out-degree of a node not more than max_out_deg and
     * sum of outdegrees not more than max_out_deg_sum, then all nodes with outdegree at most max_out_deg_to_merge
     * will be marked to be in the same cluster.
     *
     * By setting max_out_deg_sum = max_out_deg = max_out_deg_to_merge = 3 and enable_diamonds = false, we obtain the
     * almost critical clique for K4. IF we set (2,true) instead of (3,false), then we obtain almost critical clique
     * for diamonds.
     *
     * For 'diamond' check, we should rather set max_out_deg_sum = 3, not to risk 'bulging clusters'.
     */
    bool rule_heur_7( int max_out_deg_sum = 4, int max_out_deg = 2, int max_out_deg_to_merge = 1,
                      bool enable_diamonds = true );

    /**
     * Applies exhaustively rules 1-9. If at any point rule is applied, then restarts kernelization from
     * rule 1.
     */
    void fullKernelization(bool use_heuristic_rules, int additional_randomized_iterations);

    /**
     * Applies heuristic rules to current state of the graph. If heuristic rules make any difference to the graph,
     * then runs fullKernelization again with heuristic rules admitted.
     * This is used to obtain both states after kernelization (with and without heuristic rules) faster than calling
     * fullKernelization(true) again.
     */
    void improveKernelizationUsingHeuristicRules();

    /**
     * Applies heuristic rules.
     * @return true if any change was done, false otherwise
     */
    bool applyHeuristicKernelization(bool debugRuleApplicationOnTheFly = false);

    /**
     * Marks given nodes as coming from the same cluster, then removes those nodes from [V].
     * @param nodes
     */
    void makeCliqueAndRemoveFromGraph( VI nodes );

    /**
     * By empty critical cliques we mean a maximal set of nodes with N(u) == N(v) - neighborhoods need
     * not be closed.
     * @return res[i] is a pair ( size, cnt ), where cnt is the number of ACC of given size
     */
    vector< PII > getEmptyCriticalCliqueSizes();

    /**
     * By empty critical cliques we mean a maximal set of nodes with N(u) == N(v) - neighborhoods need
     * not be closed.
     *
     * @return vector contating sets of nodes from the same ACC.
     */
    VVI createEmptyCriticalCliques();

    /**
     * Adds missing edges between all nodes marked to the same cluster.
     */
    bool fillMissingEdgesBetweenClusterNodes();

    void setDisableRule(int x, bool b = true ){ assert(x <= MAX_RULES); disabled_rules[x] = b; }
    void setDisableHeurRule(int x, bool b = true ){ assert(x <= MAX_HEUR_RULES); disabled_heur_rules[x] = b; }

    /**
     * if true, then heuristic rules will be called exhaustively in function applyHeursticRules(), until no further
     * heuristic improvement can be done. Otherwise, after first improvement the function will return true.
     */
    bool use_heuristic_rules_separately = false;

//private:
    /**
     * Copy of the graph for which kernelization is made
     */
    VVI V;

    /**
     * N = V.size();
     */
    int N;

    /**
     * hashes[i] is a random long long. Used to create critical cliques
     */
    VLL hashes;

    void createHashes();

    /**
     * inV[i] is true if i is still in graph V. If not, it will not be considered for kernelization.
     */
    VB inV;

    /**
     * inCl[i] is the id of the cluster, where node i belongs.
     * Initiallly should be inCl[i] = i     *
     */
    VI inCl;

    /**
     * This is almost the same as [inCl], but uses Find and Union - therefore support removing nodes from
     * critical clique, just as in rule 9.
     */
    FAU fau;

    VB helper1;

    std::mt19937_64 drng;
    UniformIntGenerator rnd;

    static const int MAX_RULES = 100;
    static const int MAX_HEUR_RULES = 100;

    /**
     * ruleAppliedCnt[i] is the number of time a given rule applied to the graph.
     * Just for statistics.
     */
    VI ruleAppliedCnt;

    /**
     * If false, then given rule will not be executed
     */
    VB disabled_rules, disabled_heur_rules;

    /**
     * Critical cliques used during kernelization
     */
    vector<CriticalClique> ccs;

    /**
     * Creates a vector containing critical cliques of the graph [V] (passed as parameter) induced by given [nodes].
     * If nodes are empty (by default), then the whole graph V is considered
     * Critical cliques are randomly shuffled, then sorted by their size in non-ascending order.
     * @return
     */
    vector<CriticalClique> createCriticalCliques( VVI & V, VI nodes = {} );

    /**
     * Creates a vector containing critical cliques of the graph [V] induced by given [nodes].
     * If nodes are empty (by default), then the whole graph V is considered.
     *
     * ATTENTION!! Marks all nodes from the same critical clique to the same cluster (using [fau]).
     * Critical cliques are randomly shuffled, then sorted by their size in non-ascending order.
     *
     * @return
     */
    vector<CriticalClique> createCriticalCliques( VI nodes = {} );

    /**
     * Finds all nodes that are within [dst] distance from [nodes].
     * Uses [helper1] array.
     * @param nodes
     * @return
     */
    VI getNeighborhood( VI & nodes, int dst );

    /**
     * Returns the editing degree of node v with respect to critical clique K.
     * Vector [inN1] must have inN1[i] = true if and only if i belongs to N(K).
     * Vector [inK] must have inK[i] = true if and only if i belongs to K.
     *
     * v must belong to N(K) (that is N1)
     *
     * Works in time O( N(v) )
     * @param v node from N(K)  (that is from N1)
     * @param inK marker array for K, with ALREADY MARKED nodes from K
     * @param inN1 marker array fro N1, with ALREADY MARKED nodes from N1
     * @return
     */
    int getEditingDegree(int v, VB & inK, VB & inN1, int N1size );

    /**
     * Marks the state of clusters, just as it is in [fau]. In other words makes just
     * inCl[i] = fau.Find(i)
     */
    void markInClVector();
};

#endif //CESWAT_CEKERNELIZER_H
