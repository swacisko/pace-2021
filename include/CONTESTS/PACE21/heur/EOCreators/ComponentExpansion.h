//
// Created by sylwester on 3/9/21.
//

#ifndef ALGORITHMSPROJECT_COMPONENTEXPANSION_H
#define ALGORITHMSPROJECT_COMPONENTEXPANSION_H

#include <CONTESTS/PACE21/heur/ClusterGraph.h>
#include <CONTESTS/PACE21/heur/ExpansionOrder.h>
#include <CONTESTS/PACE21/heur/SwapCandidates/SwapCandidate.h>
#include "datastructures/Heap.h"
#include "Makros.h"

class ComponentExpansion{
public:
    ComponentExpansion( Cluster & cl );

    virtual ~ComponentExpansion(){}

    /**
     * Creates an expansion order, by expanding set S. It uses rules in order they are provided in [cmp_rules].
     *
     * If [use_heap] is set to true, then most data will be updated using heap. This will make the complexity of
     * O(E*logN). However, in case of dense graphs such as clusters should be, we can get O(E + N^2) by considering
     * always all nodes outside S as piercing nodes. In such a case set [use_heap] to false. NOT using heap also
     * enables us to apply easily many other rules of selecting piercing nodes (e.g. rule 10).
     * By default [use_heap] is set to false. Set it manually to true if running that on large sparse graphs.
     *
     * @param max_nodes_in_eo this is the maximal number of nodes that will be searched for. After S grows to that
     * value, expansion is terminated.
     */
    ExpansionOrder getExpansionOrder( VI A, const bool use_heap = false, const int max_nodes_in_eo = 1e9 );

    /**
     * Sets rules to be used during comparison.
     * @param rules
     */
    void setCmpRules( VI rules ){ cmp_rules = rules; }

//protected:
    /**
     * Pointer to the graph, where the expansion order is created
     */
    InducedClusterGraph * clg;

    Cluster* cl; // pointer to the cluster

    VVPII * V; // V = &clg.V

    int N; // clg.V.size()

    /**
     * sumEW[v] is the sum of weights of all edges of node v
     */
    VI sumEW;

    /**
     * eInS[v] is the SUM OF WEIGHTS OF EDGES of node v in set S
     */
    VI eInS;

    VB inS; // boolean marker for S
    int sumNWinS; // sum of node weights of nodes in S

    /**
     * Expanded set.
     */
    VI S;

    int cut_value; // sum of weights of edges with one end in S and the other outside S.

    /**
     * If true, then if a node should be added to expanded set S, and that node has less connection to S than
     * 3*sumNWinS / 4, then expansion is terminated.
     */
    bool terminate_on_cluster_violator = false;

    /**
     * Heap used to update data of vertices. It is recommended to use it in sparse graphs (if possible to use proper
     * node selection criteria - e.g. rule10 is not possible to be updated using heap).
     */
    Heap<int> heap;

    void moveNodeToS(int d, bool use_heap = false);

    //****************************************************
    /**
     * Number of comparison rules.
     * By AVERAGE we mean, that all absolute values (e.g. edges outside S) are divided by
     * the weight of the node in the cluster graph.
     * By AVERAGE2 we mean that all absolute values \sum w_i of edges (d, p_i,w_i) are divided by
     * ( nw(d) * ( \sum nw(p_i) ) ), where p_i are some neighbors of d (those that are in S in case of rule 7, those
     * that are not in S in case of rule 8).
     *
     * 1. First rule: select node with tightest (sum of edge weights) connection to S.
     * 2. Second rule: select node with least weight of edges to neighbors outside S.
     * 3. Third rule: select node that minimizes X - Y, where X is weight of edges to neighbors outside S and Y is
     * sum of weights to neighbors in S. Because X = deg - Y, we have X-Y = deg - 2*Y
     *
     * 4. Fourth rule: select node with AVERAGE tightest connection to S.
     * 5. Fifth rule: select node with AVERAGE least number of neighbors outside S.
     * 6. Sixth rule: select node that minimizes AVERAGE X - Y, where X is number of neighbors outside S and Y is number
     *  of neighbors in S. Because X = deg - Y, we have X-Y = deg - 2*Y
     *
     *  10. Tenth rule: Select the node that, when added, minimizes ratio   cut_value / all_edges, where
     *  cut_value is the number of edges between node in S and node outside S, and all_edges is the maximum number of
     *  such edges (that is exactly  Sweight * ( cluster_weight - Sweight ).
     *  The idea of this rule is to find a subset S of given cluster C that minimizes the density between S and V\S.
     *  Only such subsets for which the density is <= 1/2 can be removed from C to make a disjoint cluster (moving S
     *  from C to other cluster may still be profitable, even if density is grater than 1/2).
     *  DOES NOT WORK WITH HEAP.
     *
     *  11. Eleventh rule: Select node, that when added, minimizes value 2*cut_value - all_edges  (see rule10 for
     *  description). This is somehow similar to rule10. This time, however, we count the absolute value of swpVal - in
     *  current state we have to insert    Sweight * ( cluster_weight - Sweight ) - cut_value    edges to make cl a
     *  cluster. After making S a disjoint cluster, we will have to add   cut_value   edges. The difference is taken as
     *  a comparator. The lower, the better.
     *  DOES NOT WORK WITH HEAP.
     */

    /**
     * When selecting node to add, we compare first nodes by rule cmp_rules[0], then by cmp_rules[1], etc.
     */
    VI cmp_rules;

    /**
     * Function used to compare two 'piercing' nodes [a] and [b]. The node that is smallest with respect to this
     * function will be added to S, that is if cmpFun(a,b) == true, then a is better than b.
     */
    virtual bool cmpFun( int a, int b );
    //****************************************************

};

#endif //ALGORITHMSPROJECT_COMPONENTEXPANSION_H
