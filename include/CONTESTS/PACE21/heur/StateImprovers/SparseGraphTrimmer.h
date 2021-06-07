//
// Created by sylwester on 5/18/21.
//

#ifndef ALGORITHMSPROJECT_SPARSEGRAPHTRIMMER_H
#define ALGORITHMSPROJECT_SPARSEGRAPHTRIMMER_H

#include <CONTESTS/PACE21/heur/Config.h>
#include "../State.h"
#include "../SwapCandidates/SwapCandidate.h"

class SparseGraphTrimmer{
public:
    SparseGraphTrimmer( State * st );

    /**
     * Creates the state. It uses the following observations:
     * 1. If there are two node u and v with deg(u) = deg(v) = 1 with common neighbor d with deg(d) <= 4, then makes
     * the cluster (u,v,d) a single cluster (and 'removes' edges to proceed further)
     *
     * 2. If there exists an edge (u,v) with deg(u) = deg(v) = 2 and both u and v are incident to d with deg(d) <= 4,
     * then make (u,v,d) a single cluster.
     *
     * 3. If none of the above applies, then takes a leaf v with N(v) = u and makes (u,v) a single cluster.
     *
     * @return
     */
    State& createState();


//private:

    /**
     * Pointer to the cluster graph
     */
    ClusterGraph * clg;

    /**
     * Pointer to the state that is to be created - it should aleardy have initialized all arrays
     */
    State* st;


    /**
     * Considers nodes in non-descending degrees, then for each considered node u makes u \cup N(u) a cluster and
     * removes those nodes from graph.
     * Works on unweighted graph structure
     */
    void contractNodesNondescendingDegree();

    /**
     * For each induced path, it takes greedily 3 consequtive nodes on that path, then makes it a single cluster. Then
     * those nodes are removed from graph and procedure is repeated.
     *
     * #CAUTION! With simple (current) implementation it may work in O(N^2) time. Works only for VVI structure
     * [clg->origV].
     */
    void greedyPathClustering();


    void cherryRemover();
};

#endif //ALGORITHMSPROJECT_SPARSEGRAPHTRIMMER_H
