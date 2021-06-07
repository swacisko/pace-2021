//
// Created by sylwester on 3/17/21.
//

#ifndef ALGORITHMSPROJECT_GRAPHTRIMMER_H
#define ALGORITHMSPROJECT_GRAPHTRIMMER_H

#include "Makros.h"

/**
 * Class responsible for finding ('trimming') some structure in a graph.
 */
class GraphTrimmer{
public:
    GraphTrimmer(VVI& V);

    /**
     * Finds all maximal subgraphs that are paths.
     * Each of the returned vectors res[i] contains nodes on a path.
     *
     * Each node on each path has degree <= 2.
     * If a graph contains cycle component, than all nodes of the cycle are on the path exactly once
     *
     * If a path has only one end node with degree 1, then it will be as the first element.
     *
     * @return
     */
    VVI findMaximalPathsAndCycles();


//    static void test();

private:

    VB was; // helper array to mark visited nodes

    VVI *V;
    int N; // V.size();

    /**
     * Starts creating a path from given node.
     * Node [num] must be of degree 2.
     *
     * It moves along the path until node with degree != 2 is met (or [num] again in case of a cycle).
     *
     * @param num
     * @param path
     * @param was
     */
    void addToPath( int num, VI & path, VB & was );

    /**
     * Finds the path on which node num lies. [num] must be of degree 2.
     * @param num
     * @return path that contains node [num]
     */
    VI getPathWithNode( int num );
};


#endif //ALGORITHMSPROJECT_GRAPHTRIMMER_H
