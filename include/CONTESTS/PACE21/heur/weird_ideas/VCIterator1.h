//
// Created by sylwester on 12/20/20.
//

#ifndef PACE21_VCITERATOR1_H
#define PACE21_VCITERATOR1_H

#include "Makros.h"

class VCIterator1{
public:


    /**
    * Creates helper graph and runs VC solver on it. Then it extracts information, and edits edges in [V] correspondingly.
    * [V] is not modified, a new graph is created and returned.
    * @param g graph to transform
    * @return
    */
    VVI nextIteration( VVI g );


private:

    VVI helperGraph;

    /**
     * Creates a graph in which each node represents a change of one edge in [V].
     * @param V,
     */
    void createHelperGraph(VVI & V);

    /**
     * nodeToEdgeRemapper[i] is an edge e, such that i-th node in [helperGraph] corresponds to edge e in [V].
     */
    VPII nodeToEdgeRemapper;

    /**
     * edgeToNodeMapper[ {a,b} ] is the id of node in graph [helperGraph] to which edge {a,b} is mapped
     */
    unordered_map<PII,int, pairhash> edgeToNodeMapper;


};

#endif //PACE21_VCITERATOR1_H
