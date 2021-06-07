//
// Created by sylwester on 3/13/21.
//

#ifndef ALGORITHMSPROJECT_CRITICALCLIQUE_H
#define ALGORITHMSPROJECT_CRITICALCLIQUE_H

#include "Makros.h"

class CriticalClique{
public:

    // default constructor
    CriticalClique() : V(nullptr), id(-1){}

    CriticalClique( VVI& V, VI nodes, int id );

    CriticalClique( const CriticalClique& oth ) { if( this != &oth ){V = oth.V; nodes = oth.nodes; id = oth.id;} }

    static void test();

    int size(){ return nodes.size(); }

    friend ostream& operator<<( ostream& str, CriticalClique& cc );

//private:
    /**
     * Pointer to the graph, where [nodes] form a critical clique
     */
    VVI * V;

    /**
     * nodes in the critical clique
     */
    VI nodes;

    /**
     * id of the clique
     */
    int id;
};

#endif //ALGORITHMSPROJECT_CRITICALCLIQUE_H
