//
// Created by sylwester on 3/8/21.
//

#ifndef ALGORITHMSPROJECT_EXPANSIONORDER_H
#define ALGORITHMSPROJECT_EXPANSIONORDER_H

#include "Makros.h"
#include "Cluster.h"

/**
 * Represents an expansion order
 */
class ExpansionOrder{
public:

    /**
     * Crates an expansion order - just stores arrays v and pointer to cluster c
     */
    ExpansionOrder( VI v, Cluster* c ){
        cl = c;
        swap(ord,v);
    }

    friend ostream& operator<<( ostream &str, ExpansionOrder& eo ){
        str << "[ord: " << eo.ord << ", cl: " << eo.cl->id << "]";
        return str;
    }

    VI ord; // just the order
    Cluster* cl; // pointer to the cluster for which the EO was created.
};

#endif //ALGORITHMSPROJECT_EXPANSIONORDER_H
