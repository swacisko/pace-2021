//
// Created by sylwester on 3/8/21.
//

#include "CONTESTS/PACE21/heur/Cluster.h"

VI Cluster::helper(1,-1);

Cluster::Cluster( ClusterGraph & v, VI nodes, int id ) {
    if( !nodes.empty() ){
        g = induce(v,nodes, helper);
        cluster_weight = accumulate( ALL(g.node_weights),0 );
    }
    this->id = id;
}

ostream& operator<<( ostream& str, Cluster& cl ){
    str << "[ Cluster, id: " << cl.id << endl << "weight: " << cl.cluster_weight << endl
        << "Induced cluster graph: " << cl.g << endl;
    return str;
}