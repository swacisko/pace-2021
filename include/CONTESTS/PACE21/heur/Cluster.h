//
// Created by sylwester on 3/8/21.
//

#ifndef ALGORITHMSPROJECT_CLUSTER_H
#define ALGORITHMSPROJECT_CLUSTER_H

#include "Makros.h"
#include "ClusterGraph.h"


/**
 * Represents a cluster in the graph
 */
class Cluster{
public:
    static const int invalid_id = -1;

    /**
     * CAUTION!!  This is single-threaded only in such implementation.
     * Static helper array used to quickly induce graph [g].
     */
    static VI helper;

//    Cluster() : cluster_weight(0), id(invalid_id){}
    Cluster() = default;

    /**
     * Constructs Cluster - that is it induces the cluster from given ClusterGraph.
     */
    Cluster( ClusterGraph & v, VI nodes, int id );

    friend ostream& operator<<( ostream& str, Cluster& cl );

    /**
     * Number of nodes in the cluster - size of the induced graph. Weights are NOT taken into account.
     * @return
     */
    int size(){ return g.V.size(); }

    /**
     * Sum of weights of all nodes in the cluster
     */
    int cluster_weight = 0;

    /**
     * Structure of the graph par induced by node only from this cluster (passed in constructor)
     */
    InducedClusterGraph g;

    /**
     * id of the cluster
     */
    int id = invalid_id;
};

#endif //ALGORITHMSPROJECT_CLUSTER_H
