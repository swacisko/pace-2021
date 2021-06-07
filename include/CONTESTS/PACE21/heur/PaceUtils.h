//
// Created by sylwester on 12/20/20.
//

#ifndef PACE21_PACEUTILS_H
#define PACE21_PACEUTILS_H


#include "ClusterGraph.h"
#include "State.h"

namespace PaceUtils{

    /**
     * Calculates the number of edges that have to be edited, to create clusters given in [partition]
     * @param V graph in the original state
     * @param partition
     * @return
     */
    extern long long evaluateSolution(VVI & V, VI & partition );

    /**
     * For given graph and partition returns the numbers (insertions,deletions,total_modifications).
     */
    extern tuple<int,int,int> getEdgeModificationStatistics( VVI & V, VI & partition );
    
    /**
     * Calculates result for cluster given by state [st], based on [st.partition]
     */
    extern long long evaluateState( State & st );

    /**
     * For given partition into clusters, returns the VVI object containing the nodes in
     * the clusters
     * @param partition
     * @return
     */
    extern VVI partitionToClusters( VI & partition );

    /**
     * Remaps partition, such that there are no gaps in numeration.
     */
    extern VI properlyRemapPartition( VI part );

    /**
     * For a partition [part] of [clg] craetes and returns a partition of clg.origV
     */
    extern VI mapClgPartitionToOriginalPartition( ClusterGraph& clg, VI & part );

    /**
     * Given a partition [part] of [clg->origV], it creates a corresponding partition of [clg]
     */
    VI mapOriginalPartitionToClgPartition( ClusterGraph & clg, VI & part );

    /**
     * 'Reverse' of the [partitionToClusters] function.
     * @param clusters
     * @return
     */
    extern VI clustersToPartition( VVI & clusters );

    /**
     * From given graph in [V] and partition to clusters in [inCl], creates a cluster graph, where
     * every nodes represents a cluster,
     * @param V
     * @param inCl
     * @return
     */
    extern ClusterGraph convertToClusterGraph( VVI & V, VI & inCl );

    /**
     * Creates all paths of length 2 that occur as induced subgraph in [V]
     * @param V
     * @return
     */
    extern vector< tuple<int,int,int> > getInducedP2Paths(VVI & V);

    /**
     * Runs tests
     */
    extern void test();

    /**
     * Reads test in DIMACS format and outputs graph in VVI.
     */
    void convertDIMACSTestToVVI();
}

#endif //PACE21_PACEUTILS_H
