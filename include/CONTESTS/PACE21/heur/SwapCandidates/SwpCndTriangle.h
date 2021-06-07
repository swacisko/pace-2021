//
// Created by sylwester on 3/20/21.
//

#ifndef ALGORITHMSPROJECT_SWPCNDTRIANGLE_H
#define ALGORITHMSPROJECT_SWPCNDTRIANGLE_H

#include "SwapCandidate.h"

/**
 * Swap Candidate - moving triangles.
 */
class SwpCndTriangle : public SwapCandidateAdapter{
public:
    /**
     * Creates a swap candidate - moving three nodes u,v,w, all to cluster [trg_cl].
     */
    SwpCndTriangle( int swpval, int u, int v, int w, int trg_cl );

    friend ostream& operator<<(ostream& str, SwpCndTriangle& cnd);
};

/**
 * This class checks each triangle (a,b,c) in the cluster graph created for given state.
 * It creates (depending on parameters) all/nonpositive swap candidates from those triangles, by moving it to
 * some cluster from the neighborhood of {a,b,c} (or an empty cluster).
 *
 * For each found triangle checks some number of neighboring clusters. Hence, it works in time O(E * sqrt(E) * X),
 * where X is some value dependent on the maximum number of neighboring clusters of a triangle. Pessimistically it can
 * be O(N) if we consider checking all neighboring clusters.
 */
class SwpCndTriangleCreator : public SwapCandidateCreatorAdapter{
public:
    SwpCndTriangleCreator(State &s);

    vector<SwapCandidate *> createSwapCandidatesRaw() override;

    vector<SwpCndTriangle> createSwapCandidates();

    /**
     *
     * @param only_empty_cluster if true, then we consider only as swap candidates moving triangles to empty cluster.
     * This works in time O(E * sqrt(E)). Otherwise we consider moving triangles to any neighboring clusters. This
     * may work in O( E * sqrt(E) * state->clusters.size() ).
     * @return
     */
    vector<SwpCndTriangle> create_MoveTo_SwapCandidates(const bool only_empty_cluster = true );
};

#endif //ALGORITHMSPROJECT_SWPCNDTRIANGLE_H
