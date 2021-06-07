//
// Created by sylwester on 4/22/21.
//

#ifndef ALGORITHMSPROJECT_NODEEDGEGREEDYW1_H
#define ALGORITHMSPROJECT_NODEEDGEGREEDYW1_H

#include <CONTESTS/PACE21/heur/Config.h>
#include "../State.h"
#include "../SwapCandidates/SwapCandidate.h"
#include "NodeEdgeGreedyNomap.h"

/**
 * Algorithm works in iterations.
 * In each iteration a random permutation of nodes is chosen. Then, in that order, nodes are tried to be moved to some
 * other cluster. If, at any point, a node with NONPOSITIVE swap value is found, then it is instantaneously moved.
 *
 * The same goes for checking edges.
 *
 * If during whole iteration no improvement was found, algorithm stops.
 */
class NodeEdgeGreedyW1 : public NodeEdgeGreedyNomap{
public:

    NodeEdgeGreedyW1(State & st );

    virtual void initializeForState(State & st);
    virtual void improve() override;

//private:


    template<class _T> void localShuffle(_T & v);

    virtual void shuffleClg() override;
    VVI V;

    virtual tuple<int,int,int> getBestNodeMoveForRange( VI & perm, int a, int b ) override;
    virtual SwapCandidateAdapter getBestTriangleAll(int v) override;
    VB helper_was4; // used in getBestTriangleAll

    virtual void moveNodeTo(int v, int to) override;
    virtual tuple<PII, PII, int> getBestEdgeMoveForRange(VI & perm, int a, int b ) override;

    virtual int swapValueForNode(int v, int trg_cl, int edges_clv, int tot_clv_possible_edges, int edges_trg) override;
    virtual bool compareCurrentResultWithBruteResult() override;

    virtual void createEdgesToCluster(int v, bool use_sort = true) override;

};

#endif //ALGORITHMSPROJECT_NODEEDGEGREEDYW1_H
