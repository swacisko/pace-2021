//
// Created by sylwester on 3/9/21.
//

#include "CONTESTS/PACE21/heur/SwapCandidates/SwapCandidate.h"


vector<Cluster *> SwapCandidate::getAffectedClusters(VPII & toSwap, vector<Cluster> &clusters, VI& cl) {
    set<int> zb;
    for( PII p : toSwap ){
        zb.insert( cl[p.first] );
        zb.insert( p.second );
    }

    vector<Cluster*> res;
    res.reserve(zb.size());
    for( auto p : zb ) res.push_back( &clusters[p] );

    return res;
}

void SwapCandidate::test() {

    VVI V = {
            {1, 2},
            {0, 4},
            {0, 3, 4},
            {2},
            {1, 2, 5},
            {4, 6},
            {5}
    };

    int N = V.size();
    VI cl(N, 0);
    cl[0] = cl[1] = 0;
    cl[2] = cl[4] = 1;
    cl[3] = 2;
    cl[5] = cl[6] = 3;

    VPII toSwap = {
            {3,1}, {6,1}
    };

    ClusterGraph g(&V, cl);

    vector<Cluster> clusters;
    clusters.emplace_back( g, VI({0,1}),0 );
    clusters.emplace_back( g, VI({2,4}),1 );
    clusters.emplace_back( g, VI({1}),2 );
    clusters.emplace_back( g, VI({5,6}),3 );

    vector<Cluster*> affCl = getAffectedClusters(toSwap, clusters, cl);

    assert( affCl[0]->id == 1 );
    assert( affCl[1]->id == 2 );
    assert( affCl[2]->id == 3 );

    clog << "SwapCandidate getAffectedClusters test passed" << endl;
}

SwapCandidateCreatorAdapter::SwapCandidateCreatorAdapter(State &s) {
    state = &s;
    clg = state->clg;
    edges_to_cluster = VI(state->clusters.size()+1, 0);
    was = VB(2*s.N,false);
    was2 = VB(2*s.N,false);
    was3 = VB(2*s.N,false);
}

void SwapCandidateCreatorAdapter::updateEdgesInNeighboringClustersForNode(int d, int in_cl_d, VI &neigh) {
    if( state->cl_neigh_graph.empty() ) state->createClNeighGraph(); // creating structure if necessary

    { // using state->cl_neigh_graph structure for speedup
        for (auto & [cl, w] : state->cl_neigh_graph[d]) {
            if (in_cl_d == cl) continue; // we want p and d to be in different clusters
            if (edges_to_cluster[cl] == 0) neigh.push_back(cl);
            edges_to_cluster[cl] += w;
        }
    }
}

int SwapCandidateCreatorAdapter::getSwpValForMove(int swp_val1, int swpcnd_w, int deg_in_cl, int swp_trg_cl) {
    int swp_val2 = swpcnd_w * state->clusters[swp_trg_cl].cluster_weight - edges_to_cluster[swp_trg_cl] + deg_in_cl;
//    DEBUG(swp_val2);
    swp_val1 += edges_to_cluster[swp_trg_cl];
//    DEBUG(swp_val1);
    return swp_val2 - swp_val1;
}
