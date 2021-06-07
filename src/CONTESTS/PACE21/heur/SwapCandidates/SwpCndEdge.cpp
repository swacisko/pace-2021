//
// Created by sylwester on 3/20/21.
//

#include <graphs/GraphUtils.h>
#include <datastructures/Heap.h>
#include <Constants.h>
#include "CONTESTS/PACE21/heur/SwapCandidates/SwpCndEdge.h"
#include "CollectionOperators.h"

SwpCndEdgeCreator::SwpCndEdgeCreator(State &s) : SwapCandidateCreatorAdapter(s) {}

SwpCndEdge::SwpCndEdge( int swpval, int u, int v, int trg_cl ){
    swap_value = swpval;
//    nodes = {u,v};
//    move_node_to = {trg_cl, trg_cl};
    this->u = u;
    this->v = v;
    move_node_to = trg_cl;
}


vector<SwapCandidate *> SwpCndEdgeCreator::createSwapCandidatesRaw() {
    return vector<SwapCandidate *>();
}

vector<SwpCndEdge> SwpCndEdgeCreator::createSwapCandidates() {
    return vector<SwpCndEdge>();
}

vector<SwpCndEdge> SwpCndEdgeCreator::create_MoveTo_SwapCandidatesForCluster(Cluster &cl, const bool only_common_neighbors) {
    const bool debug = false;

    vector<SwpCndEdge> res;
    int N = clg->N;
    VI ord(ALL(cl.g.nodes));
    sort(ALL(ord), [&](int a, int b){ return clg->V[a].size() > clg->V[b].size(); } ); // sorting in non-ascending order

    for( int d : ord ){
        int in_cl_d = state->inCl[d];
        int deg_in_cl_d = state->degInCl[d];
        int nw_d = clg->node_weights[d];
        int clw = cl.cluster_weight;
        was[d] = true;

        if(debug) clog << "Considering node d: " << d << endl;

        VI neigh_d;
        updateEdgesInNeighboringClustersForNode(d, in_cl_d, neigh_d);
        for( int c : neigh_d) was2[c] = true; // cluster c was visited as a neighbor of d

        if(debug){ DEBUG(neigh_d); DEBUG(edges_to_cluster); }

        for( auto& [p,w] : clg->V[d] ){
            if( was[p] ) continue; // pair (p,d) was already considered, no need to do it again
            if( state->inCl[p] != in_cl_d ) continue; // p and d need to be in the same cluster

            if(debug){
                clog << "unvisited neighbor p: " << p << endl;
                clog << "before update:"; DEBUG(edges_to_cluster);
            }

            int deg_in_cl_p = state->degInCl[p];
            int nw_p = clg->node_weights[p];
            int deg_in_cl_dp = deg_in_cl_d + deg_in_cl_p - (w<<1); // number of edges between p and d to the rest of cluster
            int nw_dp = nw_d + nw_p;

            VI neigh_dp( 1, state->getIdOfEmptyCluster() ); // clusters to which both d and p have edges
            if( !only_common_neighbors ) neigh_dp += neigh_d;

            auto mark = [&, p=p](const bool pos){
                for( auto& [q,w2] : clg->V[p] ){
                    int in_cl_q = state->inCl[q];
                    if( in_cl_q == in_cl_d ) continue;

                    if(only_common_neighbors) {
                        if (was2[in_cl_q]) { // q needs to be a neighbor of both d and p
                            if (pos) edges_to_cluster[in_cl_q] += w2;
                            else edges_to_cluster[in_cl_q] -= w2;

                            if (pos && !was3[in_cl_q]) {
                                was3[in_cl_q] = true;
                                neigh_dp.push_back(in_cl_q);
                            }
                        }
                    }else{
                        if (pos) edges_to_cluster[in_cl_q] += w2;
                        else edges_to_cluster[in_cl_q] -= w2;

                        if (pos && !was3[in_cl_q] && !was2[in_cl_q]) {
                            was3[in_cl_q] = true;
                            neigh_dp.push_back(in_cl_q);
                        }
                    }
                }
            };

            mark(true);

            if(debug){ DEBUG(neigh_dp);clog << "After update for p: "; DEBUG(edges_to_cluster); }

            if( find_possible_movements_to_all_clusters ){
                clog << "SwpCndEdge::create_MoveTo_SwapCandidatesForCluster, NOT IMPLEMENTED YET" << endl;
                exit(1);
            }else{ // creating candidates
                int swpval1 = nw_dp * ( clw - nw_dp ) - deg_in_cl_dp;
                createAllCandidatesForSwpCnd( swpval1, nw_dp, deg_in_cl_dp, neigh_dp, res,
                    [&d,pp=p]( int swpval, int trg_cl ){
                        return SwpCndEdge( swpval, d,pp, trg_cl );
                    } );
            }

            { // clearing section
                mark(false);
                for( auto a : neigh_dp ) was3[a] = false;
                neigh_dp.clear();
            }

            if(debug) ENDL(1);
        }

        for(int a : neigh_d) edges_to_cluster[a] = 0; // clearing array
        for( int c : neigh_d) was2[c] = false; // clearing

        if(debug) ENDL(3);
    }

    for( int d : ord ) was[d] = false;

    // if there is a swap candidate that contains whole cluster and is moved with swap value 0 to an empty cluster, then
    // we remove it from the set of swap candidates
    for( int i=(int)res.size()-1; i>=0; i-- ) {
        VI nd = res[i].getNodes();
        if( state->inCl[nd[0]] != state->inCl[nd[1]] ) continue;
        // now, all nodes in the same cluster
        if( state->clusters[ state->inCl[nd[0]] ].size() == 2
            &&
            VI(2,res[i].move_node_to) == VI(2, state->getIdOfEmptyCluster()) ){
            swap( res[i], res.back() );
            res.pop_back();
        }
    }

    return res;
}

vector<SwpCndEdge> SwpCndEdgeCreator::create_MoveTo_SwapCandidates_DifferentClusters(const bool only_common_neighbors) {
    const bool debug = false;

    vector<SwpCndEdge> res;
    int N = clg->N;

    VVPII Vout(N); // graph that contains only edges that go out of its own cluster
    for( int i=0; i<N; i++ ){
        int in_cl_i =  state->inCl[i];
        for( auto & [p,w] : clg->V[i] ){
            if( in_cl_i == state->inCl[p] ) continue;
            Vout[i].emplace_back(p,w);
        }
    }

    VI ord( N );
    iota(ALL(ord),0);
    sort(ALL(ord), [&](int a, int b){ return Vout[a].size() > Vout[b].size(); } ); // sorting in non-ascending order

    for( int d : ord ){
        int in_cl_d = state->inCl[d];
        int deg_in_cl_d = state->degInCl[d];
        int nw_d = clg->node_weights[d];
        int clw_d = state->clusters[in_cl_d].cluster_weight;
        was[d] = true;

        if(debug) clog << "Considering node d: " << d << endl;

        VI neigh_d;
        updateEdgesInNeighboringClustersForNode(d, in_cl_d, neigh_d);
        for( int c : neigh_d) was2[c] = true; // cluster c was visited as a neighbor of d

        if(debug){ DEBUG(neigh_d); DEBUG(edges_to_cluster); }

        for( auto& [p,w] : Vout[d] ){
            if( was[p] ) continue; // pair (p,d) was already considered, no need to do it again

            // p and d are certainly in different clusters, because we iterate over Vout structure.

            if(debug){
                clog << "unvisited neighbor p: " << p << endl;
                clog << "before update:"; DEBUG(edges_to_cluster);
            }

            int in_cl_p = state->inCl[p];
            int deg_in_cl_p = state->degInCl[p];
            int nw_p = clg->node_weights[p];
            int deg_in_cl_dp = deg_in_cl_d + deg_in_cl_p; // number of edges between p and d to the rest of cluster
            int clw_p = state->clusters[in_cl_p].cluster_weight;
            int nw_dp = nw_d + nw_p;

            VI neigh_dp( 1, state->getIdOfEmptyCluster() ); // clusters to which both d and p have edges
            if( !only_common_neighbors ){
                for( int x : neigh_d ) if( state->inCl[x] != in_cl_p ) neigh_dp.push_back(x); // do not add cluster in_cl_p
            }

            auto mark = [&, p=p](const bool pos){
                for( auto& [q,w2] : clg->V[p] ){
                    int in_cl_q = state->inCl[q];
                    if( in_cl_q == in_cl_d || in_cl_q == in_cl_p ) continue; // q must not be in d nor p in cluster

                    if(only_common_neighbors) {
                        if (was2[in_cl_q]) { // q needs to be a neighbor of both d and p
                            if (pos) edges_to_cluster[in_cl_q] += w2;
                            else edges_to_cluster[in_cl_q] -= w2;

                            if (pos && !was3[in_cl_q]) {
                                was3[in_cl_q] = true;
                                neigh_dp.push_back(in_cl_q);
                            }
                        }
                    }else{
                        if (pos) edges_to_cluster[in_cl_q] += w2;
                        else edges_to_cluster[in_cl_q] -= w2;

                        if (pos && !was3[in_cl_q] && !was2[in_cl_q]) {
                            was3[in_cl_q] = true;
                            neigh_dp.push_back(in_cl_q);
                        }
                    }
                }
            };

            mark(true);

            if(debug){ DEBUG(neigh_dp);clog << "After update for p: "; DEBUG(edges_to_cluster); }

            if( find_possible_movements_to_all_clusters ){
                clog << "SwpCndEdge::create_MoveTo_SwapCandidatesForCluster, NOT IMPLEMENTED YET" << endl;
                exit(1);
            }else{ // creating candidates
                int swpval1 = nw_d * (clw_d - nw_d) + nw_p * (clw_p - nw_p) - deg_in_cl_d - deg_in_cl_p;

                // ***** !!!! ******
                // to the final swap value, we need to add the number of missing edges between d and p. We can do that
                // by increasing swpval1 by the following value
                swpval1 += w - ( nw_d * nw_p - w);
                // ***** !!!! ******

                createAllCandidatesForSwpCnd( swpval1, nw_dp, deg_in_cl_dp, neigh_dp, res,
                                              [&d,pp=p]( int swpval, int trg_cl ){
                                                  return SwpCndEdge( swpval, d,pp, trg_cl );
                                              } );
            }

            { // clearing section
                mark(false);
                for( auto a : neigh_dp ) was3[a] = false;
                neigh_dp.clear();
            }

            if(debug) ENDL(1);
        }

        for(int a : neigh_d) edges_to_cluster[a] = 0; // clearing array
        for( int c : neigh_d) was2[c] = false; // clearing

        if(debug) ENDL(3);
    }

    for( int d : ord ) was[d] = false;

    return res;
}


ostream& operator<<(ostream& str, SwpCndEdge& cnd){
//    str << "(nodes --> move_to: " << cnd.nodes << "  -->  " << cnd.move_node_to
    str << "(nodes --> move_to: " << PII(cnd.u, cnd.v) << "  -->  " << cnd.move_node_to
        << ", swp_val: " << cnd.swpVal() << ")";
    return str;
}


