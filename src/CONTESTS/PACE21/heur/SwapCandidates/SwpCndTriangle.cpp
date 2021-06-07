//
// Created by sylwester on 3/20/21.
//

#include <graphs/GraphUtils.h>
#include <CONTESTS/PACE21/heur/Global.h>
#include "CONTESTS/PACE21/heur/SwapCandidates/SwpCndTriangle.h"



SwpCndTriangle::SwpCndTriangle( int swpval, int u, int v, int w, int trg_cl ){
    swap_value = swpval;
    nodes = VI({u,v,w});
    move_node_to = {trg_cl, trg_cl, trg_cl};
}

SwpCndTriangleCreator::SwpCndTriangleCreator(State &s) : SwapCandidateCreatorAdapter(s) {}

vector<SwpCndTriangle> SwpCndTriangleCreator::createSwapCandidates() {
    return vector<SwpCndTriangle>();
}

vector<SwapCandidate *> SwpCndTriangleCreator::createSwapCandidatesRaw() {
    return vector<SwapCandidate *>();
}

vector<SwpCndTriangle> SwpCndTriangleCreator::create_MoveTo_SwapCandidates(const bool only_empty_cluster) {
    const bool debug = false;

    vector<SwpCndTriangle> res;

    VVPII V = clg->V; // making a copy
    int N = (int)V.size();

    VI order(N);
    iota(ALL(order),0);
    sort(ALL(order), [&]( int a, int b ){
        if(V[a].size() != V[b].size()) return V[a].size() < V[b].size();
        else return a < b;
    } );
    VI inOrder(N);
    for(int i=0; i<N; i++) inOrder[order[i]] = i;

    int sq = (int)ceil( sqrt( 2 * GraphUtils::countEdges(V) ) ); // just to check, can be removed later
    for( int i=0; i<N; i++ ){ // removing all edges (a,b) such that b is before a in order
        for( int k = (int)V[i].size() - 1; k>=0; k-- ){
            if( inOrder[V[i][k].first] < inOrder[i] ){
                swap( V[i][k], V[i].back() );
                V[i].pop_back();
            }
        }
        sort( ALL(V[i]), [&](PII a, PII b){
            if(inOrder[a.first] != inOrder[b.first]) return inOrder[a.first] < inOrder[b.first];
            else return a < b;
        } );

        assert( V[i].size() <= sq );
    }

    if(debug) DEBUG(V);

    fill(ALL(was),false);
    fill(ALL(was2),false);
    VI weight_ac(N,0); // weight of edge (a,c)

    if(state->cl_neigh_graph.empty() ) state->createClNeighGraph(); // creating cl_neigh_graph if not created earlier

    /**
     * This function shall be called only if we consider other clusters to move the  triangle to. Otherwise, if the
     * swap candidate triangle is supposed to be moved to empty cluster, we do not need to create neighborhoods, nor
     * update edges_to_cluster array.
     *
     * Note, that this function may increase running time significantly, probably by   O( sqrt(E) )
     *
     * CAUTION! We do not unmark array was2 !!! It needs to be unmarked only by exclusive neighbors of nodes b or c,
     * manually.
     */
    auto markNeighbors = [&]( int x, VI & neigh, const bool pos ) {
        // considering all neighboring clusters to move the triangle to, so we create neighbors.
        for (auto &[cl, w] : state->cl_neigh_graph[x]) {

            if(pos) edges_to_cluster[cl] += w; // we increase even if in_cl_p == in_cl_a - it does not matter
            else edges_to_cluster[cl] -= w;

            if( pos && !was2[cl]) {
                was2[cl] = true;
                neigh.push_back(cl);
            }
        }
    };

    if(debug) DEBUG(order);

    const int empty_cluster_id = state->getIdOfEmptyCluster();

    for(int a : order){
        if( Global::checkTle() ) continue;

        int in_cl_a = state->inCl[a];
        int deg_in_cl_a = state->degInCl[a];
        int nw_a = clg->node_weights[a];
        int clw_a = state->clusters[in_cl_a].cluster_weight;

        for( auto & [d,w] : V[a] ){ was[d] = true; weight_ac[d] = w; }

        // marking edges to cluster. Now we do not care, whether the node p belongs to the same cluster as
        VI neigh_a; // cluster neighbors of node a, INCLUDING CLUSTER in_cl_a

        if(!only_empty_cluster) markNeighbors( a, neigh_a, true );

        if(debug){
            ENDL(10);
            DEBUG(a);
            DEBUG3(in_cl_a, deg_in_cl_a, nw_a);
            DEBUG2(clw_a, neigh_a);
            DEBUG(edges_to_cluster);
        }

        for( auto & [b,wab] : V[a] ){
            if( Global::checkTle() ) continue;

            int in_cl_b = state->inCl[b];
            int deg_in_cl_b = state->degInCl[b];
            int nw_b = clg->node_weights[b];
            int clw_b = state->clusters[in_cl_b].cluster_weight;

            VI neigh_b; // neigh_b is the 'set' of neighboring cluster of node b, that are NOT in neigh_a
            if(!only_empty_cluster) markNeighbors( b, neigh_b, true );

            if(debug){
                ENDL(6);
                DEBUG(b);
                DEBUG3(in_cl_b, deg_in_cl_b, nw_b);
                DEBUG3(clw_b, wab, neigh_b);
                DEBUG(edges_to_cluster);
            }

            for(auto & [c,wbc] : V[b]){
                if(!was[c]) continue; // #TEST #CAUTION if commented should still work, but may be much slower! original uncommented

                // now he have triangle a,b,c
                int wac = weight_ac[c];
                int in_cl_c = state->inCl[c];
                int deg_in_cl_c = state->degInCl[c];
                int nw_c = clg->node_weights[c];
                int clw_c = state->clusters[in_cl_c].cluster_weight;

                // neigh_c is the 'set' of neighboring cluster of node c, that are NEITHER in neigh_a, NOR in neigh_b
                VI neigh_c;
                if(!only_empty_cluster) markNeighbors( c, neigh_c, true );

                int nw_abc = nw_a + nw_b + nw_c;
                int wabc = wab + wbc + wac;

                if(debug){
                    ENDL(3);
                    DEBUG(c);
                    DEBUG3(in_cl_c, deg_in_cl_c, nw_c);
                    DEBUG2(clw_c, neigh_c);
                    DEBUG2(wbc, wac);
                    DEBUG2(wabc, nw_abc);
                    DEBUG(edges_to_cluster);
                    clog << "Triangle {" << a << "," << b << "," << c << "}" << endl;
                }
//                if(!debug) clog << "Triangle {" << a << "," << b << "," << c << "}" << endl;

                int swpval1 = 0;
                int deg_in_cl_abc = 0;

                if( in_cl_a == in_cl_b && in_cl_a == in_cl_c ){ // all nodes in the same cluster
                    // weight in clusters is just the sum of weights decreased by two times the wight of edges between
                    deg_in_cl_abc = deg_in_cl_a + deg_in_cl_b + deg_in_cl_c - ( wabc << 1 );
                    int clw = clw_a; // this also clw_a = clw_b = clw_c
                    swpval1 = nw_abc * (clw - nw_abc) - deg_in_cl_abc;
                }else if( in_cl_a != in_cl_b && in_cl_a != in_cl_c && in_cl_b != in_cl_c ){ // three different clusters
                    // just the sum, all edges ab,bc,ca are between clusters
                    deg_in_cl_abc = deg_in_cl_a + deg_in_cl_b + deg_in_cl_c;
                    swpval1 = nw_a * (clw_a - nw_a) +   nw_b * (clw_b - nw_b) +   nw_c * (clw_c - nw_c) - deg_in_cl_abc;
                }else if( in_cl_a == in_cl_b ){ // a and b in the same cluster
                    deg_in_cl_abc = deg_in_cl_c + (deg_in_cl_a + deg_in_cl_b - ( wab << 1 ) );
                    int nw_ab = nw_a + nw_b;
                    swpval1 = nw_c * (clw_c - nw_c) +   nw_ab * (clw_a - nw_ab) - deg_in_cl_abc;
                }else if( in_cl_a == in_cl_c ){ // a and c in the same cluster
                    deg_in_cl_abc = deg_in_cl_b + (deg_in_cl_a + deg_in_cl_c - ( wac << 1 ) );
                    int nw_ac = nw_a + nw_c;
                    swpval1 = nw_b * (clw_b - nw_b) +   nw_ac * (clw_a - nw_ac) - deg_in_cl_abc;
                } else if( in_cl_b == in_cl_c ){ // b and c in the same cluster
                    deg_in_cl_abc = deg_in_cl_a + (deg_in_cl_b + deg_in_cl_c - ( wbc << 1 ) );
                    int nw_bc = nw_b + nw_c;
                    swpval1 = nw_a * (clw_a - nw_a) +   nw_bc * (clw_b - nw_bc) - deg_in_cl_abc;
                }

                if(debug) DEBUG(swpval1);

                if( in_cl_a != in_cl_b ) swpval1 += wab - ( nw_a * nw_b - wab );
                if( in_cl_a != in_cl_c ) swpval1 += wac - ( nw_a * nw_c - wac );
                if( in_cl_b != in_cl_c ) swpval1 += wbc - ( nw_b * nw_c - wbc );

//                VI neigh_abc( 1, empty_cluster_id ); // original
                VI neigh_abc;
                if( !only_empty_cluster ) neigh_abc.reserve( neigh_a.size() + neigh_b.size() + neigh_c.size() );
                neigh_abc.push_back(empty_cluster_id);

                if( !only_empty_cluster ){
                    auto fun = [&]( int d ){ return d != in_cl_a && d != in_cl_b && d != in_cl_c; };
                    copy_if( ALL(neigh_a), back_inserter(neigh_abc), fun );
                    copy_if( ALL(neigh_b), back_inserter(neigh_abc), fun );
                    copy_if( ALL(neigh_c), back_inserter(neigh_abc), fun );
                }

                if(debug){
                    DEBUG3( swpval1, deg_in_cl_abc, neigh_abc );
                }

                createAllCandidatesForSwpCnd(swpval1, nw_abc, deg_in_cl_abc, neigh_abc, res,
                                             [aa=a,bb=b,cc=c]( int swpval, int trg_cl ){
                                                  return SwpCndTriangle( swpval, aa,bb,cc, trg_cl );
                                              } );

//                if(debug) DEBUG(res);

                if(!only_empty_cluster){
                    markNeighbors( c, neigh_c, false );
                    for( int d : neigh_c ) was2[d] = false; // clearing those clusters, to which only node c is incident
                }
            }

            if(!only_empty_cluster){
                markNeighbors( b, neigh_b, false );
                for( int d : neigh_b ) was2[d] = false; // clearing those clusters, to which only node b is incident
            }
        }

        if(!only_empty_cluster){
            markNeighbors( a, neigh_a, false );
            for( int d : neigh_a ) was2[d] = false; // clearing those clusters, to which node a is incident
        }

        for( auto & [d,w] : V[a] ){ was[d] = false; weight_ac[d] = 0; /* clearing*/ }
    }

    // if there is a swap candidate that contains whole cluster and is moved with swap value 0 to an empty cluster, then
    // we remove it from the set of swap candidates
    for( int i=(int)res.size()-1; i>=0; i-- ) {
        VI nd = res[i].getNodes();
        if( state->inCl[nd[0]] != state->inCl[nd[1]] ) continue;
        if( state->inCl[nd[0]] != state->inCl[nd[2]] ) continue;
        if( state->inCl[nd[1]] != state->inCl[nd[2]] ) continue;
        // now, all nodes in the same cluster
        if( state->clusters[ state->inCl[nd[0]] ].size() == 3 && res[i].move_node_to == VI(3, empty_cluster_id ) ){
            swap( res[i], res.back() );
            res.pop_back();
        }
    }

    return res;
}

ostream& operator<<(ostream& str, SwpCndTriangle& cnd){
    str << "(nodes --> move_to: " << cnd.nodes << "  -->  " << cnd.move_node_to
        << ", swp_val: " << cnd.swpVal() << ")";
    return str;
}

