//
// Created by sylwester on 4/22/21.
//

#include <combinatorics/CombinatoricUtils.h>
#include <CONTESTS/PACE21/heur/PaceUtils.h>
#include <CONTESTS/PACE21/heur/ConvexHullTrickDynamic.h>
#include <CONTESTS/PACE21/heur/Global.h>
#include <CONTESTS/PACE21/heur/SwapCandidates/ComponentExpansionRepulsion.h>
#include <graphs/GraphUtils.h>
#include "StandardUtils.h"
#include "CONTESTS/PACE21/heur/StateImprovers/NodeEdgeGreedyNomap.h"

NodeEdgeGreedyNomap::NodeEdgeGreedyNomap(State &st) : NEG(st){
    initializeForState(st);
}

void NodeEdgeGreedyNomap::initializeForState(State &st) {

    {
        edges_to_cluster = VVPII(N);

        helper_etocl = VI(2*N,0);
        helper_was_etocl = VB(2*N,false);
        helper_was4 = VB(2*N,false);

        VI v; localShuffle(v); // initializing shuffle_seq

        initializeIndependentData(st);
        return;
    }

    edges_to_cluster = VVPII(N);

    clg = st.clg;
    N = clg->V.size();
    cluster_weights = VI( st.clusters.size() + 1 );
    for( auto & cl : st.clusters ){
        for(int d : cl.g.nodes) cluster_weights[cl.id] += clg->node_weights[d];
    }

    first_free_cluster.clear();
    for( int i=st.getIdOfEmptyCluster(); i <= 3*st.getIdOfEmptyCluster() + 10; i++ ) first_free_cluster.insert(i);


    inCl = st.inCl;
    degInCl = st.degInCl;

    nonempty_clusters_cnt = countNonemptyClusters();

    move_frequency = max(1.0, ceil(sqrt(N)));

    VI partition(clg->origV->size(),-1);
    for( int i=0; i<st.clusters.size(); i++ ){
        for( int d : st.clusters[i].g.nodes ){
            for( int p : clg->clusterNodes[d] ) partition[p] = i;
        }
    }

    current_result = PaceUtils::evaluateSolution( *clg->origV, partition );
    best_result = current_result;
    best_partition = PaceUtils::properlyRemapPartition(inCl);

    helper = VI(2*N,0);
    helper2 = VI(2*N,0);
    helper_was = VB(2*N,0);
    helper_was2 = VB(2*N,0);
    helper_was3 = VB(2*N,0);

    min_w2_to_cluster_of_given_weight = VI( clg->origV->size()+1, 0);

    queue.resize(2*N); queue.clear(); // with the hope of getting a single array underneath

    {
        edges_to_cluster_triangle = weight_ac_triangle = VI(2*N,0);
    }

    helper_etocl = VI(2*N,0);
    helper_was_etocl = VB(2*N,false);
}


void NodeEdgeGreedyNomap::improve() {
    if(!Global::disable_all_logs) clog << "Nomap improve" << endl;
    use_queue_propagation = false;


    // loop'. Node interchanging can be done just like chain2Swaps, outside the 'queue loop'
    NEG::improve();
}


tuple<int,int,int> NodeEdgeGreedyNomap::getBestNodeMoveForRange(VI &perm, int a, int b) {
    if(!use_node_swaps) return {-1,-1,1e9};
    const bool debug = false;

    if(debug){
        clog << "Checking moves of nodes from range perm[" << a << ":" << b << "] = "
             << StandardUtils::getSubarray(perm,a,b-1) << endl;
    }

    int best_swpval = 1e9;
    int best_v = -1;
    int best_to = -1;

    best_node_move_results.clear();

    for( int i=a; i<=b; i++ ){
        int d = perm[i];
        int cl_d = inCl[d];
        int nw_d = clg->node_weights[d];

        int w_0 = degInCl[d];

        int tot_cld_possible_edges = (cluster_weights[cl_d] - clg->node_weights[d]) * clg->node_weights[d];

       /* createEdgesToCluster(d);
        for( auto & [c,w] : edges_to_cluster[d] ){
            if( c == cl_d ) continue;
            int swpval = swapValueForNode(d, c, w_0, tot_cld_possible_edges, w);

            if(swpval < best_swpval){
                best_swpval = swpval;
                best_v = d;
                best_to = c;

                best_node_move_results.clear(); // #TEST
            }

            if( swpval <= best_swpval ) best_node_move_results.emplace_back( d,c,swpval ); // #TEST
        }*/

        { // alternative way, without creating edges_to_cluster
            for( auto & [d2,w2] : clg->V[d] ){
                int c = inCl[d2];
                if( c == cl_d ) continue;

                helper_etocl[c] += w2;
                int swpval = swapValueForNode(d, c, w_0, tot_cld_possible_edges, helper_etocl[c]);

                if(swpval < best_swpval){
                    best_swpval = swpval;
                    best_v = d;
                    best_to = c;

                    best_node_move_results.clear(); // #TEST
                }

                if( swpval <= best_swpval ) best_node_move_results.emplace_back( d,c,swpval ); // #TEST
            }

            for( auto & [d2,w2] : clg->V[d] ) helper_etocl[ inCl[d2] ] = 0;
        }


        if( cluster_weights[cl_d] != nw_d ){
            // if there is more than one node in cluster contatinig d, then we try to move it to an empty cluster
            int c = *first_free_cluster.begin();
            int w = 0;
//            while(cluster_weights.size() <= c) cluster_weights.push_back(0);
            resizeStructuresForEmptyCluster(c);

            int swpval = swapValueForNode(d, c, w_0, tot_cld_possible_edges, w);

            if(swpval < best_swpval){
                best_swpval = swpval;
                best_v = d;
                best_to = c;

                best_node_move_results.clear(); // #TEST
            }

            if( swpval <= best_swpval ) best_node_move_results.emplace_back( d,c,swpval ); // #TEST
        }
    }

    const bool USE_RANDOM_BEST_MOVE = true; // #TEST - getting random one of all best possible moves // originally false
    if(USE_RANDOM_BEST_MOVE && !best_node_move_results.empty() ){
        auto get_clw_sum = [&](tuple<int,int,int > & a){
            int to_a = get<1>(a);
            return cluster_weights[to_a] + cluster_weights[to_a];
        };

        if( prefer_cluster_mode != 0 ){
            localShuffle(best_node_move_results);
            sort(ALL(best_node_move_results), [&]( auto& a, auto& b ){ return get_clw_sum(a) < get_clw_sum(b); });
            if( prefer_cluster_mode == 2 ) reverse(ALL(best_node_move_results));

            for( int i=1; i<best_node_move_results.size(); i++ ){
                if( get_clw_sum( best_node_move_results[i] ) != get_clw_sum(best_node_move_results[i-1]) ){
                    best_node_move_results.resize(i);
                    break;
                }
            }
        }

        int ind = shuffle_seq[shuffle_ind++] % best_node_move_results.size();
        if(shuffle_ind == shuffle_seq.size()) shuffle_ind = 0;
        return best_node_move_results[ind];
    }

    return { best_v, best_to, best_swpval };
}

int NodeEdgeGreedyNomap::swapValueForNode(int v, int trg_cl, int edges_clv, int tot_clv_possible_edges, int edges_trg) {
    return cluster_weights[trg_cl] * clg->node_weights[v] - tot_clv_possible_edges - ((edges_trg - edges_clv)<<1);
}

tuple<PII, PII, int> NodeEdgeGreedyNomap::getBestEdgeMoveForRange(VI &perm, int a, int b) {
    if(!use_edge_swaps) return {{-1,-1},{-1,-1},(int)1e9};

    const bool debug = false;

    if(debug){
        clog << "************Trying edge swaps, perm[" << a << ":" << b << "]: " << StandardUtils::getSubarray(perm,a,b) << endl;
        DEBUG(current_result);
        DEBUG(inCl);
        DEBUG(PaceUtils::partitionToClusters(inCl));
        DEBUG(cluster_weights);
        DEBUG(first_free_cluster);
    }

    int best_swpval = 1e9;
    PII best_e = {-1,-1};
    PII best_to = {-1,-1};

    // this 2* factor is just to increase variability whilst preserving good complexity
    // selecting FACTOR 1.5 or 2 should be ok
    const double FACTOR = 1.5; // original 1.5
    const double ADD = 5;  // original 0


    const int empty_cluster = *first_free_cluster.begin();

    auto itt = first_free_cluster.begin(); itt++;
    const int empty_cluster_2 = *itt;

    assert(is_empty_cluster[empty_cluster]);
    assert(is_empty_cluster[empty_cluster_2]);

    resizeStructuresForEmptyCluster(empty_cluster_2);

    best_edge_move_results.clear();

    for( int i=a; i<=b; i++ ) {
        int u = perm[i];
        int cl_u = inCl[u];
        int nw_u = clg->node_weights[u];
        int etc_u_clu = degInCl[u];

        createEdgesToCluster(u);

        /**
         * best_u_moves is a pair (cl,val) of at most 2 best moves of node u to cluster cl with swap value val.
         * It is used for edge repulsion.
         */
        VPII best_u_moves;
        const bool USE_EDGE_REPULSION = use_edge_repulsion;
        if(USE_EDGE_REPULSION){
            best_u_moves.emplace_back( empty_cluster,
            swapValueForNode(u, empty_cluster, etc_u_clu, (cluster_weights[cl_u] - nw_u) * nw_u, 0) );
        }

        ConvexHullTrickDynamic hull(false); // minimizing
        hull.insert_line(0,0, empty_cluster);
        for (auto & [c, w2] : edges_to_cluster[u]) {
            if( c == cl_u ) continue;
            helper[c] += w2; // helper is just edges_to_cluster of both nodes d and d2

            int cwc = cluster_weights[c];
            if( min_w2_to_cluster_of_given_weight[cwc] > -(w2 << 1) ) { // #TEST - originally was >=
                hull.insert_line(cwc, -(w2 << 1), c);
                min_w2_to_cluster_of_given_weight[cwc] = -(w2 << 1);
            }

            if(USE_EDGE_REPULSION){
                int swpval = swapValueForNode( u, c, etc_u_clu, (cluster_weights[cl_u] - nw_u) * nw_u, w2  );
                if( best_u_moves.size() == 2 ){
                    if( best_u_moves[0].second > best_u_moves[1].second ){
                        swap( best_u_moves[0], best_u_moves[1] );
                    }
                    if( best_u_moves[1].second > swpval ) best_u_moves[1] = {c,swpval};
                }else{
                    best_u_moves.emplace_back(c,swpval);
                }
            }

            if(debug) clog << "Adding line for cluster " << c << ": " << cluster_weights[c] << " * x - " << (w2<<1) << endl;
        }
        for (auto & [c, w2] : edges_to_cluster[u]) min_w2_to_cluster_of_given_weight[cluster_weights[c]] = 0; // clearing


        for( auto & [v,w0] : clg->V[u] ){

            if( clg->V[v].size() > FACTOR*clg->V[u].size() + ADD ) continue;

            int cl_v = inCl[v];
            int nw_v = clg->node_weights[v];

            int residue = 0, s = nw_u + nw_v;

            int etc_v_clv = degInCl[v];

            if( cl_v == cl_u ){
                int w1 = etc_u_clu + etc_v_clv - (w0<<1);
                residue = (w1<<1) - s * ( cluster_weights[cl_u] - s );
            }else{
                int w1 = etc_u_clu + etc_v_clv;
                residue = nw_u * ( nw_u - cluster_weights[cl_u] ) + nw_v * ( nw_v - cluster_weights[cl_v] )
                        + nw_u * nw_v - 2*( w0-w1 );
            }


            int best_c, best_c_val;
            best_c = hull.get_line_id(s);
            best_c_val = 1e9 + 1;

//            while(cluster_weights.size() <= best_c) cluster_weights.push_back(0);
            resizeStructuresForEmptyCluster(best_c);
            if(debug){
                ENDL(1);
                DEBUG(v); DEBUG(w0); DEBUG(cl_v); DEBUG(nw_v); DEBUG(residue);
            }

            if(best_c != cl_u && best_c != cl_v){
                int e_to_bestc = 0;
                if( best_c != empty_cluster ) e_to_bestc = findEdgesToCluster(u,best_c);

                best_c_val = cluster_weights[best_c] * s - ( e_to_bestc << 1 )  + residue;
            }
            else if(debug){
                clog << "Trying to move edge to cluster containing one of its ends" << endl;
                DEBUG(best_c); DEBUG(cl_u); DEBUG(cl_v);
            }


            if( cl_u == cl_v && nw_u + nw_v == cluster_weights[cl_u] && cluster_weights[best_c] == 0 ){
                // if we would like to move a whole two-node cluster to an empty cluster
                best_c_val = 1e9;
            }

            if(debug){ DEBUG(best_c); DEBUG(best_c_val); }

            VPII best_v_moves;
            if (cl_v == cl_u && USE_EDGE_REPULSION) {
                best_v_moves.emplace_back(empty_cluster_2,
                     swapValueForNode(v, empty_cluster_2, etc_v_clv,
                                      (cluster_weights[cl_v] - nw_v) * nw_v, 0));
            }

            const int MODE = 0;  // #TEST - originally 0
            if(MODE == 0) { // alternative without creating edges_to_cluster
                for( auto & [d2,w3] : clg->V[v] ){
                    int c = inCl[d2];
                    if( c == cl_u || c == cl_v ) continue;

                    helper_etocl[c] += w3;
                    int w2 = helper[c] + helper_etocl[c];

                    int c_val = cluster_weights[c] * s - (w2<<1) + residue;

                    if( c_val < best_c_val ){
                        best_c_val = c_val;
                        best_c = c;
                    }

                    if( cl_v == cl_u && USE_EDGE_REPULSION){
                        int swpval2 =  swapValueForNode(v, c, etc_v_clv, (cluster_weights[cl_v] - nw_v) * nw_v, helper_etocl[c]);
                        if( best_v_moves.size() == 2 ){
                            if( best_v_moves[0].second > best_v_moves[1].second ){
                                swap( best_v_moves[0], best_v_moves[1] );
                            }
                            if( best_v_moves[1].second > swpval2 ) best_v_moves[1] = {c,swpval2};
                        }else{
                            best_v_moves.emplace_back(c,swpval2);
                        }
                    }
                }

                for( auto & [d2,w3] : clg->V[v] ) helper_etocl[ inCl[d2] ] = 0;
            }
            else{ // creating edges_to_cluster[v]
                createEdgesToCluster(v, false);
                for (auto[c, w2] : edges_to_cluster[v]) {
                    if (c == cl_u || c == cl_v) continue;
                    w2 += helper[c];
                    int c_val = cluster_weights[c] * s - (w2 << 1) + residue;
                    if (debug) { DEBUG(c);DEBUG(w2);DEBUG(c_val); }

                    if (c_val < best_c_val) { best_c_val = c_val;best_c = c; }

                    if (cl_v == cl_u && USE_EDGE_REPULSION) {
                        int swpval2 = swapValueForNode(v, c, etc_v_clv, (cluster_weights[cl_v] - nw_v) * nw_v,
                                                       w2 - helper[c]);
                        if (best_v_moves.size() == 2) {
                            if (best_v_moves[0].second > best_v_moves[1].second) {
                                swap(best_v_moves[0], best_v_moves[1]);
                            }
                            if (best_v_moves[1].second > swpval2) best_v_moves[1] = {c, swpval2};
                        } else {
                            best_v_moves.emplace_back(c, swpval2);
                        }
                    }
                }
            }

            if(cl_v == cl_u) {
                for (PII umove : best_u_moves) {
                    for (PII vmove : best_v_moves) {
                        if( umove.first == vmove.first ) continue;

                        int swpval1 = umove.second;
                        int swpval2 = vmove.second;
                        int swpval = swpval1 + swpval2 + nw_u * nw_v - (w0 << 1);
                        if (swpval < best_swpval) {
                            best_swpval = swpval;
                            best_e = {u, v};
                            best_to = {umove.first, vmove.first};
                            best_edge_move_results.clear();
                        }
                        if (swpval <= best_swpval) {
                            best_edge_move_results.emplace_back(PII(u, v), PII(umove.first, vmove.first), swpval);
                        }
                    }
                }
            }

            if(best_c_val < best_swpval){
                best_swpval = best_c_val;
                best_e = {u,v};
                best_to = {best_c,best_c};
                best_edge_move_results.clear();
            }
            if( best_c_val <= best_swpval ) best_edge_move_results.emplace_back( PII(u,v), PII(best_c,best_c), best_c_val );
        }

        for (auto & [c, w2] : edges_to_cluster[u]) helper[c] = 0; // clearing helper;

    }

    const bool USE_RANDOM_BEST_MOVE = true; // #TEST - getting random one of all best possible moves // originally false
    if(USE_RANDOM_BEST_MOVE && !best_edge_move_results.empty() ){
        auto get_clw_sum = [&](tuple< PII, PII,int > & a){
            PII to_a = get<1>(a);
            return cluster_weights[to_a.first] + cluster_weights[to_a.second];
        };

        if( prefer_cluster_mode != 0 ){
            localShuffle(best_edge_move_results);
            sort(ALL(best_edge_move_results), [&]( auto& a, auto& b ){ return get_clw_sum(a) < get_clw_sum(b); });
            if( prefer_cluster_mode == 2 ) reverse(ALL(best_edge_move_results));
            for( int i=1; i<best_edge_move_results.size(); i++ ){
                if( get_clw_sum( best_edge_move_results[i] ) != get_clw_sum(best_edge_move_results[i-1]) ){
                    best_edge_move_results.resize(i);
                    break;
                }
            }
        }

        int ind = shuffle_seq[shuffle_ind++] % best_edge_move_results.size();
        if(shuffle_ind == shuffle_seq.size()) shuffle_ind = 0;
        return best_edge_move_results[ind];
    }


    return {best_e, best_to, best_swpval};
}


void NodeEdgeGreedyNomap::moveNodeTo(int v, int to) {
    const bool debug = false;

    int cl_v = inCl[v];
    int nw_v = clg->node_weights[v];
    cluster_weights[cl_v] -= nw_v;

    if(cluster_weights[cl_v] == 0) nonempty_clusters_cnt--;

    if( is_empty_cluster[to] ){
        first_free_cluster.erase(to);
        is_empty_cluster[to] = false;
    }

    if(cluster_weights[cl_v] == 0){
        first_free_cluster.insert(cl_v);
        is_empty_cluster[cl_v] = true;
    }

    degInCl[v] = 0;
    for( auto & [p,w] : clg->V[v] ){

        if( inCl[p] == cl_v ) degInCl[p] -= w;

        if( inCl[p] == to ){
            degInCl[p] += w;
            degInCl[v] += w;
        }
    }

    inCl[v] = to;

//    while(cluster_weights.size() <= to) cluster_weights.push_back(0);
    resizeStructuresForEmptyCluster(to);
    cluster_weights[to] += nw_v;
    if( cluster_weights[to] == nw_v ) nonempty_clusters_cnt++;
}


SwapCandidateAdapter NodeEdgeGreedyNomap::getBestTriangleAll(int v) {
    const bool debug = false;
    const bool only_empty_cluster = (!use_triangle_swaps_to_other_clusters); // #TEST

    SwapCandidateAdapter res;
    res.swap_value = 1e9;
    res.move_node_to = {-1};
    res.nodes = {-1};

    VI to_clear;

    auto markNeighbors = [&]( int x, VI & neigh, const bool pos ) {
        if(!helper_was4[x]){
            helper_was4[x] = true;
            createEdgesToCluster(x,false);
            to_clear.push_back(x);
        }
        // considering all neighboring clusters to move the triangle to, so we create neighbors.
        for (auto &[cl, w] : edges_to_cluster[x]) {

            if(pos) edges_to_cluster_triangle[cl] += w; // we increase even if in_cl_p == in_cl_a - it does not matter
            else edges_to_cluster_triangle[cl] -= w;

            if( pos && !helper_was2[cl]) {
                helper_was2[cl] = true;
                neigh.push_back(cl);
            }
        }
    };

    const int MAX_BEST_CL_SIZE = max_best_cl_size_triangle_swaps;
    auto insertToBestClusters = [&]( int cl, int val, VPII & best_cl ){
        if(best_cl.size() == MAX_BEST_CL_SIZE ){
            auto it = max_element( ALL(best_cl), [](auto & a, auto & b){
                return a.second < b.second;
            } );
            auto ind = it - best_cl.begin();
            if(it->second > val) best_cl[ind] = {cl,val};
        }else{
            best_cl.emplace_back(cl,val);
        }
    };

    const int empty_cluster_id = *first_free_cluster.begin();
    resizeStructuresForEmptyCluster(empty_cluster_id);

    const double FACTOR = 1.33; // original 1.5
    const double ADD = 5;  // original 0

    VI neigh_a, neigh_b, neigh_c, neigh_abc;
    /**
     * a vector containing best cluster candidates for given nodes a and b. For those candidates node c will be checked
     * This should be kept as small as possible to keep good complexity
     */
    VPII best_clusters;
    VPII best_cl_a;


    for(int a : {v}){
        if( Global::checkTle() ) continue;

        int in_cl_a = inCl[a];
        int deg_in_cl_a = degInCl[a];

        int nw_a = clg->node_weights[a];
        int clw_a = cluster_weights[in_cl_a];

        for( auto & [d,w] : clg->V[a] ){ helper_was[d] = true; weight_ac_triangle[d] = w; }

        neigh_a.clear();

        if(!only_empty_cluster){
            neigh_a.clear();
            best_cl_a.clear();
            markNeighbors( a, neigh_a, true );

            for( auto & [cl,w] : edges_to_cluster[a] ){
                if( cl == in_cl_a ) continue;
                int val = swapValueForNode(a,cl,deg_in_cl_a, (clw_a - nw_a) * nw_a, edges_to_cluster_triangle[cl]);
                insertToBestClusters(cl,val,best_cl_a);
            }
        }

        for( auto & [b,wab] : clg->V[a] ){
            if( Global::checkTle() ) continue;
            if( clg->V[b].size() > FACTOR * clg->V[a].size() + ADD ) continue; // #TEST

            int in_cl_b = inCl[b];
            int deg_in_cl_b = degInCl[b];

            int nw_b = clg->node_weights[b];
            int clw_b = cluster_weights[in_cl_b];

            if(!only_empty_cluster){
                neigh_b.clear();
                best_clusters = best_cl_a;
                markNeighbors( b, neigh_b, true );

                for( auto & [cl,w] : edges_to_cluster[b] ){
                    if( (cl == in_cl_a) || (cl == in_cl_b) ) continue;
                    int val = swapValueForNode(b,cl,deg_in_cl_b, (clw_b - nw_b) * nw_b, edges_to_cluster_triangle[cl]);
                    insertToBestClusters(cl,val,best_clusters);
                }
            }

            for(auto & [c,wbc] : clg->V[b]){
                if( c == a ) continue;

                if( clg->V[c].size() > FACTOR * clg->V[a].size() + ADD  ||
                 clg->V[c].size() > FACTOR * clg->V[b].size() + ADD ){
                    continue;
                }

                // now he have triangle a,b,c or path P3 in the form a,b,c with no edge (a,c)
                int wac = weight_ac_triangle[c];

                int in_cl_c = inCl[c];
                int deg_in_cl_c = degInCl[c];

                int nw_c = clg->node_weights[c];
                int clw_c = cluster_weights[in_cl_c];

                int nw_abc = nw_a + nw_b + nw_c;
                int wabc = wab + wbc + wac;

                LL swpval1 = 0;
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

                if( in_cl_a != in_cl_b ) swpval1 += wab - ( nw_a * nw_b - wab );
                if( in_cl_a != in_cl_c ) swpval1 += wac - ( nw_a * nw_c - wac );
                if( in_cl_b != in_cl_c ) swpval1 += wbc - ( nw_b * nw_c - wbc );

                neigh_abc.clear();
                neigh_abc.push_back(empty_cluster_id);

                const bool ADD_CLUSTERS = (!only_empty_cluster) && helper_was[c] && (b<c); // we do not need to consider b>c, since (a,b,c) is a triangle
                if( ADD_CLUSTERS ) {
                    if(!helper_was4[c]){ // creating edges to cluster, if not created earlier
                        createEdgesToCluster(c);
                        helper_was4[c] = true; to_clear.push_back(c);
                    }

                    assert(best_clusters.size() <= MAX_BEST_CL_SIZE);
                    for (auto p : best_clusters) {
                        int cl = p.first;
                        if( (cl == in_cl_a) || (cl == in_cl_b) || (cl == in_cl_c) ) continue;

                        if(!helper_was3[cl]){
                            neigh_abc.push_back(cl); // adding best clusters to neigh_abc to check
                            helper_was3[cl] = true;
                        }
                    }

                    for( auto cl : neigh_abc ) edges_to_cluster_triangle[cl] += findEdgesToCluster(c, cl);
                }

                {
                    auto getSwpValForMove = [&](LL swp_val1, int swpcnd_w, int deg_in_cl, int swp_trg_cl) {
                        LL swp_val2 = swpcnd_w * cluster_weights[swp_trg_cl] - edges_to_cluster_triangle[swp_trg_cl] + deg_in_cl;
                        swp_val1 += edges_to_cluster_triangle[swp_trg_cl];
                        return swp_val2 - swp_val1;
                    };

                    for( int trg_cl : neigh_abc ){
                        LL final_swpval = getSwpValForMove( swpval1, nw_abc, deg_in_cl_abc,  trg_cl ); // swpval after moving d to cluster trg_cl

                        if( final_swpval < res.swpVal() ){
                            res.swap_value = final_swpval;
                            res.nodes = {a,b,c};
                            res.move_node_to = {trg_cl, trg_cl, trg_cl};
                        }
                    }
                }

                if( ADD_CLUSTERS ) {
                    for( auto cl : neigh_abc ){
                        edges_to_cluster_triangle[cl] -= findEdgesToCluster(c, cl);
                        helper_was3[cl] = false;
                    }
                }
            }

            if(!only_empty_cluster){
                markNeighbors( b, neigh_b, false );
                for( int d : neigh_b ) helper_was2[d] = false; // clearing those clusters, to which only node b is incident
            }
        }

        if(!only_empty_cluster){
            markNeighbors( a, neigh_a, false );
            for( int d : neigh_a ) helper_was2[d] = false; // clearing those clusters, to which node a is incident
        }

        for( auto & [d,w] : clg->V[a] ){ helper_was[d] = false; weight_ac_triangle[d] = 0; /* clearing*/ }
    }

    for( int x : to_clear ) helper_was4[x] = false;

    return res;
}

void NodeEdgeGreedyNomap::createEdgesToCluster(int v, bool use_sort) {
    edges_to_cluster[v].clear();

    VI neigh_cl; neigh_cl.reserve(clg->V[v].size());

    for( auto & [u,w] : clg->V[v] ){
        int cl_u = inCl[u];
        helper_etocl[cl_u] += w;
        if( !helper_was_etocl[cl_u] ){
            helper_was_etocl[cl_u] = true;
            neigh_cl.push_back(cl_u);
        }
    }

    for( int cl : neigh_cl ){
        edges_to_cluster[v].emplace_back( cl, helper_etocl[cl] );
        helper_was_etocl[cl] = helper_etocl[cl] = 0; // clearing in the same time
    }

    localShuffle(edges_to_cluster[v]); // #TEST - using local shuffle to increase randomness

    bool use_binary_search =
            use_sort &&
            use_binary_search_etocl &&
                    (edges_to_cluster[v].size() > MAGIC_SORT_THR); // #TEST - sorting - this takes time! but anables faster findEdgesToCluster queries using binary search
    if(use_binary_search){ // sorting edges_to_cluster[v] to enable binary search
        sort( ALL(edges_to_cluster[v]) );
    }
}

int NodeEdgeGreedyNomap::findEdgesToCluster(int v, int cl) {
    if( cl == inCl[v] ) return degInCl[v];
    if(is_empty_cluster[cl]) return 0;

    const bool use_binary_search =
            use_binary_search_etocl &&
            (edges_to_cluster[v].size() > MAGIC_SORT_THR);

    if(use_binary_search){
        auto it = lower_bound( ALL(edges_to_cluster[v]), PII(cl,-1) );
        if( it != edges_to_cluster[v].end() && it->first == cl ) return it->second;
        else return 0;
    }else{
        for( PII & p : edges_to_cluster[v] ){
            if( p.first == cl ){
                return p.second;
            }
        }
        return 0;
    }
}

VPII NodeEdgeGreedyNomap::getEdgesToCluster(int v) {
    return edges_to_cluster[v];
}

void NodeEdgeGreedyNomap::createClusterNodes() {
    int C = maxClusterId();
    cluster_nodes = VVI(C+1);
    for( int i=0; i<N; i++ ) cluster_nodes[ inCl[i] ].push_back(i);
}

VI NodeEdgeGreedyNomap::getClusterNodes(int c) {
    return cluster_nodes[c];
}

int NodeEdgeGreedyNomap::maxClusterId() {
    return *max_element(ALL(inCl));
}

template<class _T>
void NodeEdgeGreedyNomap::localShuffle(_T &v) {
    if( shuffle_seq.empty() ){
        shuffle_seq = VI(N);
        for( int i=0; i<N; i++ ) shuffle_seq[i] = rnd.nextInt(1e9);
    }

    for( int i=(int)v.size()-1; i>=0; i-- ){
        int ind = shuffle_seq[shuffle_ind] % (i+1);
        if( ind != i ) swap( v[i],v[ind] );
        shuffle_ind++;
        if(shuffle_ind == shuffle_seq.size()) shuffle_ind = 0;
    }
}
