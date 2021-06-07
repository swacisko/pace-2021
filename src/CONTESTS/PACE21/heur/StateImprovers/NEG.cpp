//
// Created by sylwester on 5/24/21.
//

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
#include "CONTESTS/PACE21/heur/StateImprovers/NEG.h"

NEG::NEG(State &st) : N(st.N), rnd(UniformIntGenerator(0,1e9)){
    clg = st.clg;
}


void NEG::initializeIndependentData(State & st){
    cluster_weights = VI( st.clusters.size() + 1 );
    for( auto & cl : st.clusters ){
        for(int d : cl.g.nodes) cluster_weights[cl.id] += clg->node_weights[d];
    }


    first_free_cluster.clear();
    int A = st.getIdOfEmptyCluster();
    int B = 3*(st.getIdOfEmptyCluster() + 10);
//    int B = st.clg->origV->size() + 10;
    for( int i=A; i <= B; i++ ) first_free_cluster.insert(i);
    is_empty_cluster = VB(B+1);
    for(int d : first_free_cluster) is_empty_cluster[d] = true;


    inCl = st.inCl;
    degInCl = st.degInCl;
//    move_frequency = max(1.0, ceil(sqrt(N)));
    move_frequency = 2; // #TEST

    VI partition(clg->origV->size(),-1);
    for( int i=0; i<st.clusters.size(); i++ ){
        for( int d : st.clusters[i].g.nodes ){
            for( int p : clg->clusterNodes[d] ) partition[p] = i;
        }
    }

    current_result = PaceUtils::evaluateSolution( *clg->origV, partition );
    best_result = current_result;
    initial_result = current_result;

    best_partition = PaceUtils::properlyRemapPartition(inCl);

    helper = VI(2*N,0);
    helper2 = VI(2*N,0);
    helper_was = VB(2*N,0);
    helper_was2 = VB(2*N,0);
    helper_was3 = VB(2*N,0);

    min_w2_to_cluster_of_given_weight = VI( clg->origV->size()+1, 0);

    checked_for_v = VB(N,false);

    queue.resize(2*N); queue.clear(); // with the hope of getting a single array underneath

    {
        edges_to_cluster_triangle = weight_ac_triangle = VI(2*N,0);
    }

    nonempty_clusters_cnt = countNonemptyClusters();
}

void NEG::shuffleClg(){
    for( int i=0; i<clg->N; i++ ){
        localShuffle( clg->V[i] );
    }
};

void NEG::improve() {
    const bool debug = false;

    VI v; localShuffle(v);

    auto debug_all = [&](){
        clog << endl;
        {
            DEBUG(inCl);
            map<int, VI> clusters;
            clog << "Clusters: " << endl;
            for (int i = 0; i < inCl.size(); i++) clusters[inCl[i]].push_back(i);
            for( auto & [id,nd] : clusters ) clog << id << " --> " << nd << endl;
            DEBUG(cluster_weights);
        }
        DEBUG(current_result);
    };

    UniformIntGenerator rnd(0,1e9);

    int nn_iter = 0;
    int iter = 0;
    iterations_done = &iter;
    perturbations_done = 0; //    int perturbations_done = 0;
    int last_nn_iter = (max_nonnegative_iters-1);
    int iters_since_last_perturbation = 0;

    int TRIANGLE_SWAPS_FREQUENCY = triangle_swaps_frequency;
    int EDGE_SWAPS_FREQUENCY = edge_swaps_frequency;
    int NODE_INTERCHANING_FREQUENCY = node_interchanging_frequency;
    int CHAIN2_FREQUENCY = chain2_swaps_frequency;
    int JOIN_CLUSTERS_FREQUENCY = join_clusters_frequency;

    {
        double E = GraphUtils::countEdges(clg->V);
        double avg_deg = 2*E / clg->V.size();

        int MIN_DENSITY = 30; // #TEST!! originally this scaling was present!! - original value 30
        while( avg_deg > MIN_DENSITY ){
            TRIANGLE_SWAPS_FREQUENCY *= 1.5;
            EDGE_SWAPS_FREQUENCY *= 1.5;
            NODE_INTERCHANING_FREQUENCY *= 1.5;
            avg_deg /= 2;
        }

//        DEBUG(EDGE_SWAPS_FREQUENCY);
//        DEBUG(TRIANGLE_SWAPS_FREQUENCY);

    }

    VI iter_nodes_to_check;
    VI perm = CombinatoricUtils::getRandomPermutation(N, rnd.rand());

    queue.clear();
    in_queue = VB(N,false);

    while( perturbations_done <= max_perturbations ){
        if(iter == max_iterations_to_do) break; // terminate early

        iter++;
        iters_since_last_perturbation++;

        if(debug){
            clog << endl << "\t***********************iter: " << iter << ", nn_iter: " << nn_iter << ", best: "
                 << best_result << " current: " << current_result << endl;
        }else if(!Global::disable_all_logs){ // #TEST
            cerr << "\rNEG: " << iter << ", perturbations_done: " << perturbations_done << ", nn_iter: " << nn_iter
                 << ", best: " << best_result << " current: " << current_result << ", #clusters: "
                 << nonempty_clusters_cnt << ", time_since_start: " << Global::secondsFromStart() << flush;
        }

        bool improved = false;

        {
            queue.clear();
            localShuffle(perm);
            if(debug){ DEBUG(perm);DEBUG(*clg); }
            for (int d : perm){ in_queue[d] = true; }

            const double PERM_FRACTION = perm_fraction;
            queue.assign( perm.begin(), perm.begin() + ( PERM_FRACTION * N) );
        }

        int step = move_frequency;

        while( !queue.empty() ){
            if(Global::checkTle()){
                if( current_result <= best_result ){
                    best_result = current_result;
                    best_partition = PaceUtils::properlyRemapPartition(inCl);
                }
                return;
            }

            iter_nodes_to_check.clear();
            while( !queue.empty() && iter_nodes_to_check.size() < move_frequency ){
                int d = queue.front();
                queue.pop_front();
                in_queue[d] = false;
                iter_nodes_to_check.push_back(d);
            }

            // this code only to use with queue propagation
            int a = 0;
            int b = iter_nodes_to_check.size() - 1;

            //*****************  NODE MOVE
            auto [v,to,val] = getBestNodeMoveForRange(iter_nodes_to_check, a,b);

            if(debug){
                debug_all();
                DEBUG3(v,to,val);
            }

            if( v != -1 ){ // is use_node_swaps is false, then v == -1
                if( val < 0 ){
                    improved = true;
                }
                if(val <= perturb_swp_thr){
                    if(debug){ clog << "Moving " << v << " to " << to << ", swpval: " << val << endl; }
                    moveNodeTo( v, to );
                    current_result += val;
                }
                if(val < 0) addClusterNodesToQueue(to);
            }

            //*****************  EDGE MOVE
            bool use_edges_swaps_in_iteration =
                    use_edge_swaps &&
                    //                    ((nn_iter == last_nn_iter-2) ||
                    ((nn_iter == min(max(1,last_nn_iter-2), 3)) || // #TEST - early edge swaps
                     ((iters_since_last_perturbation % EDGE_SWAPS_FREQUENCY) == EDGE_SWAPS_FREQUENCY-1));
            // #TEST - using edges swaps every 20 node swap iterations
            /**
             * using edge swaps only if we are in nonnegative iterations, because it is rather slow
             * compared to node moves
             */
            if(use_edges_swaps_in_iteration){
                auto [e,to,val] = getBestEdgeMoveForRange(iter_nodes_to_check, a,b);

                if( val < 0 ){
                    improved = true;
                    addClusterNodesToQueue(inCl[e.first]);
                    addClusterNodesToQueue(inCl[e.second]);
                }
                if( val <= perturb_swp_thr ){

                    moveNodeTo(e.first, to.first);
                    moveNodeTo(e.second, to.second);

//                    if( to.first != to.second && val < 0 ) clog << "Applying edge repulsion with swpval " << val << endl;

                    current_result += val;
                }
                if(val < 0){
                    addClusterNodesToQueue(to.first);
                    addClusterNodesToQueue(to.second);
                }
            }

            //*****************  TRIANGLE MOVE
            bool use_triangle_swaps_in_iteration =
                    (use_triangle_swaps > 0) &&
                    //                    ((nn_iter == last_nn_iter) ||
                    ((nn_iter == min(max(3,last_nn_iter-1), 5)) || // #TEST
                     ((iters_since_last_perturbation % TRIANGLE_SWAPS_FREQUENCY) == TRIANGLE_SWAPS_FREQUENCY-1));
            // #TEST - using triangle swaps
            if(use_triangle_swaps_in_iteration){
                auto best = getBestTriangleDiffClToMove(iter_nodes_to_check, a,b);
                if( best.swpVal() < 0 ){
                    improved = true;
                }
                if( best.swpVal() <= perturb_swp_thr ){
//                    if(!Global::disable_all_logs && !is_empty_cluster[best.getMoveNodesTo()[0]] && best.swpVal() < 0){
//                        clog << "Found triangle swap with swpval: " << best.swpVal() << " to "
//                        << (is_empty_cluster[best.getMoveNodesTo()[0]] ? "" : "NON") << "empty cluster" << endl;
//                    }

                    current_result += best.swpVal();
                    VI nodes_before_swap = best.nodes;
                    for( auto [v,c] : best.getNodesToSwap() ) moveNodeTo(v,c);
                }
                if(best.swpVal() < 0){
                    for( auto c : best.getNodes() ) addClusterNodesToQueue(inCl[c]); // #TEST - moved from above
                    addClusterNodesToQueue( best.move_node_to[0] );
                }

//                if(!Global::CONTEST_MODE) assert( compareCurrentResultWithBruteResult() ); // FIXME:remove
            }

//            //*****************  NODE INTERCHANGING
            if(use_node_interchanging_in_inner_loop){ // this should happen only in NEG_map
                bool use_node_interchanging_in_iteration =
                        use_node_interchanging &&
                        //                    ((nn_iter == last_nn_iter)||
                        ((nn_iter == min(max(2, last_nn_iter - 2), 5)) || // #TEST
                         ((iters_since_last_perturbation % NODE_INTERCHANING_FREQUENCY) == NODE_INTERCHANING_FREQUENCY - 1));
                if (use_node_interchanging_in_iteration) {

                    auto best = getBestInterchangeNodePairForInterval(iter_nodes_to_check, a, b);
                    val = best.swpVal();
                    VPII to_swap = best.getNodesToSwap();

                    if (val < 0) {
                        improved = true;
                        for (auto[x, t] : to_swap) addClusterNodesToQueue(t);
                        if (!Global::disable_all_logs)
                            clog << endl << endl << "Improved solution by val: " << val << " using node interchange"
                                 << endl << endl;
                    }
                    if (val <= perturb_swp_thr) {
                        if (debug) {
                            clog << endl << endl << "interchanging nodes " << to_swap[0].first << " and "
                                 << to_swap[1].first << ", to_swap: " << to_swap << endl << endl;
                        }
                        for (auto[x, t] : best.getNodesToSwap()) moveNodeTo(x, t);
                        current_result += val;
                    }
                }
            }
        }
        // *********************************************************  end of queue

        if(Global::checkTle()) return;

        if(!use_node_interchanging_in_inner_loop){ // using node interchanging in outer loop
            //*****************  NODE INTERCHANGING
            bool use_node_interchanging_in_iteration =
                    use_node_interchanging &&
//                                        ((nn_iter == last_nn_iter-1)||
                    ((nn_iter == min(max(2,last_nn_iter-2), 5)) || // #TEST
                     ((iters_since_last_perturbation % NODE_INTERCHANING_FREQUENCY) == NODE_INTERCHANING_FREQUENCY-1)  );
            if(use_node_interchanging_in_iteration){
                createClusterNodes();
                for( int i=0; i<N; i++ ) createEdgesToCluster(i);

                iter_nodes_to_check = VI(ALL(perm));
                auto best = getBestInterchangeNodePairForInterval(iter_nodes_to_check, 0, iter_nodes_to_check.size()-1);
                auto val = best.swpVal();
                VPII to_swap = best.getNodesToSwap();

                if( val < 0 ){
                    improved = true;
                    for( auto [x,t] : to_swap ) addClusterNodesToQueue(t);
                    if(!Global::disable_all_logs)
                        clog << endl << endl << "Improved solution by val: " << val << " using node interchange" << endl << endl;
                }
                if(val <= perturb_swp_thr){
                    if(debug){
                        clog << endl << endl << "interchanging nodes " << to_swap[0].first << " and "
                             << to_swap[1].first << ", to_swap: " << to_swap << endl << endl;
                    }
                    for( auto [x,t] : best.getNodesToSwap() ) moveNodeTo( x,t );
                    current_result += val;
                }
            }
        }

        if(Global::checkTle()) return;

        bool use_chain2_swaps_in_iteration =
                use_chain2_swaps &&
//                        (( nn_iter == last_nn_iter-1)  || //#TEST - use only in almost last iteration -as it is rather time-consuming
                    ((nn_iter == min(max(3,last_nn_iter-2), 4)) || // #TEST
                  ((iters_since_last_perturbation % CHAIN2_FREQUENCY) == CHAIN2_FREQUENCY-1)  );
        if( use_chain2_swaps_in_iteration ){
//            clog << "\tUsing chain2 swaps" << endl;
            int prev_current = current_result;
            makeChain2Swaps();
            int diff = current_result - prev_current;
            if(diff < 0){
//                clog << endl << endl << "chain2 improved by " << diff << endl << endl;
                improved = true;
            }
        }

        bool use_join_clusters_in_iteration =
                use_join_clusters &&
                        (nn_iter == last_nn_iter); // || //#TEST - do not use here as it is rather time-consuming
//                 (((iters_since_last_perturbation % JOIN_CLUSTERS_FREQUENCY) == JOIN_CLUSTERS_FREQUENCY-1));
        if( use_join_clusters_in_iteration ){
//            clog << "\tUsing cluster joins" << endl;
            auto to_join = getBestClustersToJoin();
            if( !to_join.empty() ){
//                clog << endl << "Joining clusters!, current_value: " << current_result << flush;
                improved = true;

                createClusterNodes();
                VB affected( 2*N,false );
                sort( ALL(to_join), [](auto & a, auto & b){
                    return get<2>(a) < get<2>(b); // sorting by least swap value - the least, the better
                } );

                for( auto [c1,c2,val] : to_join ){
                    if( affected[c1] || affected[c2] ) continue;
                    affected[c1] = affected[c2] = true;
                    current_result += val;
                    if( c1 < c2 ) swap(c1,c2); // we want to move to cluster with smaller id
                    VI to_move(ALL(cluster_nodes[c1]));
                    for( int d : to_move) moveNodeTo(d,c2);
                }
                cluster_nodes.clear();
//                clog << "  |  After joining clusters current_value: " << current_result << endl;
            }
        }

        if(Global::checkTle()) return;

        if(!Global::CONTEST_MODE) assert( compareCurrentResultWithBruteResult() ); // FIXME:remove

        if(!improved){
            nn_iter++;
        }
        else nn_iter = 0;

        perturb_swp_thr = 0;


        if(nn_iter == max_nonnegative_iters ) {
            if(!allow_perturbations) break;

            perturb(perturb_mode);
            if( do_not_perturb_if_improved && best_result < initial_result ) return;

            perturbations_done++;
            nn_iter = 0;
            iters_since_last_perturbation = 0;

            shuffleClg();

            if((perturbations_done % max(1,(max_perturbations/3))) == 0){
                prefer_cluster_mode++;
                prefer_cluster_mode %= 3;
            }
        }
    }

    if( current_result <= best_result ){
        best_result = current_result;
        best_partition = PaceUtils::properlyRemapPartition(inCl);
    }

    if(!Global::disable_all_logs){
        auto mapped_part = PaceUtils::mapClgPartitionToOriginalPartition(*clg, best_partition);
        tuple<int, int, int> mods = PaceUtils::getEdgeModificationStatistics(*clg->origV,mapped_part);
        clog << ", best_res_mods: " << mods << flush;
    }

}

tuple<int,int,int> NEG::getBestNodeMoveForRange(VI &perm, int a, int b) {
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

        int w_0 = findEdgesToCluster(d,cl_d);

        int tot_cld_possible_edges = (cluster_weights[cl_d] - clg->node_weights[d]) * clg->node_weights[d];

        VPII etocld = getEdgesToCluster(d);
        for( auto & [c,w] : etocld ){
            if( c == cl_d ) continue;
            int swpval = swapValueForNode(d, c, w_0, tot_cld_possible_edges, w);

            if(swpval < best_swpval){
                best_swpval = swpval;
                best_v = d;
                best_to = c;

                best_node_move_results.clear(); // #TEST
            }

            if( swpval <= best_swpval ) best_node_move_results.emplace_back( d,c,swpval ); // #TEST
        }


        if( cluster_weights[cl_d] != nw_d ){
            // if there is more than one node in cluster contatinig d, then we try to move it to an empty cluster
            int c = *first_free_cluster.begin();
            int w = 0;
            resizeStructuresForEmptyCluster(c);

            int swpval = swapValueForNode(d, c, w_0, tot_cld_possible_edges, w);

            if(swpval < best_swpval){
                best_swpval = swpval;
                best_v = d;
                best_to = c;

                best_node_move_results.clear(); // #TEST
            }

            if( swpval <= best_swpval ) best_node_move_results.emplace_back( d,c,swpval ); // #TEST

            if(return_on_first_negative_swap && best_swpval < 0) break;
        }

        if(return_on_first_negative_swap && best_swpval < 0) break;
    }

    const bool USE_RANDOM_BEST_MOVE = true; // #TEST - getting random one of all best possible moves // originally false
    if(USE_RANDOM_BEST_MOVE && !best_node_move_results.empty()){
        int ind = (perm[0] + perm.back()) % best_node_move_results.size();

        return best_node_move_results[ind];
    }

    return { best_v, best_to, best_swpval };
}

int NEG::swapValueForNode(int v, int trg_cl, int edges_clv, int tot_clv_possible_edges, int edges_trg) {
    int before = tot_clv_possible_edges - edges_clv + edges_trg;
    int after = cluster_weights[trg_cl] * clg->node_weights[v] - edges_trg + edges_clv;
    return after-before;
}

tuple<PII, PII, int> NEG::getBestEdgeMoveForRange(VI &perm, int a, int b) {
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


    for( int i=a; i<=b; i++ ) {
        int u = perm[i];
        int cl_u = inCl[u];
        int nw_u = clg->node_weights[u];

        if(debug){
            ENDL(2);
            DEBUG(u); DEBUG(cl_u); DEBUG(nw_u);
        }

        ConvexHullTrickDynamic hull(false); // minimizing
        const int empty_cluster = *first_free_cluster.begin();

        createEdgesToCluster(u);
        VPII etoclu = getEdgesToCluster(u);

        if(USE_HULL_TRICK_IN_EDGE_SWAPS){ // insert to hull
            hull.insert_line(0,0, empty_cluster);

            for (auto & [c, w2] : etoclu ) {
                if( c == cl_u ) continue;
                helper[c] += w2; // helper is just edges_to_cluster of both nodes d and d2

                int cwc = cluster_weights[c];
                if( min_w2_to_cluster_of_given_weight[cwc] > -(w2 << 1) ) { // #TEST - originally was >=
                    hull.insert_line(cwc, -(w2 << 1), c);
                    min_w2_to_cluster_of_given_weight[cwc] = -(w2 << 1);
                }

                if(debug) clog << "Adding line for cluster " << c << ": " << cluster_weights[c] << " * x - " << (w2<<1) << endl;
            }

            for (auto & [c, w2] : etoclu) min_w2_to_cluster_of_given_weight[cluster_weights[c]] = 0; // clearing
        }else{
            for (auto & [c, w2] : etoclu) {
                if (c == cl_u) continue;
                helper[c] += w2; // helper is just edges_to_cluster of both nodes d and d2
            }
        }

        int etc_u_clu = findEdgesToCluster(u,cl_u);

        for( auto & [v,w0] : clg->V[u] ){
            createEdgesToCluster(v);
            auto etoclv = getEdgesToCluster(v);

            // this 2* factor is just to increase variability whilst preserving good complexity
            // selecting FACTOR 1.5 or 2 should be ok
            const double FACTOR = 1.2; // original 1.5
            const double ADD = 5;  // original 0
            if( etoclv.size() > FACTOR*etoclu.size() + ADD ) continue;

            int cl_v = inCl[v];
            int nw_v = clg->node_weights[v];

            int residue = 0, s = nw_u + nw_v;

            int etc_v_clv = findEdgesToCluster(v, cl_v);

            if( cl_v == cl_u ){
                int w1 = etc_u_clu + etc_v_clv - (w0<<1);
                residue = (w1<<1) - s * ( cluster_weights[cl_u] - s );
            }else{
                int w1 = etc_u_clu + etc_v_clv;
                residue = nw_u * ( nw_u - cluster_weights[cl_u] ) + nw_v * ( nw_v - cluster_weights[cl_v] )
                          + nw_u * nw_v - 2*( w0-w1 );
            }


            int best_c, best_c_val;

            if(USE_HULL_TRICK_IN_EDGE_SWAPS){
                best_c = hull.get_line_id(s);
                best_c_val = 1e9;
            }else{
                cerr << "Hull trick is not used, clsuters to check is not set, terminating" << endl;
                assert( false );
                best_c = empty_cluster;
                best_c_val = 1e9;
            }

            resizeStructuresForEmptyCluster(best_c);

            if(debug){
                ENDL(1);
                DEBUG(v); DEBUG(w0); DEBUG(cl_v); DEBUG(nw_v); DEBUG(residue);
            }

            if(best_c != cl_u && best_c != cl_v){
                best_c_val = cluster_weights[best_c] * s - ( ( findEdgesToCluster(u,best_c) )<<1 )  + residue;
            }
            else if(debug){
                clog << "Trying to move edge to cluster containing one of its ends" << endl;
                DEBUG(best_c); DEBUG(cl_u); DEBUG(cl_v);
            }


            if( cl_u == cl_v && nw_u + nw_v == cluster_weights[cl_u] && cluster_weights[best_c] == 0 ){
                // if we would like to move a whole two-noded cluster to an empty cluster
                best_c_val = 1e9;
            }

            if(debug){
                DEBUG(best_c); DEBUG(best_c_val);
            }

            for( auto [c,w2] : etoclv ){ // original #CAUTION this will not work is HULL TRICK IS NOT USED! But should iterate faster
                if( c == cl_u || c == cl_v ) continue;
                w2 += helper[c];
                int c_val = cluster_weights[c] * s - (w2<<1) + residue;
                if(debug){
                    DEBUG(c); DEBUG(w2); DEBUG(c_val);
                }

                if( c_val < best_c_val ){
                    if(debug){
                        clog << "after updating neighbors of " << v << ", c_val: " << c_val << " < " << best_c_val
                             << " = best_c_val" << endl;
                    }
                    best_c_val = c_val;
                    best_c = c;
                }
            }

            if(best_c_val < best_swpval){
                best_swpval = best_c_val;
                best_e = {u,v};
                best_to = {best_c,best_c};
            }

            if(best_swpval < 0 && return_on_first_negative_swap) break;
        }

        for (auto & [c, w2] : etoclu) helper[c] = 0; // clearing helper;
        if(best_swpval < 0 && return_on_first_negative_swap) break;
    }

    return {best_e, best_to, best_swpval};
}


void NEG::perturb( int perturbation_mode ) {
    if( current_result <= best_result ){
        best_result = current_result;
        best_partition = PaceUtils::properlyRemapPartition(inCl);
    }

    if( do_not_perturb_if_improved && best_result < initial_result ) return;

    if(Global::checkTle()) return;

    if(perturbation_mode == 0){
        splitClustersIntoTwo();
    }
    else if(perturbation_mode == 1){
        joinClustersInPairs();
    }
    else if( perturbation_mode == 2 ){
        applyExpOrdREP();
    }

    {
        VI part = PaceUtils::properlyRemapPartition( inCl );
        VI partition(clg->origV->size());
        for( int i=0; i<part.size(); i++ ) for( int d : clg->clusterNodes[i] ) partition[d] = part[i];
        current_result = PaceUtils::evaluateSolution(*clg->origV, partition);
    }

    nonempty_clusters_cnt = countNonemptyClusters();
}

void NEG::splitClustersIntoTwo() {
    const bool debug = false;

    VVI clusters = StandardUtils::partitionToLayers(inCl);
    unordered_set<int> zb(ALL(inCl));
    assert(zb.size() <= clusters.size()+1);

    if(debug){
        int sz = 0;
        for( VI & v : clusters ) if( !v.empty() ) sz++;
        clog << "Splitting clusters into two, total clusters:" << sz << endl;
        DEBUG(first_free_cluster.size());
    }

    int MIN_CLUSTER_SIZE_TO_SPLIT = 3; // #TEST - originally 2, value 4 seems to be ok

//    UniformDoubleGenerator rnd(0.25, 0.45);

    for( VI & v : clusters ){
        if(v.size() < MIN_CLUSTER_SIZE_TO_SPLIT) continue;

//        int S = (v.size())>>1;
//        int S = (v.size()) / 5;
//        int S = (v.size())>>2; // original value, seems to work very good
//        int S = v.size() * rnd.rand();

        int S; // #TEST - perturbations split size based on iteration
        int X = 4; // #TEST oriignal 3
        if( (perturbations_done % X ) < X-1 ) S =  (v.size())>>2;
        else S = v.size() * 0.4;

        S = max(S,1);

        VI subs = CombinatoricUtils::getRandomSubset(v.size()-1, S );
        int c = *first_free_cluster.begin();

        for( int ind : subs ){
            int d = v[ind];
            moveNodeTo(d,c);
        }
    }

    if(debug){
        VVI cl = StandardUtils::partitionToLayers(inCl);
        int sz = 0;
        for( VI & v : cl ) if( !v.empty() ) sz++;

        clog << "After splitting, clusters size:" << sz << endl;
        DEBUG(first_free_cluster.size());
    }
}

void NEG::joinClustersInPairs() {

    VVI clusters = StandardUtils::partitionToLayers(inCl);

    vector< pair<double,PII> > edges;

    VB was(helper.size(),false);

    for( int i=0; i<clusters.size(); i++ ){
        if( clusters[i].empty() ) continue;

        VI cl_neigh;
        for( int d : clusters[i] ){
            createEdgesToCluster(d);
            VPII etocld = getEdgesToCluster(d);
            for( auto & [c,w] : etocld ){
                if(c != inCl[d]){
                    helper[c] += w;
                    if(!was[c]){
                        cl_neigh.push_back(c);
                        was[c] = true;
                    }
                }
            }
        }

        for( int c : cl_neigh ) {
            double density = (double)helper[c] / (cluster_weights[ inCl[clusters[i][0]] ] * cluster_weights[c]);
            edges.push_back({density, {i, c}});
        }

        for( auto c : cl_neigh ){
            helper[c] = 0; // clearing
            was[c] = false;
        }
        cl_neigh.clear();
    }

    UniformIntGenerator rnd(0,1);
    StandardUtils::shuffle( edges, rnd.getRNG() );

    sort(ALL(edges), [](auto& a, auto& b){
        return a.first > b.first; // sorting by density
    });

    for( auto & [dens, e] : edges ){
        int a = e.first, b = e.second;
        if( was[a] || was[b] ) continue;

        was[a] = was[b] = true;

        for( int d : clusters[a] ){
            moveNodeTo(d,b);
        }
    }

}


//void NEG::moveEdgeTo(int v, int w, int to) {
//    const bool debug = false;
//    if(debug){
//        clog << "moving edge (" << v << "," << w << ") to " << to << endl;
//    }
//
//    moveNodeTo(v,to);
//    moveNodeTo(w,to);
//}

void NEG::setConfigurations(Config &cnf) {
    max_perturbations = cnf.neg_max_perturb;
    max_nonnegative_iters = cnf.neg_max_nonneg_iters;
    use_edge_swaps = cnf.neg_use_edge_swaps;
    use_triangle_swaps = cnf.neg_use_triangle_swaps;
    use_queue_propagation = cnf.neg_use_queue_propagation;
    use_node_interchanging = cnf.neg_use_node_interchange;
    use_join_clusters = cnf.neg_use_join_clusters;
    use_chain2_swaps = cnf.neg_use_chain2_swaps;
    use_two_node_swaps = cnf.neg_use_two_node_swaps;
    perm_fraction = cnf.neg_perm_fraction;

    move_frequency = cnf.neg_move_frequency; // #TEST - this was not here earlier

    edge_swaps_frequency = cnf.neg_edge_swaps_frequency;
    triangle_swaps_frequency = cnf.neg_triangle_swaps_frequency;
    do_not_perturb_if_improved = cnf.neg_do_not_perturb_if_improved;
    chain2_swaps_frequency = cnf.neg_chain2_swaps_frequency;
    node_interchanging_frequency = cnf.neg_node_interchanging_frequency;
    use_edge_repulsion = cnf.neg_use_edge_repulsion;

    max_best_cl_size_triangle_swaps = cnf.neg_max_best_cl_size_triangle_swaps;
    use_triangle_swaps_to_other_clusters = cnf.neg_use_triangle_swaps_to_other_clusters;
}

void NEG::addClusterNodesToQueue(int cl_id) {
    if(!use_queue_propagation) return; // do not use queue propagation

    for( int v : cluster_nodes[cl_id] ){
        if( !in_queue[v] ){
            in_queue[v] = true;
            queue.push_front(v); // pushing to front of the queue
        }
    }
}

int NEG::countNonemptyClusters() {
    return unordered_set<int>( ALL(inCl) ).size();
}

void NEG::applyExpOrdREP() {

    State st( *clg, SINGLE_NODES );
    VI remapped = PaceUtils::properlyRemapPartition(inCl);
    st.applyPartition( remapped );

    VI remapper;
    {
        remapper = VI( st.clusters.size()+1, -1 );
        for( int i=0; i<N; i++ ){
            remapper[ st.inCl[i] ] = inCl[i];
        }
    }

    int C = st.clusters.size();
    VI perm = CombinatoricUtils::getRandomPermutation(C-1); // we take C-1 because the last cluster is empty

    VB affected_clusters(C,false);

    int clusters_REPed = 0;

    for( int c : perm ){

        ComponentExpansionRepulsion eo_rep( st );
        eo_rep.keep_only_nonpositive_candidates = false; // we want the middle or the last candidate

        auto [ord, candidates] = eo_rep.createSwapCandidates( st.clusters[c] );

        if( candidates.empty() ){
            delete ord;
            continue;
        }

        auto last_cand = candidates[ candidates.size()-1 ];

        { // applying candidate
            clusters_REPed++;
            VPII to_swap = last_cand.getNodesToSwap();

            unordered_map<int,int> id_mapper;

            for( PII p : to_swap ){
                int v = p.first;
                int to = p.second;

                if( to >= st.getIdOfEmptyCluster() ){
                    auto it = id_mapper.find(to);
                    if( it == id_mapper.end() ){
                        int old_to = to;
                        to = *first_free_cluster.begin();
                        id_mapper[old_to] = to;
                    }
                    else to = it->second;
                }else{
                    assert( to < remapper.size() );
                    to = remapper[to];
                }

                moveNodeTo( v, to);

                while(affected_clusters.size() <= to) affected_clusters.push_back(false);
                affected_clusters[to] = true;
            }
        }

        delete ord;
    }

}

SwapCandidateAdapter NEG::getBestInterchangeNodePair(int v) {
    const bool debug = false;

    if(debug){
        ENDL(3);
        clog << "Checking node " << v << " for interchanging" << endl;
//        DEBUG(cluster_nodes);
    }

    SwapCandidateAdapter best;
    best.swap_value = 1e9;

    int cl_v = inCl[v];
    int clw_v = cluster_weights[cl_v];
//    int e_v_tocl_v = edges_to_cluster[v][cl_v];
    int e_v_tocl_v = findEdgesToCluster(v,cl_v);
    int nw_v = clg->node_weights[v];

    for( auto & [u,w] : clg->V[v] ) { // checking all edges
        int cl_u = inCl[u];
        if(cl_u == cl_v) continue;
        int clw_u = cluster_weights[cl_u];

        if(debug) clog << "Checking edge " << PII(v,u) << " with weight " << w << endl;

//        int val = getInterchangeValue(v, u, w);
        int val = getInterchangeValue(v, u, w, cl_v, nw_v, clw_v, cl_u, clw_u,
//                                      e_v_tocl_v, edges_to_cluster[v][cl_u] ); // #TEST - simplified version
                                      e_v_tocl_v, findEdgesToCluster(v,cl_u) ); // #TEST - simplified version

        if(debug) DEBUG(val);

        if (val < best.swap_value) {
            best.swap_value = val;
            best.nodes = {v, u};
            best.move_node_to = {cl_u, cl_v};
        }

        checked_for_v[u] = true;
    }

    auto etoclv = getEdgesToCluster(v);
    for( auto & [c,w] : etoclv ) {
        if( c == cl_v ) continue;
        int clw_c = cluster_weights[c];


//        **************
        int cl_u = c;
        int clw_u = cluster_weights[cl_u];
        int ev_to_clu = w;
//        **************

        if( w < PERC * clw_c ) continue;

        if(debug) clog << "Checking cluster " << c << endl;

        auto clnodesc = getClusterNodes(c);
        for( int d : clnodesc ){
            if( checked_for_v[d] ) continue; // edge (v,d) was checked earlier

            if(debug) clog << "Checking non-edge " << PII(v,d) << endl;

//            int val = getInterchangeValue(v, d, 0); // 0 because there is no edge between v and d
            int val = getInterchangeValue(v, d, 0, cl_v, nw_v, clw_v, cl_u, clw_u,
                                          e_v_tocl_v, ev_to_clu ); // #TEST - simplified version

            if(debug) DEBUG(val);

            if (val < best.swap_value) {
                best.swap_value = val;
                best.nodes = {v, d};
                best.move_node_to = {c, cl_v};
            }
        }

    }

    for( auto & [u,w] : clg->V[v] ) checked_for_v[u] = false;

    if(debug){
        clog << "best interchaning pair: " << best.getNodesToSwap() << ", swpval: " << best.swpVal() << endl;
    }

    return best;
}

int NEG::getInterchangeValue(int &v, int &u, int w_vu, int &cl_v, int &nw_v, int &clw_v, int &cl_u, int &clw_u,
                                        int &ev_to_clv, int ev_to_clu) {
    int nw_u = clg->node_weights[u];
    int eu_to_clu = findEdgesToCluster(u,cl_u);
    int eu_to_clv = findEdgesToCluster(u,cl_v);

    int before = 0;
    before += (clw_v - nw_v) * nw_v - ev_to_clv;
    before += (clw_u - nw_u) * nw_u - eu_to_clu;
    before += ev_to_clu - w_vu;
    before += eu_to_clv - w_vu;
    before += w_vu;

    int after = 0;
    after += ( clw_v - nw_v ) * nw_u - ( eu_to_clv - w_vu );
    after += ( clw_u - nw_u ) * nw_v - ( ev_to_clu - w_vu );
    after += ev_to_clv;
    after += eu_to_clu;
    after += w_vu;

    return after - before;


    { // simplified version
        int swpval = (ev_to_clv<<1) + (eu_to_clu<<1) - ( (eu_to_clv + ev_to_clu) << 1) + (w_vu<<2);
        swpval += nw_v * ( clw_u - nw_u - clw_v + nw_v ) + nw_u * ( clw_v - nw_v - clw_u + nw_u );

        assert(swpval == after - before);

        return swpval;
    }

}



int NEG::getInterchangeValue(int v, int u, int w_vu) {
    const bool debug = false;

    int cl_v = inCl[v];
    int nw_v = clg->node_weights[v];
    int clw_v = cluster_weights[cl_v];

    int cl_u = inCl[u];
    int nw_u = clg->node_weights[u];
    int clw_u = cluster_weights[cl_u];

//    int ev_to_clv = edges_to_cluster[v][cl_v];
    int ev_to_clv = findEdgesToCluster(v,cl_v);
//    int eu_to_clu = edges_to_cluster[u][cl_u];
    int eu_to_clu = findEdgesToCluster(u,cl_u);
//    int ev_to_clu = edges_to_cluster[v][cl_u];
    int ev_to_clu = findEdgesToCluster(v,cl_u);
//    int eu_to_clv = edges_to_cluster[u][cl_v];
    int eu_to_clv = findEdgesToCluster(u,cl_v);

    int before = 0;
    before += (clw_v - nw_v) * nw_v - ev_to_clv;
    before += (clw_u - nw_u) * nw_u - eu_to_clu;
    before += ev_to_clu - w_vu;
    before += eu_to_clv - w_vu;
    before += w_vu;

    int after = 0;
    after += ( clw_v - nw_v ) * nw_u - ( eu_to_clv - w_vu );
    after += ( clw_u - nw_u ) * nw_v - ( ev_to_clu - w_vu );
    after += ev_to_clv;
    after += eu_to_clu;
    after += w_vu;

    if(debug){
        ENDL(2);
        DEBUG3(v,cl_v,nw_v);
        DEBUG3(clw_v, ev_to_clv, ev_to_clu);

        DEBUG3(u,cl_u,nw_u);
        DEBUG3(clw_u, eu_to_clu, eu_to_clv);

        DEBUG2(before,after);
        DEBUG(after - before);
    }

    // FIXME:optimize reduce same values to decrease constant
    return after - before;
}

SwapCandidateAdapter NEG::getBestInterchangeNodePairForInterval(VI & nodes_to_check, int a, int b) {
    SwapCandidateAdapter best;
    best.swap_value = 1e9;

    for( int i=a; i<=b; i++ ){
        int v = nodes_to_check[i];
        auto temp = getBestInterchangeNodePair(v);
        if( temp.swpVal() < best.swpVal() ) best = temp;
    }

    return best;
}

vector<tuple<int,int,int>> NEG::getBestClustersToJoin(bool allow_zero_swpval) {
    const bool debug = false;

    vector<tuple<int,int,int>> res;

    int C = maxClusterId();
    VI perm = CombinatoricUtils::getRandomPermutation(C);

    createClusterNodes();
    for(int i=0; i<N; i++) createEdgesToCluster(i);

    VI cl_neigh;
    for( int i : perm ){
        auto cluster_nodes_i = getClusterNodes(i);
        if( cluster_nodes_i.empty() ) continue;

        cl_neigh.clear();
        for( int d : cluster_nodes_i ){
            VPII e_to_cl_d = getEdgesToCluster(d);
            for( auto & [c,w] : e_to_cl_d ){
                if( c != i ){
                    helper[c] += w;
                    if( !helper_was[c] ){
                        helper_was[c] = true;
                        cl_neigh.push_back(c);
                    }
                }
            }
        }

        for( int c : cl_neigh ){

            int before = helper[c];
            int after = cluster_weights[i] * cluster_weights[c] - helper[c];
            int swpval = after - before;

            if( after < before ){
                res.emplace_back(i,c,swpval);
            }else if( after == before && allow_zero_swpval ){
                res.emplace_back(i,c,swpval);
            }
        }

        for( int c : cl_neigh ) helper_was[c] = helper[c] = 0; // clearing
    }

    cluster_nodes.clear();

    return res;
}

bool NEG::compareCurrentResultWithBruteResult() {
    VI part = PaceUtils::properlyRemapPartition(inCl);
    part = PaceUtils::mapClgPartitionToOriginalPartition(*clg, part);
    int brute_res = PaceUtils::evaluateSolution( *clg->origV, part );

    if(current_result != brute_res ){
        DEBUG(current_result);
        DEBUG(brute_res);
    }

    return brute_res == current_result;
}


SwapCandidateAdapter NEG::getBestTriangleDiffClToMove(VI &nodes_to_check, int a, int b) {
    SwapCandidateAdapter best;
    best.swap_value = 1e9;

    for( int i=a; i<=b; i++ ){
        int v = nodes_to_check[i];
        auto temp = getBestTriangleAll(v); // #TEST - trying all triangle moves, not only to empty cluster
        if(temp.swpVal() < best.swpVal()) best = temp;
    }

    return best;
}

int NEG::countClusterNumerationGaps() {
    int empty_cnt = 0;
    int res = 0;

    createClusterNodes();

    for( int i=0; i<cluster_nodes.size(); i++ ){
        if( !cluster_nodes[i].empty() ) res = empty_cnt;
        else empty_cnt++;
    }

    cluster_nodes.clear();

    return res;
}

bool NEG::makeChain2Swaps() {
    const bool debug = false;

    createClusterNodes();
    for( int i=0; i<N; i++ ) createEdgesToCluster(i);

    if(debug){
        ENDL(10);
        DEBUG(cluster_nodes);
        DEBUG(cluster_weights);
        DEBUG(first_free_cluster);
    }

    bool improved = false;
    const int EMPTY_CLUSTER = *first_free_cluster.begin();
//    while( cluster_weights.size() <= EMPTY_CLUSTER ){ cluster_weights.push_back(0);cluster_nodes.push_back( {} ); }
    resizeStructuresForEmptyCluster(EMPTY_CLUSTER);

    /**
     * possible_swaps[v] is a queue containing pairs (c,val) denoting that moving node v to cluster c has given
     * swap value
     */
    vector< deque<PII> > possible_swaps(N);

    { // creating possible swaps
        VPII temp; //enabling probably a bit quicker sorting than sorting than deque

        for( int i=0; i < cluster_nodes.size(); i++ ){
            if( cluster_nodes[i].empty() ) continue;
            int cl_u = i; // just to improve readability and debugging

            if( debug ) clog << "********* Considering cluster " << cl_u << endl;

            for( int u : cluster_nodes[i] ){
//                int eu_to_clu = edges_to_cluster[u][cl_u];
                int eu_to_clu = findEdgesToCluster(u,cl_u);

                int nw_u = clg->node_weights[u];
                int tot_eu = (cluster_weights[cl_u] - nw_u) * nw_u;

                if(debug){
                    DEBUG3(u,cl_u,nw_u); DEBUG2(eu_to_clu, tot_eu);
                }

                temp.clear();
                if( nw_u != cluster_weights[i] ){
                    // we do not add empty cluster if u is alone in a cluster
                    temp.emplace_back( EMPTY_CLUSTER, swapValueForNode( u, EMPTY_CLUSTER, eu_to_clu, tot_eu, 0 ) );
                }
                auto etoclu = getEdgesToCluster(u);
                for( auto & [c,w] : etoclu ){
                    if( c == cl_u ) continue;
                    // swap value of moving node u to cluster c
                    int swpval = swapValueForNode( u, c, eu_to_clu, tot_eu, w );
                    temp.emplace_back( c,swpval );
                }

                sort(ALL(temp), []( PII a, PII b ){
                    return a.second < b.second;
                });


                possible_swaps[u].assign( ALL(temp) );

                if(debug){
                    DEBUG(possible_swaps[u]);
                    ENDL(1);
                }
            }
        }
    }

    if(Global::checkTle()) return false;
    // possible swaps should be created now

    VB affected_clusters( cluster_weights.size(),false );
    VI weight_to_node(N,0); // helper, edges between u and given node
    VI swpval_to_cluster(cluster_weights.size(), 0); // helper, swap value of u to given cluster

    const double FACT = 0.4;
    const PII invalid(-1,-1);
    vector<SwapCandidateAdapter> to_move; // vector of all moves that are to be applied in the end

    auto makeMovesForSwpThreshold = [&](int threshold){
        VI perm = CombinatoricUtils::getRandomPermutation(N);
        for( int u : perm ){
            int cl_u = inCl[u];
            if(affected_clusters[cl_u]) continue;

            int nw_u = clg->node_weights[u];
            bool found = false;

            if(debug){
                clog << endl << "**************************" << endl;
                DEBUG3(u,cl_u,nw_u);
                DEBUG(possible_swaps[u]);
                DEBUG(affected_clusters);
            }

            for( auto & [v,w] : clg->V[u] ) weight_to_node[v] = w;
            for( auto & [c,val] : possible_swaps[u] ) swpval_to_cluster[c] = val;

            if(debug) DEBUG(swpval_to_cluster);

            auto etoclu = getEdgesToCluster(u);
            for( auto & [c1,w1] : etoclu ){
                if((c1 == cl_u) || affected_clusters[c1]) continue;
                if( w1 < cluster_weights[c1] * FACT ) continue;

                int swpval1 = swpval_to_cluster[c1];

                for( int v : cluster_nodes[c1] ){
//                    if(debug){ clog << "\t"; DEBUG(v); }

                    PII to_put_back = invalid;

                    int nw_v = clg->node_weights[v];
                    if( (nw_u == cluster_weights[cl_u]) || (nw_v == cluster_weights[c1]) ){
                        // if cl_u and c1 are single-noded clusters, then moving u to c1 and v to c2 would be equivalent
                        // to just moving v to c2, so just a node swap
                        continue;
                    }

                    int euv = weight_to_node[v]; // number of edges to cluster c1\v

                    if(debug){
                        ENDL(1);
                        DEBUG2(v,possible_swaps[v]);
                    }

                    while( !found && !possible_swaps[v].empty() ){
                        auto [c2,swpval2] = possible_swaps[v].front();
                        possible_swaps[v].pop_front();

                        if(debug){
                            DEBUG3(v,c2, swpval2);
                            clog << "Popped " << PII(c2,swpval2) << " from possible_swaps[" << v << "]" << endl;
                        }

                        if( affected_clusters[c2] ) continue;
                        if( c2 == cl_u ){ // we do not put the node back to cluster with u - that is interchanging
                            to_put_back = {c2,swpval2};
                            continue;
                        }

                        int swpval = swpval2 + swpval1 - nw_u * nw_v + ( euv << 1 );

                        if(debug){
                            ENDL(1);
                            clog << "Swap value for moving node " << u << " to cluster " << c1 << " and node " << v
                                 << " to cluster " << c2 << ": " << swpval << endl;
                            DEBUG(swpval1); DEBUG(swpval2); DEBUG(euv);
                        }

                        if( swpval <= threshold ){ // apply the move
                            SwapCandidateAdapter cnd;
//                            cnd.nodes = { u,v };
                            cnd.nodes.push_back(u); cnd.nodes.push_back(v); // #TEST
                            cnd.move_node_to = {c1, c2};
                            cnd.swap_value = swpval;
                            to_move.push_back( cnd );

                            affected_clusters[cl_u] = affected_clusters[c1] = affected_clusters[c2] = true;
                            affected_clusters[EMPTY_CLUSTER] = false;
                            if(debug) clog << "Marking affected clusters" << endl;
                            found = true;
                        }else{
                            // pushing the last entry to possible_swaps[v]
                            if(debug){
                                clog << "Pushing " << PII(c2,swpval2) << " back to possible_swaps[" << v << "]" << endl;
                            }
                            possible_swaps[v].push_front( {c2,swpval2} );
                        }

                        if(swpval < 0) improved = true;

                        // we break here - as we found the first swap candidate, it is the best swap candidate, since
                        // the swap value does not affect movinf v to c2, only moving u to c1. We will not find a better
                        // candidate than the one found here.
                        break;
                    }


                    // putting back the entry of moving v to cl_u if it occured
                    if( to_put_back != invalid ){
                        if(debug) clog << "Pushing to_put_back: " << to_put_back << " to possible_swaps[" << v << "]" << endl;
                        possible_swaps[v].push_front(to_put_back);
                    }

                    if( found ) break; // breaking, c1 can ba affected only once
                }

                if( found ) break; // breaking - cl_u can be affected only once
            }

            for( auto & [c,val] : possible_swaps[u] ) swpval_to_cluster[c] = 0; // clearing
            for( auto & [v,w] : clg->V[u] ){
                weight_to_node[v] = 0; // clearing
            }
        }
    };

    if(debug)clog << "Running with negative swpval only" << endl;
    makeMovesForSwpThreshold(-1); // now applying greedily the first found negative swap
    if(Global::checkTle()) return false;

    if(debug){ ENDL(2);clog << "Running with zero-swpval admitted" << endl; }
    makeMovesForSwpThreshold(perturb_swp_thr); // now applying possible zero-swpval (or perturb-swp_thr-val) movements (some perturbation)
    if(Global::checkTle()) return false;

    if(debug){
        if(!to_move.empty()) {
            clog << endl << "to_move.size(): " << to_move.size() << endl;
            for (auto &cnd : to_move) {
                DEBUG(cnd.nodes);
                DEBUG(cnd.move_node_to);
                DEBUG(cnd.swap_value);
                clog << "Affected clusters by cnd: (" << inCl[cnd.nodes[0]] << ", " << cnd.move_node_to[0] << ", "
                     << cnd.move_node_to[1] << ")" << endl;
                ENDL(1);
            }
        }
    }

    // applying the move
    for( auto & cnd : to_move ){
        assert(cnd.swpVal() <= 0);

        for( auto [v,c] : cnd.getNodesToSwap() ){
            if( c == EMPTY_CLUSTER ) c = *first_free_cluster.begin();
            assert( inCl[v] != c );
            moveNodeTo(v,c);
        }
        current_result += cnd.swpVal();
    }

    cluster_nodes.clear(); // clearing cluster nodes

    return improved;
}

template<class _T>
void NEG::localShuffle(_T &v) {
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

void NEG::resizeStructuresForEmptyCluster(int empty_cl) {
    while( cluster_weights.size() <= empty_cl ){
        cluster_weights.push_back(0);
        cluster_nodes.push_back( {} );
    }
}





