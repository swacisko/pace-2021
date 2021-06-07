//
// Created by sylwester on 3/20/21.
//

#include <utils/RandomNumberGenerators.h>
#include <combinatorics/CombinatoricUtils.h>
#include <CONTESTS/PACE21/heur/EOCreators/ComponentExpansion.h>
#include <CONTESTS/PACE21/heur/Global.h>
#include "CONTESTS/PACE21/heur/SwapCandidates/SwpCndEO.h"
#include "CollectionOperators.h"

SwpCndEO::SwpCndEO(ExpansionOrder &eo, int k, int move_to, int swap_value, LL hash ) {
    this->eo = &eo;
    this->k = k;
    this->move_to = move_to;
    this->swap_value = swap_value;

    this->set_hash = hash;
}


//SwpCndEO::~SwpCndEO() {
//
//}

VPII SwpCndEO::getNodesToSwap() {
    VPII res(k+1);
    for( int i=0; i<=k; i++ ) res[i] = { eo->cl->g.nodes[ eo->ord[i] ], move_to };
    return res;
}

LL SwpCndEO::swpVal() {
    return swap_value;
}



//*********************************   SwpCndEOCreator   *************************************

SwpCndEOCreator::SwpCndEOCreator(State & s) : SwapCandidateCreatorAdapter(s) {
//    min_b_coef_for_cluster_size = VI(clg->N+1,0);
    min_b_coef_for_cluster_size = VI(clg->origV->size()+1,0); // weights of clusters
    createHashes();
}

void SwpCndEOCreator::createHashes() {
    hashes = VLL(clg->N);
    LL seed = 1'718'237'121;
    int SOME_SEED = 123892313;
    UniformIntGenerator rnd(0, 1ll * 1'000'000'000 * 1'000'000'000, SOME_SEED);
    for (int i = 0; i < clg->N; i++) hashes[i] = rnd.rand();
}


vector<SwapCandidate *> SwpCndEOCreator::createSwapCandidatesRaw() {
    return vector<SwapCandidate *>();
}

vector<SwpCndEO> SwpCndEOCreator::createSwapCandidates(ExpansionOrder & eo) {
    const bool debug = false;

    vector<SwpCndEO> res;

    Cluster* cl = eo.cl; // cluster for which eo was created
    VI ord = eo.ord;
    for( int & d : ord ) d = cl->g.nodes[d];
    // now ord should contain ids of nodes in the original cluster graph clg

//    fill(ALL(min_b_coef_for_cluster_size),0);
//    fill(ALL(edges_to_cluster),0);
//    fill(ALL(was),false);

    VI cluster_neigh;
    if( state->cl_neigh_graph.empty() ) state->createClNeighGraph(); // creating structure if necessary
    hull.clear();
    hull.insert_line( 0,0, state->getIdOfEmptyCluster() );
    sum_nw = 0;
    cut_value = 0;

    unordered_set<int> min_b_coef_cluster_sizes;

    LL current_hash = 0;

    for( int i=0; i<ord.size(); i++ ){
        int d = ord[i];
        int deg_in_cl_d = state->degInCl[d];
        int nw_d = clg->node_weights[d];
        sum_nw += nw_d;

        was[d] = true; // was marks elements in X_k
        current_hash ^= hashes[d];

        if(debug){
            clog << "Adding d: " << d << "  --> now X: " << StandardUtils::getSubarray(ord,0,i) << endl;
            DEBUG(sum_nw);
        }

        {
            int in_cl_d = state->inCl[d];
            for (auto &[cl, w] : state->cl_neigh_graph[d]) {
                if (in_cl_d == cl) continue; // we want p and d to be in different clusters
                if (edges_to_cluster[cl] == 0) cluster_neigh.push_back(cl);  // adding only those cluster neighbors that were not visited yet
                edges_to_cluster[cl] += w;

                // adding new linear functions to hull.
                if( keep_only_best_cluster_to_move_to ){
                    if(debug){
                        clog << "\tConsidering neighbor (cluster): " << PII(cl,w) << endl;
                        clog << "\tAdding to hull line " << state->clusters[cl].cluster_weight << " * x  - "
                             << 2 * edges_to_cluster[cl] << ", for id: " << cl << endl;
                        clog << "\thull.size(): " << hull.size() << endl;
                    }
                    int new_b_coef = -2 * edges_to_cluster[cl];
                    int cluster_weight = state->clusters[cl].cluster_weight;
                    if( new_b_coef < min_b_coef_for_cluster_size[cluster_weight] ){
                        min_b_coef_for_cluster_size[cluster_weight] = new_b_coef;
                        min_b_coef_cluster_sizes.insert(cluster_weight);

                        hull.insert_line( cluster_weight, new_b_coef, cl );
                    }
                }
            }
        }

        if(debug){
            DEBUG(cl->g.nodes);
            DEBUG(cl->g.V);
        }

        for( auto [p,w] : cl->g.V[ state->idInCl[d] ] ){
            int q = cl->g.nodes[p]; // mapping p to its original id, since we iterate over InducedClusterGraph cl->g.V (for speedup)
            if(was[q]){
                if(debug) clog << "\tDecreasing cut_value by " << w << "  (node " << q << ")" << endl;
                cut_value -= w;
            }
            else{
                if(debug) clog << "\tIncreasing cut_value by " << w << "  (node " << q << ")" << endl;
                cut_value += w;
            }
        }

        if(debug){
            DEBUG(cut_value);
        }

        VI neigh_cl(1, state->getIdOfEmptyCluster());
        if(keep_only_best_cluster_to_move_to){
            int best_choice = hull.get_line_id(sum_nw); // sum_nw is just 'the k' from unweighted version
            if(debug){
                DEBUG(sum_nw);
                DEBUG(best_choice);
            }
            neigh_cl.push_back(best_choice);
        }
        else neigh_cl += cluster_neigh;

        if(debug){
            DEBUG(neigh_cl);
        }

        int swpval1 = sum_nw * (cl->cluster_weight - sum_nw) - cut_value;

        createAllCandidatesForSwpCnd(swpval1, sum_nw, cut_value, neigh_cl, res,
                    [&]( int final_swpval, int trg_cl ){
                        return SwpCndEO( eo, i, trg_cl, final_swpval, current_hash ^ state->hashes[trg_cl] );
                        // we need to take current_hash ^ state->hashes[trg_cl] in order to distinguish between the
                        // same sets moved to different clusters
                    }
        );

        if(debug){
            DEBUG(res);
            ENDL(2);
        }
    }

    for(int d : ord) was[d] = false;

    for(int d : cluster_neigh)edges_to_cluster[d] = 0;
    for( int d : min_b_coef_cluster_sizes ) min_b_coef_for_cluster_size[d] = 0; // clearing

    for( int i=(int)res.size()-1; i>=0; i-- ){
        if( res[i].k + 1 == res[i].eo->ord.size() ){
            if( res[i].move_to == state->getIdOfEmptyCluster() ){
                // removing those candidates, that represent whole set moved to an empty cluster
                swap( res[i], res.back() );
                res.pop_back();
            }
        }
    }

    return res;
}


vector<SwpCndEO> SwpCndEOCreator::createSwapCandidates(vector<ExpansionOrder*> &exp_orders, int beg, int end) {
    if(end == -1) end = exp_orders.size();
    unordered_set<int> eo_hash_set;

    vector<SwpCndEO> res;
    for( int j=beg; j<end; j++ ){
        if( Global::checkTle() ) return res;

        auto & eo = exp_orders[j];
        auto temp = createSwapCandidates(*eo);
        for( int i=0; i<temp.size(); i++ ){
            SwpCndEO &cnd = temp[i];

            if( eo_hash_set.count( cnd.set_hash ) == 0 ){
                eo_hash_set.insert( cnd.set_hash );
//                res.push_back( std::move(cnd) ); // original
                res.push_back( cnd );
            }
        }
    }

    return res;
}

vector<ExpansionOrder> SwpCndEOCreator::createExpansionOrders( function<VVI(Cluster*)> fun ) {
    vector<ExpansionOrder> orders;

    UniformIntGenerator rnd(0,100);

    for( auto & cl : state->clusters ){
        if( cl.size() == 0 ) continue; // do not process empty cluster
        VVI sets = fun(&cl);
        ComponentExpansion ce( cl );

        for( auto & A : sets ){

//            int cnt = 0; // original version
            int cnt = rnd.nextInt(5); // #TEST - checking different expansion orders cmp rules

            // some possibilities for expansion order rules
            if(cnt == 0) ce.setCmpRules( VI({ 11,10,3,2,1 }) ); // originally only this was done
            else if(cnt == 1) ce.setCmpRules( VI({ 10,11,3,2,1 }) );
            else if(cnt == 2) ce.setCmpRules( VI({ 3,2,1 }) );
            else if(cnt == 3) ce.setCmpRules( VI({ 2,1,10 }) );
            else if(cnt == 4) ce.setCmpRules( VI({ 1,2,10 }) );

            auto ord = ce.getExpansionOrder(A,false); // do not use heap - we use CE for cluster
            orders.push_back(ord);
        }
    }


    return orders;
}

ostream& operator<<(ostream& str, SwpCndEO& cnd){
    str << "(nodes --> move_to: " << StandardUtils::unzip( cnd.getNodesToSwap() ).first << "  -->  " << cnd.move_to
        << ", swp_val: " << cnd.swpVal() /*<< ", set_hash: " << cnd.set_hash*/ << ")";
    return str;
}

VI SwpCndEO::getAffectedClusters(State& st) {
//    clog << "Calling getAffectedClusters for SwpCndEO - polimorphism works" << endl;
//    clog << "SwpCndEOCandidate: " << getNodesToSwap() << ", swpval: " << swpVal() << endl;

    return {eo->cl->id, move_to};
}
