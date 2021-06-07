//
// Created by sylwester on 4/8/21.
//

#include <CONTESTS/PACE21/heur/SwapCandidates/ComponentExpansionAttraction.h>
#include <CONTESTS/PACE21/heur/Global.h>

ostream& operator<<( ostream& str, SwpCndEOAttraction& cnd ){
    using namespace StandardUtils;
    VPII moves = cnd.getNodesToSwap();
    str << "{swaps: " << moves << ", swpval: " << cnd.swpVal() << "}";
    return str;
}

vector<ExpansionOrder>
ExpansionOrderAttraction::induceOrders(State & st, ExpansionOrderAttraction & eo ){
    vector<ExpansionOrder> res;
    unordered_map<int, VI, fib_hash> zb;
    for( int i=0; i<eo.ord.size(); i++ ){
        int v = eo.ord[i];
        int c = st.inCl[v];
        zb[c].push_back(v);
    }

    for( auto & [c,vec] : zb ){
        res.emplace_back( vec, &st.clusters[c] );
    }

    for( auto & o : res ){
        for( int & d : o.ord ) d = st.idInCl[d]; // mapping to ids in cluster
    }
    return res;
}


//***********************************  ComponentExpansionAttraction  *****************************************
ComponentExpansionAttraction::ComponentExpansionAttraction(State &st) {
    this->st = &st;
    clg = st.clg;

    V = &clg->V;
    N = V->size();
    eToS = VI(N, 0);

    inS = VB(N,false);
    sumNWinS = 0;

    was = was2 = was3 = VB( 2*N,false);
    cluster_weights = VI(st.clusters.size() ,0); // here it is enough to have cluster.size() clusters

    int max_node_weight = 0;
    for( int nw : clg->node_weights ) max_node_weight = max(max_node_weight, nw);
    cluster_cores.resize( max_node_weight + 1 );

    edges_in_cl = VI(N,0);
}

void ComponentExpansionAttraction::initialize() {
    // cl_neigh denotes the ids of clusters that are neighbors of cluster [cl]. [was] marks this property.
    cl_neigh.clear(); // here we do not need an empty cluster - this is just cluster-neighborhood of S

    { // some initialization
        S.clear();
        S.insert(ALL(cl->g.nodes));
        for(int d : S) inS[d] = true;
        sumNWinS = 0;

        for( int d : S){
            sumNWinS += clg->node_weights[d];
            for( auto [p,w] : clg->V[d] ){
                int cl_p = st->inCl[p];
                if( cl_p != cl->id ){
                    eToS[p] += w;

                    if( !was[cl_p] ){
                        was[cl_p] = true;
                        cl_neigh.push_back(cl_p);
                    }
                }
            }
        }

        // cl_neigh should be now created, creating cluster_weights and cluster cores.
        for( int c : cl_neigh ){
            cluster_weights[c] = st->clusters[c].cluster_weight;

            for( int d : st->clusters[c].g.nodes ){
                edges_in_cl[d] = st->degInCl[d]; // initially that is the number of edges of node d in its cluster

                int nw_d =  clg->node_weights[d];

                int free_value = getFreeValue(d);
                cluster_cores[nw_d].insert( PII( free_value,d ) );

                node_neigh_weights[nw_d]++;
            }
        }

    }
}

pair<vector<ExpansionOrderAttraction *>, vector<SwpCndEOAttraction>>
ComponentExpansionAttraction::createSwapCandidates() {
    pair<vector<ExpansionOrderAttraction *>, vector<SwpCndEOAttraction>> res;

    for( auto & cl : st->clusters ){
        if( Global::checkTle() ) return res;
        if( cl.size() < min_cluster_size ) continue;
        if( cl.id == st->getIdOfEmptyCluster() ) continue; // we do not process empty cluster

        auto temp = createSwapCandidates(cl);
        res.first.push_back(temp.first);
        res.second += temp.second;
    }

    return res;
}


pair<ExpansionOrderAttraction*, vector<SwpCndEOAttraction> >
ComponentExpansionAttraction::createSwapCandidates(Cluster& cl){
    const bool debug = false;
    this->cl = &cl;

    auto debug_all = [&](){
        DEBUG(S); DEBUG(sumNWinS); DEBUG(eToS);
        DEBUG(cl_neigh); DEBUG(cluster_weights); //DEBUG(cluster_cores);
        DEBUG(edges_in_cl); DEBUG(node_neigh_weights);
        for( auto [c,cnt] : node_neigh_weights ) clog << "******   cluster_cores[" << c << "]: " << endl << cluster_cores[c] << endl;
    };

    initialize();

    if(debug) debug_all();

    // total number of nodes (not node weights) in S and all neighbors of S
    int total_size = S.size() + accumulate( ALL(cl_neigh), 0,
                 [&](int s, int c){ return s + st->clusters[c].size(); } );

    if(debug){
        DEBUG(total_size);
        clog << endl << "***************************  Proceeding to iterations  ***********" << endl;
    }

    VI ord, swpvals;
    int total_swpval = 0;


    /**
     * After this number of nodes is attracted to S, we terminate. This is due to the fact, that rarely happens that
     * we find a swap candidate with greater number of nodes
     */
    const int MAX_NODES_TO_MOVE = 15; // #TEST #CAUTION - hardcoded small value for tests
    int nodes_moved = 0;

    while( S.size() < total_size ){
        if( debug ){ // just an assertion
            // checking whether swap values and equal to those of nw_d * sumNWinS + free_value
            for( auto & zb : cluster_cores ){
                for( PII p : zb ){
                    int d = p.second;
                    int free_val = p.first;
                    assert( free_val == getFreeValue(d) );
                }
            }
        }

        if(debug) clog << "---> NEXT ITERATION" << endl;

        int node_to_move = getNodeToMove();
        if(debug){
            DEBUG(node_to_move);
//            DEBUG(getFreeValue(node_to_move));
        }

        int swpval = clg->node_weights[node_to_move] * sumNWinS + getFreeValue(node_to_move);
        if(debug) DEBUG(swpval);

        total_swpval += swpval;
        if(debug) DEBUG(total_swpval);

        modifyNodeEntryInClusterCores(node_to_move,false);

        if(debug){
            clog << "After removing neighbors and cluster-mates of " << node_to_move << ":" << endl;
            for( int c : cl_neigh ) clog << "******   cluster_cores[" << c << "]: " << endl << cluster_cores[c] << endl;
        }

        moveNodeToS(node_to_move);

        modifyNodeEntryInClusterCores(node_to_move,true);

        {
            ord.push_back(node_to_move);
            swpvals.push_back(total_swpval);
        }

        if(debug){
            clog << "After moving node" << endl;
            debug_all();
        }

        if(debug) ENDL(5);

        nodes_moved++;
        if( nodes_moved == MAX_NODES_TO_MOVE ) break;
    }

    pair<ExpansionOrderAttraction*, vector<SwpCndEOAttraction> > res;
    {
        ExpansionOrderAttraction *eoa = new ExpansionOrderAttraction(cl, ord);
        res.first = eoa;
        res.second.reserve(ord.size());
        for (int i = 0; i < ord.size(); i++) {
            if( !keep_only_nonpositive_candidates || (swpvals[i] <= 0) ){
                res.second.emplace_back(*eoa, i, swpvals[i]);
            }
        }
    }

    if(debug) DEBUG(res);

    { // clearing section
        sumNWinS = 0;
        S.clear();
        S.insert(ALL(cl.g.nodes));
        for( auto & [nw,cnt] : node_neigh_weights ) set<PII>().swap(cluster_cores[nw]);
        node_neigh_weights.clear();
        for( int d : S ){
            inS[d] = false;
            eToS[d] = 0;
        }

        for( int c : cl_neigh ){
            cluster_weights[c] = 0;
            was[c] = false;

            for( int d : st->clusters[c].g.nodes ){
                inS[d] = false;
                eToS[d] = 0;
                edges_in_cl[d] = 0;
            }
        }
        cl_neigh.clear();

        S.clear();
    }

    return res;
}

void ComponentExpansionAttraction::modifyNodeEntryInClusterCores( int node_to_move, bool insert ){
    for( auto & [d,w] : clg->V[node_to_move] ){
        int cl = st->inCl[d];

        // if node d is in some cluster from cl_neigh and it is not in S nor in the same cluster as node_to_move,
        // then it was not moved yet to S, so its 'free value' needs to be updated - inserted/erased into/from
        // [cluster_cores]
        if( was[cl] && !inS[d] && cl != st->inCl[node_to_move] ){
            int nw_d = clg->node_weights[d];
            int free_value = getFreeValue(d);

            if(insert){
                cluster_cores[nw_d].insert( PII( free_value,d ) );
            }
            else{
                cluster_cores[nw_d].erase( PII( free_value,d ) );
            }
        }
    }

    // for each node in the cluster that contains/contained node_to_move
    for( int d : st->clusters[ st->inCl[node_to_move] ].g.nodes ){
        if( d == node_to_move ) continue;
        if( !inS[d] ){
            int nw_d = clg->node_weights[d];
            int free_value = getFreeValue(d);

            if(insert){
                cluster_cores[nw_d].insert( PII( free_value,d ) );
            }
            else{
                cluster_cores[nw_d].erase( PII( free_value,d ) );
            }
        }
    }

    if(!insert){ // removing node_to_move from cluster_cores
        cluster_cores[ clg->node_weights[node_to_move] ].erase( PII( getFreeValue(node_to_move), node_to_move ) );
    }

}

int ComponentExpansionAttraction::getNodeToMove() {
    const bool debug = false;
    int best_u = -1;
    int best_val = Constants::INF;

    for( auto & [n_w,cnt] : node_neigh_weights ){
        auto it = cluster_cores[n_w].begin();
        int d = it->second;
        int free_val = it->first;
        int swpval = n_w * sumNWinS + free_val;

        if(debug){
            clog << "n_w: " << n_w << ", d: " << d << ", swpval: " << swpval << endl;
        }

        if(swpval < best_val){
            best_val = swpval;
            best_u = d;
        }
    }
    return best_u;
}

void ComponentExpansionAttraction::moveNodeToS(int u) {
    int nw_u = clg->node_weights[u];
    cluster_weights[ st->inCl[u] ] -= nw_u;
    sumNWinS += nw_u;
    S.insert(u);
    inS[u] = true;

    { // modifying node_neigh_weights
        auto it = node_neigh_weights.find(nw_u);
        assert(it != node_neigh_weights.end());

        it->second--;
        if (it->second == 0) node_neigh_weights.erase(it);
    }

    for( auto & [d,w] : clg->V[u] ){
        int cl_d = st->inCl[d];
        if( was[cl_d] && !inS[d] ){
            eToS[d] += w;
            if( cl_d == st->inCl[u] ) edges_in_cl[d] -= w;
        }
    }
}


int ComponentExpansionAttraction::getFreeValue(int d) {
    const bool debug = false;
    int nw_d = clg->node_weights[d];
    int free_value = -nw_d * cluster_weights[st->inCl[d]] + nw_d * nw_d + 2*edges_in_cl[d] - 2*eToS[d];

    if(debug){
        DEBUG(d);DEBUG(nw_d);
        DEBUG( cluster_weights[st->inCl[d]] );
        DEBUG(2*edges_in_cl[d]);
        DEBUG(2*eToS[d]);

        clog << "free_value(" << d << "): " << free_value << endl;
    }

    return free_value;
}



