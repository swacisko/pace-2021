//
// Created by sylwester on 4/8/21.
//

#include <CONTESTS/PACE21/heur/SwapCandidates/ComponentExpansionRepulsion.h>
#include <CONTESTS/PACE21/heur/Global.h>

ExpansionOrderRepulsion::ExpansionOrderRepulsion(Cluster &cl, VI v, VI to) {
    this->cl = &cl;
    ord = v;
    move_to = to;
}

ostream& operator<<( ostream& str, SwpCndEORepulsion& cnd ){
    using namespace StandardUtils;
    auto sub1 = getSubarray( cnd.eo->ord, 0, cnd.k );
    auto sub2 =  getSubarray(cnd.eo->move_to,0,cnd.k);
    VPII moves = zip( sub1, sub2 );
    str << "{swaps: " << moves << ", swpval: " << cnd.swpVal() << "}";
    return str;
}

vector<ExpansionOrder>
ExpansionOrderRepulsion::induceOrders(State &st, ExpansionOrderRepulsion &eo) {
    vector<ExpansionOrder> res;
    unordered_map<int, VI, fib_hash> zb;
    for( int i=0; i<eo.ord.size(); i++ ){
        int v = eo.ord[i];
        int c = eo.move_to[i];
        zb[c].push_back(v);
    }

    for( auto & [c,vec] : zb ){
        res.emplace_back( vec, eo.cl ); // All induced expansion orders are created for cluster [eo.cl]
    }

    for( auto & o : res ){
        for( int & d : o.ord ) d = st.idInCl[d]; // mapping to ids in cluster
    }

    return res;
}

SwpCndEORepulsion::SwpCndEORepulsion( ExpansionOrderRepulsion * eo, int k, int swap_value, LL hash ){
    this->eo = eo;
    this->k = k;
    this->swap_value = swap_value;
    this->set_hash = hash;
}

VPII SwpCndEORepulsion::getNodesToSwap(){
    VI nodes = StandardUtils::getSubarray(eo->ord,0,k);
    VI move_to = StandardUtils::getSubarray(eo->move_to,0,k);
    return StandardUtils::zip(nodes, move_to);
}


ComponentExpansionRepulsion::ComponentExpansionRepulsion( State & st ){
    this->st = &st;
    clg = st.clg;

    V = &clg->V;
    N = V->size();
    eInS = VI(N, 0);

    inS = VB(N,false);
    sumNWinS = 0;

    was = was2 = was3 = VB(2*N,false);

    // we need to have more possible cluster to move to, since many nodes may be moved to 'new empty clusters'
    cluster_weights = VI(2*N ,0);

    edges_to_cluster.resize( N );
    cluster_sets.resize(N);
    swpval_core_for_node.resize(N);
}

pair<vector<ExpansionOrderRepulsion *>, vector<SwpCndEORepulsion>> ComponentExpansionRepulsion::createSwapCandidates() {
    pair<vector<ExpansionOrderRepulsion *>, vector<SwpCndEORepulsion>> res;

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


pair<ExpansionOrderRepulsion*, vector<SwpCndEORepulsion> >
ComponentExpansionRepulsion::createSwapCandidates(Cluster &cl) {
    const bool debug = false;
    this->cl = &cl;

    assert( !cl.g.nodes.empty() );

    // cl_neigh denotes the ids of clusters that are neighbors of cluster [cl]. [was] marks this property.
    cl_neigh = VI(1, st->getIdOfEmptyCluster());
    was[cl_neigh[0]] = true;

    auto debug_all = [&](){
        DEBUG(S); DEBUG(sumNWinS); DEBUG(eInS);
        DEBUG(cl_neigh); DEBUG(cluster_weights); DEBUG(edges_to_cluster);
        DEBUG(swpval_core_for_node);
        clog << "Cluster sets: " << endl;
        for( int i=0; i<cluster_sets.size(); i++ ){
            if( !cluster_sets[i].empty() ) clog << "cluster_sets[" << i << "]: " << cluster_sets[i] << endl;
        }
    };

    if(debug) debug_all();

    if(st->cl_neigh_graph.empty()) st->createClNeighGraph();

    { // initializing some
        sumNWinS = cl.cluster_weight;

        S.insert(ALL(cl.g.nodes));
        for( int d : S ){
            inS[d] = true;
            eInS[d] = st->degInCl[d];
        }

        // S, inS, sumEW, cluster_weights, eInS
        // should be already cleared - clearing at the end of each call to createSwapCandidates

        for( int d : cl.g.nodes ){
            edges_to_cluster[d].insert(ALL(st->cl_neigh_graph[d]));
            for( auto & [c,w] : st->cl_neigh_graph[d] ){
                if(!was[c]){
                    was[c] = true; // [was] marks cluster neighbors of nodes from [cl]
                    cl_neigh.push_back(c);
                }
            }
        }
        for( int c : cl_neigh ) cluster_weights[c] = st->clusters[c].cluster_weight;

        if(debug) debug_all();

        for( int d : S ){
            updateSwpValCoreForNode(d, st->getIdOfEmptyCluster());
            for( auto [c,w] : st->cl_neigh_graph[d] ){
                updateSwpValCoreForNode(d, c);
            }
        }
    }

    if(debug) debug_all();

    VI ord, swpvals, move_to;
    int swpval = 0;

    if(debug) clog << endl << "Proceeding to repulsion loop" << endl;

    int first_left_empty_cluster = st->getIdOfEmptyCluster()+1;
    while( !S.empty() ){
        if(debug) clog << endl << "NEXT ITERATION" << endl;

        PII choice = getNodeToMove();
        int node_to_move = choice.first;
        int to = choice.second;

        if( to == st->getIdOfEmptyCluster() ){
            // if a node should be moved to an empty cluster, a new one is created instead and node is moved to the new one
            if(debug){
                clog << "should move node " << node_to_move << " to empty cluster "
                    << st->getIdOfEmptyCluster() << ". Moving it to newly created cluster "
                    << first_left_empty_cluster << " instead" << endl;
            }

            assert( !was[first_left_empty_cluster] );
            was[first_left_empty_cluster] = true;
            cl_neigh.push_back(first_left_empty_cluster);
            to = first_left_empty_cluster;

            first_left_empty_cluster++;
        }

        swpval += calculateSwpVal(node_to_move, to);
        swpvals.push_back(swpval);

        ord.push_back(node_to_move);

        if(debug){
            clog << "Moving node " << node_to_move << " to cluster " << to << " with swpval " << calculateSwpVal(
                    node_to_move, to)
                 << ". Total swpval: " << swpval << endl;
        }

        move_to.push_back(to);

        moveNodeTo( node_to_move, to );
        if(debug) debug_all();
        for( int u : S ) updateSwpValCoreForNode(u, to);

        if(debug){
            debug_all();
            clog << endl << endl;
        }
    }

    pair<ExpansionOrderRepulsion*, vector<SwpCndEORepulsion> > res;
    { // creating expansion order and swap candidates
        ExpansionOrderRepulsion* eor = new ExpansionOrderRepulsion( cl, ord, move_to );
        res.first = eor;
        assert( ord.size() == swpvals.size() );
        assert( ord.size() == move_to.size() );
        LL hash = 0;

        for( int i=0; i<ord.size(); i++ ){
            hash ^= ( st->hashes[ord[i]] + st->hashes[move_to[i]] );
            if( !keep_only_nonpositive_candidates || (swpvals[i] <= 0) ){
                res.second.emplace_back( eor, i, swpvals[i], hash );
            }
        }

        if(!res.second.empty() && res.second.back().k == ord.size()-1){
            // in this if we check whether the last swap candidate represents the whole set. If it does and all nodes
            // are moved to the same empty cluster, then we remove that candidate from the set of candidates.
            unordered_set<int> zb;
            bool ok = true;
            for( int i=0; i<ord.size(); i++ ){
                zb.insert(ord[i]);
                if( zb.size() > 1 ){ ok = false; break; }
            }
            if(ok) res.second.pop_back();
        }
    }

    if(debug) DEBUG(res.second);

    assert(S.empty());

    S.insert(ALL(cl.g.nodes));
    { // clearing section
        for (int c : cl_neigh) was[c] = false;
        for (int d : S){
            eInS[d] = 0;
            inS[d] = false;
            unordered_map<int,int>().swap(edges_to_cluster[d]); // clearing edges_to_cluster
            unordered_map<int,int>().swap(swpval_core_for_node[d]); // clearing edges_to_cluster
            set<ClusterSet>().swap(cluster_sets[d]); // clearing edges_to_cluster
        }
        for( int c : cl_neigh ) cluster_weights[c] = 0;
        S.clear();
        sumNWinS = 0;
        cl_neigh.clear();
    }

    return res;
}

PII ComponentExpansionRepulsion::getNodeToMove() {
    const bool debug = false;
    int best_swpval = 1e9;
    int best_v = -1, best_c = -1;
    for( int u : S ){
        int c = cluster_sets[u].begin()->c;

        int val = calculateSwpVal(u, c);

        if(debug){
            clog << "Best choice for node u: " << u << ": (cluster " << c << ", core: " << swpval_core_for_node[u][c]
                 << ", swpval: " << val << ")" << endl;
        }

        if(val < best_swpval){
            best_swpval = val;
            best_v = u;
            best_c = c;
        }
    }

    if(debug) clog << endl;

    assert(best_c != -1);
    return { best_v, best_c };
}


void ComponentExpansionRepulsion::moveNodeTo(int v, int c) {
    cluster_weights[c] += clg->node_weights[v];
    sumNWinS -= clg->node_weights[v];
    inS[v] = false;
    S.erase(v);
//    for( auto [d,w] : cl->g.V[v] ){ // old buggy version
    for( auto [d,w] : cl->g.V[ st->idInCl[v] ] ){ // new, probably correct version
        d = cl->g.nodes[d];
        if( inS[d] ){
            eInS[d] -= w;
            edges_to_cluster[d][c] += w;
        }
    }

    eInS[v] = 0; // we can clear that now, it will not be needed anymore
    edges_to_cluster[v].clear(); // we do not need that either
    swpval_core_for_node[v].clear(); // we do not need that either
    cluster_sets[v].clear();  // we do not need that either
}


void ComponentExpansionRepulsion::updateSwpValCoreForNode(int u, int c) {
    const bool debug = false;
    if(debug){
        clog << "u: " << u << ", c: " << c << endl;
        clog << "Before update, swpval_core_for_node[u][c]: " << swpval_core_for_node[u][c] << endl;
    }

    if(debug){
        clog << "erasing (u,c) from cluster_sets[u]:" << cluster_sets[u] << endl;
    }
    cluster_sets[u].erase( ClusterSet( this,u,c ) );


    if(debug) clog << "after erasing  cluster_sets[u]:" << cluster_sets[u] << endl;

    int nw_u = clg->node_weights[u];
    int Ec = edges_to_cluster[u][c];
    int cl_weight = cluster_weights[c];
    int core_value = nw_u * cl_weight - (Ec << 1);

    if(debug) clog << "core_value: " << core_value << endl;

    swpval_core_for_node[u][c] = core_value;

    cluster_sets[u].insert(ClusterSet( this,u,c ));

    if(debug){
        clog << "After update, swpval_core_for_node[u][c]: " << swpval_core_for_node[u][c] << " and calculateSwpVal: "
             << calculateSwpVal(u, c) << endl;
    }
}

int ComponentExpansionRepulsion::calculateSwpVal(int u, int c) {
    const bool debug = false;

    int nw_u = clg->node_weights[u];
    int nw_S = sumNWinS;
    int eu = eInS[u];
    auto it = edges_to_cluster[u].find(c);
    int Ec = ( (it == edges_to_cluster[u].end()) ? 0 : it->second  );
    int cl_weight = cluster_weights[c];

    int before = nw_u * ( nw_S - nw_u ) - eu + Ec;
    int after = nw_u * cl_weight - Ec + eu;

    if(debug) clog << " (bef:" << before <<",aft:" << after << ") ";
    int swpval = after - before;

    return swpval;
}

ostream& operator<<(ostream& str, ComponentExpansionRepulsion::ClusterSet& cs){
    str << "(v,c,core): " << "(" << cs.v << "," << cs.c << ","
         << cs.cer->swpval_core_for_node[cs.v][cs.c] << ")";
    return str;
}