//
// Created by sylwester on 3/8/21.
//

#include <datastructures/FAU.h>
#include <graphs/GraphUtils.h>
#include <utils/RandomNumberGenerators.h>
#include <combinatorics/CombinatoricUtils.h>
#include <CONTESTS/PACE21/heur/PaceUtils.h>
#include <CONTESTS/PACE21/heur/EOCreators/ComponentExpansion.h>
#include "CollectionOperators.h"
#include "CONTESTS/PACE21/heur/State.h"

State::State(ClusterGraph &clg, StateInitializationType init_type) {
    this->clg = &clg;
    this->N = clg.N;
    initializeStateData(init_type);
    createClNeighGraph();

    { // creating hashes
        LL seed = 1'718'111'121;
        UniformIntGenerator rnd(0, 1ll * 1'000'000'000 * 1'000'000'000, seed);
        hashes = VLL(2*N);
        for(int i=0; i<hashes.size(); i++) hashes[i] = rnd.rand();
    }
}

void State::clearState(){
    fill(ALL(inCl), -1); fill(ALL(degInCl),0); fill(ALL(idInCl),-1);
    fill(ALL(was),false); fill(ALL(was2),false); fill(ALL(was2),false);
    fill(ALL(marker),0);
//    for( int i=0; i<N; i++ ) VI().swap(marker2[i]);
    clusters.clear(); cl_neigh_graph.clear();
}


void State::applySwap( VPII & to_swap ) {
    for( auto & cl : clusters ){
        for( int d : cl.g.nodes ) marker[d] = cl.id;
    }

    /**
     * If i >= getIdOfEmptyCluster() then empty_cl_mapper[i] is the id of the newly created cluster
     * empty_cl_mapper[i] is the i
     */
    unordered_map<int,int> empty_cl_mapper;
    int first_empty = getIdOfEmptyCluster();
    for( PII & p : to_swap ){
        if( p.second < getIdOfEmptyCluster() ) marker[p.first] = p.second;
        else{
            auto it = empty_cl_mapper.find( p.second );
            if( it == empty_cl_mapper.end() ){ // cluster p.second is NOT present - creating it
                if( first_empty == N ){
                    cerr << "In State::applySwap(): id of the first empty cluster first_empty = " << first_empty
                         << " >= " << N << " = N. This should not happen. If this happen, then it means that all "
                         "clusters are single-noded (e.g. just after State initialization) and some node is about "
                         "to be moved to a new cluster (what makes no difference, swpval=0). We do not need to make "
                         "that move." << endl << "IF YOU SEE THIS MESSAGE THEN IT MEANS, THAT PROBABLY SOMETHING SHOULD "
                         "BE MADE TO PREVENT FROM SUCH SITUATIONS - e.g. disallowing whole clusters to be moved to "
                         "an empty cluster" << endl;
                    continue;
                }
                empty_cl_mapper[p.second] = first_empty;
                marker[p.first] = first_empty;
                first_empty++;
            }else{ // cluster p.second is PRESENT, with id empty_cl_mapper[p.second]
                marker[p.first] = empty_cl_mapper[p.second];
            }
        }
    }

    VI first_occurrence(N,-1); // the first node that is to be in cluster i is at position first_occurence[i]


    VPII to_merge;
    to_merge.reserve(N);
    for( int i=0; i<N; i++ ){
        if( first_occurrence[ marker[i] ] == -1 ) first_occurrence[marker[i]] = i;
        else{
            to_merge.emplace_back( first_occurrence[marker[i]], i );
        }
    }

//    DEBUG(empty_cl_mapper);
//    DEBUG(first_occurrence);
//    DEBUG(to_merge);

    initializeStateData(SINGLE_NODES); // this resets all arrays. Now node i is in cluster i.
    mergeClusters(to_merge);
}

void State::initializeStateData(StateInitializationType init_type) {
    vector<Cluster>().swap(clusters);

    was = was2 = was3 = VB(2*N,false); // creating helper arrays
    marker = VI(N,false);
//    marker2 = VVI(N);

    inCl = VI(N, 0);
    idInCl = VI(N,0);
    degInCl = VI(N,0);

    int cnt = 0;
    if(init_type == SINGLE_NODES) {
        for (int c = 0; c < N; c++) {
            inCl[c] = cnt;
            clusters.emplace_back(*clg, VI({c}), cnt);
            cnt++;
        }
    }else if(init_type == RANDOM_MATCHING){
        vector<tuple<int,int,int>> edges;
        edges.reserve( GraphUtils::countEdges(clg->V) );
        int REPS = 1;
        VB was(N,false);
        vector<tuple<int,int,int>> matching;
        for (int i = 0; i < N; i++) {
            for (auto&[d, w] : clg->V[i]) if (d > i) edges.emplace_back(i, d, w);
        }

        while(REPS--) {
            if (!edges.empty()) {
                auto rnd = UniformIntGenerator(0, edges.size() - 1);
                shuffle(ALL(edges), rnd.getRNG());
            }

            fill(ALL(was),false);
            vector<tuple<int,int,int>> temp; temp.reserve(N);
            for( auto & [a,b,w] : edges ) {
                if(was[a] || was[b]) continue;
                was[a] = was[b] = true;
                temp.emplace_back(a,b,w);
            }
            if( temp.size() > matching.size() ){
                matching = temp;
//                clog << "Found matching of size: " << temp.size() << endl;
            }
        }

        fill(ALL(was),false);
//        for( auto & [a,b,w] : edges ){
        for( auto & [a,b,w] : matching ){
            if( !was[a] && !was[b] ){
                was[a] = was[b] = true;
                inCl[a] = inCl[b] = cnt;
                clusters.emplace_back(*clg, VI({a,b}), cnt);
                if( clusters.back().g.nodes[0] == a){ idInCl[a] = 0; idInCl[b] = 1; }
                else{
                    // this 'else' should not happen, unless inducing method is changed (now it preserves order of nodes
                    // passed in the constructor)
                    idInCl[a] = 1; idInCl[b] = 0;
                }
                degInCl[a] = degInCl[b] = w;
                cnt++;
            }
        }

        for( int i=0; i<N; i++ ){
            if(was[i]) continue;
            inCl[i] = cnt;
            clusters.emplace_back(*clg, VI({i}), cnt);
            cnt++;
        }
    }else if( init_type == SQRT_RANDOM ){
        clog << "SQRT RANDOM not tested yet" << endl;
        UniformIntGenerator rnd(0,1e9);
        VI perm = CombinatoricUtils::getRandomPermutation(N, rnd.rand());
        int sq = ceil(sqrt(N));

        VI current;

        auto add_current = [&](){
            clusters.emplace_back(*clg, current, cnt);
            for( int d : current ) was[d] = true;
            for(int d : current){
                inCl[d] = cnt;
                for( auto & [p,w] : clg->V[d] ){
                    if(d < p && was[p]){
                        degInCl[d] += w;
                        degInCl[p] += w;
                    }
                }
            }
            for( int d : current ) was[d] = false;

            for( int j=0; j<clusters.back().g.nodes.size(); j++ ) idInCl[ clusters.back().g.nodes[j] ] = j;

            current.clear();
            cnt++;
        };

        for( int i=0; i<N; i++ ){
            current.push_back( perm[i] );
            if( ((i+1) % sq) == 0 ){
                add_current();
            }
            if(!current.empty()) add_current();
        }

        clog << "FINISHED" << endl;
    }
    else if( init_type == MAXIMUM_MATCHING ){

    }else if( init_type == RANDOM_STATE_PERM ){
        VI perm = CombinatoricUtils::getRandomPermutation(N);

//        { // #TEST
//            sort( ALL(perm), [&](int a, int b){
//                return clg->V[a].size() < clg->V[b].size();
//            });
//        }

        VI ind_in_perm(N,0);
        for(int i=0; i<perm.size(); i++) ind_in_perm[perm[i]] = i;

        cnt = 0;
        VVI cl_nodes;
        VI cl_w(N+1,0), e_to_cl(N+1,0);

        for( int u : perm ){
            int nw_u = clg->node_weights[u];
            int ind_u = ind_in_perm[u];
            int tot_w = 0;

            for( auto & [v,w] : clg->V[u] ){
                if( ind_in_perm[v] > ind_u ) continue;
                tot_w += w;
                e_to_cl[ inCl[v] ] += w;
            }

            int best_val = tot_w, best_cl = cnt; // initially the best cluster is empty cluster

            for( auto & [v,w] : clg->V[u] ){
                if( ind_in_perm[v] > ind_u ) continue;

                int clv = inCl[v];
                int val = tot_w - e_to_cl[clv];
                val += nw_u * cl_w[clv] - e_to_cl[clv];
                if(val < best_val){
                    best_val = val;
                    best_cl = clv;
                }
            }

            { // move node u to cluster best_cl
                if(best_cl == cnt) cl_nodes.push_back(VI());

                idInCl[u] = cl_nodes[best_cl].size();
                cl_nodes[best_cl].push_back(u);
                inCl[u] = best_cl;
                cl_w[best_cl] += clg->node_weights[u];

                for( auto & [v,w] : clg->V[u] ){ // increasng degInCl
                    if( ind_in_perm[v] > ind_u ) continue;
                    if( inCl[u] == inCl[v] ){
                        degInCl[u] += w;
                        degInCl[v] += w;
                    }
                }
            }

            if(best_cl == cnt) cnt++; // increase number of clusters

            for( auto & [v,w] : clg->V[u] ) e_to_cl[ inCl[v] ] = 0; // clearing
        }

        for( int i=0; i<cl_nodes.size(); i++ ) clusters.emplace_back( *clg, cl_nodes[i],i );
    }
    else if( init_type == LEAF_TRIMMING ) sparseGraphTrimming();
    else if( init_type == EXPANSION_ORDER ){
        VI nodes(N); iota(ALL(nodes),0);
        Cluster cl( *clg, nodes,0 );
        ComponentExpansion ce(cl);
        ce.terminate_on_cluster_violator = true;

        VI perm = CombinatoricUtils::getRandomPermutation(N);
        sort( ALL(perm), [&]( int a, int b ){ return clg->V[a].size() < clg->V[b].size(); } );

        int cnt = 0;
        VVI cl_nodes;
        VB was(N,false);

        for( int v : perm ){
            if(was[v]) continue;

            auto ord = ce.getExpansionOrder( {v}, true, 1e9 );
//            DEBUG( ord.ord.size() );

            clusters.emplace_back( *clg, ord.ord, cnt );
            int cnt2 = 0;
            for( int d : ord.ord ){
                was[d] = true;
                inCl[d] = cnt;
                idInCl[d] = cnt2++;
            }

            for( int d : ord.ord ){
                for( auto & [p,w] : clusters.back().g.V[ idInCl[d] ] ){
                    degInCl[d] += w;
                }
            }

            cl_nodes.push_back(ord.ord);

            GraphUtils::removeNodesFromGraph( cl.g.V, ord.ord );
            cnt++;
        }
    }

    clusters.emplace_back( *clg, VI(), cnt++ ); // adding empty cluster

}

void State::mergeClusters(VPII &part) {
    int C = clusters.size();
    FAU fau( C );
    VB inPart(C,false);
    for( auto [a,b] : part ){
        fau.Union(a,b);
        inPart[a] = inPart[b] = true;
    }

    VVI newSets(C);
    for( int i=0; i<C-1; i++ ){ // C-1 is an empty cluster
        if( inPart[i] ){
            int c = fau.Find(i);
            newSets[c] += clusters[i].g.nodes;
        }
    }

    vector<Cluster> newClusters;

    for( int i=0; i<C-1; i++ ){ // inserting unmerged clusters, // C-1 is an empty cluster
        if( !inPart[i] ){
            newClusters.push_back( Cluster() ); // adding empty cluster
            swap( newClusters.back(), clusters[i] ); // swapping, not to copy every edge
        }
    }

    int invalid_id = -1;
    for(int i=0; i<C-1; i++){ // inserting new clusters after merging, // C-1 is an empty cluster
        if(!newSets[i].empty()){
            sort(ALL(newSets[i]));
//            newClusters.emplace_back( *clg, newSets[i], invalid_id );
            newClusters.emplace_back( *clg, newSets[i], invalid_id );
        }
    }

//    DEBUG(newClusters); ENDL(5);

    swap( clusters, newClusters );

    vector<Cluster>().swap(newClusters);

    C = clusters.size();
    for(int i=0; i<C; i++){
        clusters[i].id = i;
        for( int j=0; j<clusters[i].g.nodes.size(); j++ ){
            int d = clusters[i].g.nodes[j];
            inCl[d] = i;
            idInCl[ d ] = j;
            degInCl[d] = 0;
            for( auto& [p,w] : clusters[i].g.V[j] ) degInCl[d] += w;
        }
    }

    clusters.emplace_back( *clg, VI(), clusters.size() ); // adding empty cluster

    int zero_node_clusters = 0;
    for( auto & cl : clusters ) if(cl.g.nodes.empty()) zero_node_clusters++;
    assert(zero_node_clusters == 1);

    VVPII().swap(cl_neigh_graph);
}

void State::mergeClusters(VVI to_merge) {
    VPII temp;
    for(auto & v : to_merge ){
        for( int i=1; i<v.size(); i++ ) temp.emplace_back( v[0], v[i] );
    }
//    DEBUG(temp);
    mergeClusters(temp);
}


ostream& operator<<(ostream& str, State& st){
    str << "Clusters:" << endl;
//    for(auto cl : st.clusters) str << cl << endl;
    for( auto & cl : st.clusters ){
        clog << "[id: " << cl.id << ", weight: " << cl.cluster_weight << ", nodes: " << cl.g.nodes << "]" << endl;
    }
    str << "inCl: " << st.inCl << endl;
    str << "degInCl: " << st.degInCl << endl;
    str << "idInCl: " << st.idInCl << endl;
    return str;
}

void State::createClNeighGraph() {
    if( !cl_neigh_graph.empty() ) clog << "Creating cl_neigh_graph, even though it is already not empty!!" << endl;

    cl_neigh_graph = VVPII(N);
    VB was(N,false);
    VI neigh;
    VI weight(N,0);

    for( int i=0; i<N; i++ ){
        neigh.clear();
        for( auto & [p,w] : clg->V[i] ){
            int in_cl_p = inCl[p];
            weight[in_cl_p] += w;
            if( !was[in_cl_p] ){
                was[in_cl_p] = true;
                neigh.emplace_back(in_cl_p);
            }
        }

        // FIXME: possible optimization - we do not need to sort that - but tests will need to be changed to pass if
        // we remove that.
        sort(ALL(neigh), [&]( auto a, auto b ){
            if( weight[a] != weight[b] ) return  weight[a] > weight[b];
            else return a < b;
        });

        for( int d : neigh ){
            if( d != inCl[i] ) cl_neigh_graph[i].emplace_back( d, weight[d] );
        }
        for( int d : neigh ){
            was[d] = false;
            weight[d] = 0;
        }
    }
}

int State::calculateResultForState() {
    int res = 0;
    for( auto & cl : clusters ){
        for( auto d : cl.g.nodes ){
            res += ((cl.cluster_weight - clg->node_weights[d]) * clg->node_weights[d]);
            for( auto & [p,w] : clg->V[d] ){
                if( inCl[p] == inCl[d] ) res -= w; // this was calculated two lines above, but it is not a modification
                else res += w; // edges between clusters need to be deleted
            }
        }
    }

    res >>= 1; // each edges modification was calculated twice.
    return res;
}

void State::applyPartition(VI &part) {
    assert( part.size() == clg->N );

    initializeStateData(SINGLE_NODES);
    VVI to_merge = PaceUtils::partitionToClusters(part);
    mergeClusters(to_merge);
}

void State::sparseGraphTrimming() {

    VI deg(N,0);
    was = VB(N,false);
    deque<int> deg1_nodes;

    for(int i=0; i<N; i++){
        deg[i] = clg->V[i].size();
        if(deg[i] == 1) deg1_nodes.push_back(i);
    }

    int cnt = 0;
    VVI cluster_nodes;

    while( !deg1_nodes.empty() ){

    }







}
