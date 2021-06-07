//
// Created by sylwester on 5/18/21.
//

#include <CONTESTS/PACE21/heur/StateImprovers/SparseGraphTrimmer.h>
#include <datastructures/Heap.h>
#include <CONTESTS/PACE21/heur/PaceUtils.h>
#include <graphs/GraphTrimmer.h>
#include <graphs/GraphUtils.h>

SparseGraphTrimmer::SparseGraphTrimmer(State *st) {
    this->st = st;
    clg = st->clg;
}

State &SparseGraphTrimmer::createState() {
//    contractNodesNondescendingDegree();
    greedyPathClustering();

    return *st;
}

void SparseGraphTrimmer::contractNodesNondescendingDegree() {
    const bool debug = false;

    int N = st->N;

    VI deg(N,0);
    for( int i=0; i<N; i++ ) deg[i] = clg->V[i].size();


    auto cmp = [&](PII a, PII b){
        if(a.first != b.first) return a.first < b.first;
        else return a.second < b.second;
    };

    set<PII, decltype(cmp)> zb(cmp);
    for(int i=0; i<N; i++) zb.insert(PII(deg[i],i));

//    if(debug) DEBUG(zb);


    VB was(N,false);

    int cnt = 0;
    VI partition(N,-1);

    auto contract = [&](int v){
        if(debug) clog << "Contracting node " << v << ", cnt = " << cnt << endl;

        partition[v] = cnt;
        for( auto & [u,w] : clg->V[v] ){
            if(was[u]) continue;

            if(debug) clog << "Marking node " << u << " to partition " << cnt << endl;
            partition[u] = cnt;

            for( auto & [d,ww] : clg->V[u] ){
                if(was[d]) continue;
                zb.erase( PII(deg[d],d) );
                deg[d]--;
                zb.insert(PII(deg[d],d));
            }
        }

        was[v] = true;
        for( auto & [u,w] : clg->V[v] ){
            was[u] = true;
            zb.erase(PII(deg[u],u));
        }

        zb.erase(PII(deg[v],v));
        cnt++;
    };

    while( !zb.empty() ){
        int v = zb.begin()->second;
//        clog << "\rdeg[v] = " << deg[v] << flush;

        zb.erase(zb.begin());

        contract(v);
    }

    st->applyPartition(partition);

    VI part = PaceUtils::mapClgPartitionToOriginalPartition(*clg, partition);
    DEBUG(PaceUtils::evaluateSolution( *clg->origV, part));
    exit(1);
}

void SparseGraphTrimmer::greedyPathClustering() {
    VVI V = *clg->origV;
    GraphTrimmer trm(V);
    auto paths = trm.findMaximalPathsAndCycles();

    int cnt = 0;
    int N = V.size();
    VI partition(N,-1);

    while( !paths.empty() ){
        for( int i=0; i<paths.size(); i++ ){
//            if( paths[i].size() == 1 ) assert( V[paths[i][0]].size() > 0 );

            int first = paths[i][0];
            if( V[first].size() > 1 ){
                int a = -1;
                for( int d : V[first] ) if( d != paths[i][1] ){ a = d; break; }
                if(a != -1) paths[i].insert(paths[i].begin(),a);
            }

            int last = paths[i].back();
            if( V[last].size() > 1 ){
                int a = -1;
                for( int d : V[last] ) if( paths[i].size() == 1 || d != paths[i][ paths[i].size()-2 ] ){ a = d; break; }
                if(a != -1) paths[i].push_back(a);
            }


            VI temp;
            for( int j=0; j<paths[i].size(); j++ ){

                if( temp.size() == 3 ){
                    for( int x=0; x<temp.size(); x++ ) partition[ temp[x] ] = cnt;
                    cnt++;
                    GraphUtils::removeNodesFromGraph(V,temp);
                    temp.clear();
                }
                temp.push_back(paths[i][j]);
            }

            if( !temp.empty() ){
                for( int x=0; x<temp.size(); x++ ) partition[ temp[x] ] = cnt;
                cnt++;
                GraphUtils::removeNodesFromGraph(V,temp);
                temp.clear();
            }
        }

        int E = GraphUtils::countEdges(V);
//        DEBUG( E );

        GraphTrimmer trm(V);
        paths = trm.findMaximalPathsAndCycles();

        {
//            DEBUG(paths.size());

            int L = accumulate(ALL(paths), 0, [](int s, VI &v) { return s + v.size(); });
//            DEBUG(L);

            int L2 = 0;
            for( int i=0; i<paths.size(); i++ ){
                int a = paths[i][0];
                int b = paths[i].back();
                if(V[a].size() >= 2 && V[b].size() >= 2) L2++;
            }
//            DEBUG(L2);
        }
    }

    for( int i=0; i<N; i++ ) if( partition[i] == -1 ) partition[i] = cnt++;

    st->applyPartition(partition);

//    DEBUG(PaceUtils::evaluateSolution( *clg->origV, partition));
//    exit(1);
}

void SparseGraphTrimmer::cherryRemover() {
    map<int,VI> zb;
    VVI V = *clg->origV;
    int N = V.size();

    VI partition(N);
    int cnt = 0;
    iota(ALL(partition),0);

    do{
        zb.clear();
        for(int i=0; i<N; i++){
            if( V[i].size() == 1 ){
                int d = V[i][0];
                zb[d].push_back(i);
            }
        }

        int cherries = 0;
        for( auto & [v,vec] : zb ){
            if(vec.size() >= 2){
                int a = vec[0], b = vec[1];
                partition[a] = partition[b] = partition[v] = cnt;
                cnt++;
                VI temp = {a,b,v};
                GraphUtils::removeNodesFromGraph(V,temp);

                cherries++;
            }
        }
        DEBUG(cherries);

        if(cherries == 0) break;


    }while( !zb.empty() );

    for( int i=0; i<N; i++ ) if( partition[i] == -1 ) partition[i] = cnt++;

    st->applyPartition(partition);
    DEBUG(PaceUtils::evaluateSolution( *clg->origV, partition));

    exit(2);
}

