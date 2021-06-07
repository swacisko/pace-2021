//
// Created by sylwester on 8/8/19.
//


#include <graphs/GraphUtils.h>
#include <graphs/trees/Tree.h>
#include <graphs/components/ConnectedComponents.h>

#include "graphs/GraphUtils.h"
#include "combinatorics/CombinatoricUtils.h"
//#include <graphs/VertexCover/approximation/PartitionSVC.h>
//#include <graphs/VertexCover/SolutionVC.h>
#include <utils/StandardUtils.h>
#include "graphs/cliques/CliqueExtension.h"
#include "graphs/components/BridgesAndArtPoints.h"

VI GraphUtils::getComplimentaryNodes( VVI & V, VI & nodes ){
    VB inNodes( V.size(),false );
    for(auto p : nodes) inNodes[p] = true;
    VI res;
    for( int i=0; i<V.size(); i++ ){
        if( !inNodes[i] ) res.push_back(i);
    }
    return res;
}

VI GraphUtils::getRandomColoring( VVI & V ){
    VI color(V.size(),-1);
    VI perm = CombinatoricUtils::getRandomPermutation( V.size() );

//    cerr << "Radnom permutation: " << perm << endl;

    VB colorTaken(V.size()+1);
    for( int p : perm ){
        for( int d : V[p] ){ // marking colors taken by neighbors
            int c = color[d];
            if( c != -1 ) colorTaken[c] = true;
        }

        for( int i=0; i<colorTaken.size(); i++ ){ // choosing first free color and assigning that
            if( !colorTaken[i] ){
                color[p] = i;
//                cerr << "assigning color " << i << " to " << p << endl;
                break;
            }
        }

        for( int d : V[p] ){ // resetting taken colors
            int c = color[d];
            if( c != -1 ) colorTaken[c] = false;
        }
    }
    return color;
}

VPII GraphUtils::getGraphEdges( VVI & V ){
    VPII res;
    res.reserve( countEdges(V) );
    for( int i=0; i<V.size(); i++ ){
        for( int d : V[i] ){
            if(d>i) res.push_back( {i,d} );
        }
    }
    return res;
}

VI GraphUtils::getNeighborhood(VVI &V, VI &S, bool useSet) {
    VI res;
    if( useSet ){
        unordered_set<int> zb;
        for( int p : S ){
            for(int d : V[p]){
                zb.insert(d);
            }
        }
        res = VI(ALL(zb));
    }else{
        VB was(V.size(),false);
        for( int p : S ){
            for(int d : V[p]){
                if( !was[d] ){
                    was[d] = true;
                    res.push_back(d);
                }
            }
        }
    }
    return res;
}

VI GraphUtils::getNeighborhoodOfAInB(VVI &V, VI &A, VI &B, bool useSet) {
    unordered_set<int> inBSet;
    VB inB;
    if( useSet ){
        inBSet.insert( ALL(B) );
    }else{
        inB = VB(V.size(),false);
        for(int d : B) inB[d] = true;
    }

    function< bool (int d) > inBFunction = [&inBSet, &inB, &useSet]( int d ) -> bool {
        if( useSet ) return (inBSet.count(d) > 0);
        else return inB[d];
    };

    VI res;
    if( useSet ){
        unordered_set<int> zb;
        for( int p : A ){
            for(int d : V[p]){
                if( inBFunction(d) ) zb.insert(d);
            }
        }
        res = VI(ALL(zb));
    }else{
        VB was(V.size(),false);
        for( int p : A ){
            for(int d : V[p]){
                if( inBFunction(d) && !was[d] ){
                    was[d] = true;
                    res.push_back(d);
                }
            }
        }
    }
    return res;
}


VVI GraphUtils::transposeGraph(VVI &v) {
    VVI g( SIZE(v) );
    REP( i,SIZE(v) )	REP( k,SIZE(v[i]) )	g[ v[i][k] ].PB(i);
    return g;
}

bool GraphUtils::isMaximalIndependentSet(VVI &V, VI &mis) {
    VB was(V.size(),false);
    for( int d : mis ){
        was[d] = true;
        for(int p : V[d]) was[p] = true;
    }

    for( int i=0; i<V.size(); i++ ){
        if( !was[i] ){
            return false;
        }
    }
    return true;
}

void GraphUtils::addEdge(VVI &V, int a, int b, bool directed) {
    V[a].push_back(b);
    if( !directed ) V[b].push_back(a);
}

void GraphUtils::removeEdge(VVI &V, int a, int b, bool directed) {
    auto rem = [ &V ](int a, int b){
        for( int i=(int)V[a].size()-1; i>=0; i-- ){
            if( V[a][i] == b ){
                swap(V[a][i], V[a].back() );
                V[a].pop_back();
                break;
            }
        }
    };


    rem(a,b);
    if(!directed) rem(b,a);

}

void GraphUtils::removeNeighborsOfNode(VVI &V, int a, unordered_set<int> neighbors) {
    for( int i=(int)V[a].size()-1; i>=0; i-- ){
        if( neighbors.count( V[a][i] ) > 0 ){
            swap(V[a][i], V[a].back() );
            V[a].pop_back();
        }
    }
}

void GraphUtils::removeNodeFromGraph(VVI &V, int a) {
    for( int d : V[a] ) removeEdge( V,d,a,true );
    V[a].clear();
}

int GraphUtils::countNodesWithDegree(VVI &V, int d, int D) {
    int res = 0;
    for( auto & v : V ){
        if(v.size() >= d && v.size() <= D) res++;
    }
    return res;
}

bool GraphUtils::containsEdge(VVI &V, int a, int b) {
    return find( ALL(V[a]),b ) != V[a].end();
}

int GraphUtils::countEdges(VVI &V) {
    int res = 0;
    for(auto& v : V) res += v.size();
    return res >> 1;
}

VVI GraphUtils::sortNodeNeighborhoods( VVI & V ){
    VVI V2(V.size());
    for(int i=0; i<V.size(); i++) V2[i].reserve(V[i].size());
    for( int i=0; i<V.size(); i++ ){
        for( int d : V[i] ) V2[d].push_back(i);
    }
    return V2;
}

VVI & GraphUtils::makeCompliment(VVI &V, VI nodes) {
    /*unordered_set<PII, pairhash> present;

    for( int p : nodes ){
        for(int d : V[p]){
            if( p < d ) present.insert({p,d});
        }

        unordered_set<int> zb(ALL(nodes));
        for( int j=(int)V[p].size()-1; j>=0; --j ){
            if( zb.count( V[p][j] ) ){
                swap( V[p][j], V[p].back() );
                V[p].pop_back();
            }
        }

    }

    for( int p : nodes ){
        for( int q : nodes ){
            if( p == q ) continue;
            if( present.count( { min(p,q), max(p,q) } ) > 0 ) continue;

            V[p].push_back(q);
        }
    }*/

   /* sort(ALL(nodes));
    for( int p : nodes ) sort( ALL(V[p]) );
    for( int p : nodes ){
        VI diff;
        set_difference( ALL(nodes), ALL(V[p]), back_inserter(diff) );
        diff.erase( find( ALL(diff),p ) );

        VI newNodes = diff;
        set_difference( ALL(V[p]), ALL(nodes), back_inserter(newNodes) );
        V[p] = newNodes;
    }*/

    return V;
}

VVI GraphUtils::getComplimentaryGraph(VVI &V) {
    VVI V2(V.size());
    for( int i=0; i<V.size(); i++ ){
        VI tmp = V[i];
        tmp.push_back(i);
        sort(ALL(tmp));
        V2[i] = getComplimentaryNodes( V,tmp );
    }

    bool correct = true;
    for( int i=0; i<V.size(); i++ ){
        VB was(V.size(),false);
        for( int d : V[i] ) was[d] = true;
        for(int d : V2[i]) was[d] = true;
        if( V[i].size() + V2[i].size() != V.size() -1 ) correct = false;
        for(int k=0; k<V.size(); k++) if( k != i && !was[k]) correct = false;
    }

    assert(correct);

    return V2;


//    VVI V2 = V;
//    return makeCompliment(V2, CombinatoricUtils::getRandomPermutation(V.size()) );
}

void GraphUtils::removeNodesFromGraph(VVI &V, VI nodes) {

    unordered_set<int> nd(ALL(nodes));
    VI neigh = getNeighborhood( V,nodes,true );
    for( int t : neigh ){
        if( nd.count(t) ) continue;

        for( int i = (int)V[t].size()-1; i>=0; i-- ){
            if( nd.count(V[t][i]) ){
                swap( V[t][i], V[t].back() );
                V[t].pop_back();
            }
        }
    }

    for( int d : nodes ) V[d].clear();
}

void GraphUtils::removeNodesFromGraph(VVI &V, VVI &W, VI nodes) {
    unordered_set<int> nd(ALL(nodes));
    VI neigh = getNeighborhood( V,nodes,true );
    for( int t : neigh ){
        if( nd.count(t) ) continue;

        for( int i = (int)V[t].size()-1; i>=0; i-- ){
            if( nd.count(V[t][i]) ){
                swap( V[t][i], V[t].back() );
                V[t].pop_back();

                swap( W[t][i], W[t].back() );
                W[t].pop_back();
            }
        }
    }

    for( int d : nodes ){
        V[d].clear();
        W[d].clear();
    }
}

void GraphUtils::writeGraphHumanReadable(VVI &V) {
    clog << "****" << endl;
    for( int i=0; i<V.size(); i++ ){
        clog << i << ": ";
        VI neigh = V[i];
        sort(ALL(neigh));
        for(int d : neigh) clog << d << "  ";
        clog << endl;
    }
}

void GraphUtils::mergeNodeToNode(VVI &V, int a, int b) {
//    assert( containsEdge( V,a,a ) == false );
//    assert( containsEdge( V,b,b ) == false );

    VI Na = V[a];
    Na.push_back(a);

    VI Nb = V[b];
    Nb.push_back(b);


    sort(ALL(Na));
    sort(ALL(Nb));

//    DEBUG(VI({a,b}));
//    DEBUG(Na);
//    DEBUG(Nb);

    VI diff;
    set_difference( ALL(Na), ALL(Nb), back_inserter(diff) ); // diff is Na \ Nb

//    DEBUG(diff);

    removeNodeFromGraph(V,a);
    for( int d : diff ){

//        assert( containsEdge( V,d,b ) == false );
        V[d].push_back(b);

//        assert( containsEdge( V,b,d ) == false );
        V[b].push_back(d);
    }

//    cerr << "GraphUtils::mergeNodeToNode not tested yet" << endl;
//    exit(1);
}

void GraphUtils::contractEdge(VVI &V, int a, int b) {
    mergeNodeToNode(V,a,b);
}

#ifndef QUICK_BUILD
void GraphUtils::writeBasicGraphStatistics(VVI &V) {
    clog << "V.size(): " << V.size() << endl;
    clog << "|E(V)|: " << countEdges(V) << endl;
    clog << "average degree: " << ( (double)2 * countEdges(V) ) / V.size() << endl;
    clog << "largest degree: " << max_element( ALL(V), []( VI& v, VI& w ){ return v.size() < w.size(); } )->size() << endl;
    int m = min( (int)V.size()-1,3 );

    VVI comps = ConnectedComponents::getConnectedComponents(V);
    clog << "There are " << comps.size() << " connected components" << endl;
    sort( ALL(comps), [](VI& v1, VI& v2){ return v1.size() > v2.size(); } );
    m = min( (int)comps.size(),7 );
    clog << "first " << m << " largest component sizes: "; for( int i=0; i<m; i++ ) clog << comps[i].size() << " "; clog << endl;
    clog << "first " << m << " smallest component sizes: "; for( int i=(int)comps.size()-m; i< comps.size(); i++ ) clog << comps[i].size() << " "; clog << endl;

    m = 3;
    VI degs(V.size(),0);
    for( VI& v : V ) degs[v.size()]++;
    for(int i=0; i<=m; i++) clog << "nodes with deg " << i << ": " << degs[i] << endl;
    clog << endl;

//    cerr << "isTree: " << Tree::isTree(V) << endl;

    auto artsAndBridges = BridgesAndArtPoints::getBridgesAndArtPoints(V);
    clog << "There are " << artsAndBridges.first.size() << " articulation points and " << artsAndBridges.second.size() << " bridges" << endl;
    clog << endl;



//    VI clq = CliqueExtension::findMaximalNodeCliqueExtension(V, true);
//    clog << "maximal clique found via sparse greedy extension: " << clq.size() << endl;

//    PartitionSVC vcCreator(V);
//    vcCreator.setSupressAllLogs(true); vcCreator.setMaxIterations( 200 ); vcCreator.setMaxRunTime(5'000);vcCreator.getSvcParams().alpha = 0.5;
//    vcCreator.getSvcParams().iterationsPerSubgraph = 200;vcCreator.getSvcParams().setInitialSolutionForSubgraph = false;vcCreator.setTakeFirstBestSolution(false);
//    vcCreator.getSvcParams().initialSolutionIterations = 1;vcCreator.run();
//    VI vc = ((SolutionVC*) vcCreator.getBestSolution())->getVC();
//    cerr << "Quickly found vertex cover size: " << vc.size() << "  and corresponding independent set size: " << V.size() - vc.size() << endl;

}
#endif // QUICKBUILD

void GraphUtils::removeEdges(VVI &V, VPII &edges, bool directed) {

//    cerr << "removeEdges not tested yet" << endl;

    if( directed == false ){
        int E = edges.size();
        for( int i=0; i<E; i++ ) edges.emplace_back( edges[i].second, edges[i].first ); // adding reverse edges to remove
    }

    sort( ALL(edges) );

    for( int i=0; i<edges.size(); i++ ){
        int p = i;
        unordered_set<int> toRemove;
        while( p < edges.size() && edges[p].first == edges[i].first ){
            toRemove.insert( edges[p].second );
            p++;
        }

        int t = edges[i].first;
        for( int k=(int)V[t].size()-1; k>=0; k-- ){
            int d = V[t][k];
            if( toRemove.count(d) ){
//                cerr << "Removing edge " << t << " -> " << V[t][k] << endl;
                swap( V[t][k], V[t].back() );
                V[t].pop_back();
            }
        }

        i = p-1;
    }

}

bool GraphUtils::isConnected(VVI &V) {
//    cerr << "isConnected not tested yet!" << endl; exit(1);
    VB was(V.size(),false);
    int cnt = 0;
    function< void(int) > dfs = [&V,&was, &dfs, &cnt](int num){
        was[num] = true;
        cnt++;
        for( int d : V[num] ) if( !was[d] ) dfs(d);
    };
    dfs(0);
    return (cnt == V.size());
}

VI GraphUtils::getNeighborhoodExclusive(VVI &V, VI &A, VB &helper) {
    VI res;
    for( int a : A ) helper[a] = true;
    for( int a : A ){
        for( int w : V[a] ){
            if(!helper[w]){
                res.push_back(w);
                helper[w] = true;
            }
        }
    }

    for( int a : A ) helper[a] = false;
    for(int d : res) helper[d] = false;
    return res;
}

