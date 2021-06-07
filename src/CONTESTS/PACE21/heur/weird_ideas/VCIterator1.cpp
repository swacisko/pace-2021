//
// Created by sylwester on 12/20/20.
//

#include <CONTESTS/Pace21/heur/weird_ideas/VCIterator1.h>
#include <graphs/GraphUtils.h>
#include <graphs/VertexCover/approximation/PartitionSVC.h>
#include <graphs/VertexCover/SolutionVC.h>


void VCIterator1::createHelperGraph(VVI & V) {
    helperGraph.clear();
    edgeToNodeMapper.clear();
    nodeToEdgeRemapper.clear();

    int N = V.size();

    VB neigh(N,false);
    int cnt = 0;


    for( int a=0; a<N; a++ ){

        for( int b : V[a] ) neigh[b] = true;

        for( int b : V[a] ){
            PII ab = { min(a,b), max(a,b) };

            for( int c : V[b] ){
                if( a >= c ) continue;
                if( neigh[c] ) continue; // if abc is a triangle, nothing to add

                PII ac = {a,c};
                PII bc = { min(b,c), max(b,c) };

                if( edgeToNodeMapper.count( ab ) == 0 ){
                    nodeToEdgeRemapper.push_back(ab);
                    edgeToNodeMapper[ab] = cnt++;
                    helperGraph.push_back(VI());
                }

                if( edgeToNodeMapper.count( ac ) == 0 ){
                    nodeToEdgeRemapper.push_back(ac);
                    edgeToNodeMapper[ac] = cnt++;
                    helperGraph.push_back(VI());
                }
                if( edgeToNodeMapper.count( bc ) == 0 ){
                    nodeToEdgeRemapper.push_back(bc);
                    edgeToNodeMapper[bc] = cnt++;
                    helperGraph.push_back(VI());
                }

                int abNode = edgeToNodeMapper[ab];
                int bcNode = edgeToNodeMapper[bc];
                int acNode = edgeToNodeMapper[ac];

                GraphUtils::addEdge( helperGraph, abNode, acNode );
                GraphUtils::addEdge( helperGraph, bcNode, acNode );
            }

        }

        for( int b : V[a] ) neigh[b] = false;

    }



}

VVI VCIterator1::nextIteration( VVI g ) {

    createHelperGraph(g);

    PartitionSVC partSVC(helperGraph);

    DEBUG(helperGraph);
    cerr << "edgeToNodeMapper:" << endl;
    for( auto t : edgeToNodeMapper ) cerr << t.first << " -> " << t.second << endl;
    GraphUtils::writeBasicGraphStatistics(helperGraph);
//    exit(1);

    partSVC.getSvcParams().alpha = 0.75;
    partSVC.getSvcParams().initialSolutionMaxRunTime = 5'000;
    partSVC.getSvcParams().maxRunTimeForSubgraph = 1'000;
    partSVC.getSvcParams().iterationsPerSubgraph = 100;
    partSVC.getSvcParams().setInitialSolutionForSubgraph = false;

    partSVC.getSvcParams().maxDeviationFromBestResult = 0;
    partSVC.getSvcParams().localOptimumViolationFrequency = 100;

    partSVC.setTakeFirstBestSolution(false);

    int SECONDS = 30;
    partSVC.setMaxRunTime(SECONDS * 1000);
    partSVC.getSvcParams().initialSolutionIterations = 1;
    partSVC.run();

    SolutionVC *s = (SolutionVC*) partSVC.getBestSolution();
    VI vc = s->getVC();

    VPII edg = GraphUtils::getGraphEdges(g);
    for(auto& e : edg) if( e.first > e.second ) swap( e.first, e.second );

    vector< unordered_set<int> > edges(g.size());
    for( PII e : edg ) edges[e.first].insert(e.second);

    VPII toRemove;
    for( auto v : vc ){
        auto edgeToEdit = nodeToEdgeRemapper[v];
        int a = edgeToEdit.first, b = edgeToEdit.second; // it already should be a < b
        bool addEdge = true;
        if( edges[a].count( a ) ) addEdge = false;

        if( addEdge ) GraphUtils::addEdge( g, a,b );
        else toRemove.push_back( {a,b} );
    }

    GraphUtils::removeEdges(g, toRemove);

    return g;
}


