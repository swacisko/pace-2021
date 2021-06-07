//
// Created by sylwester on 12/20/20.
//

#include <Makros.h>
#include <graphs/GraphUtils.h>
#include <utils/StandardUtils.h>
#include "CONTESTS/PACE21/heur/PaceUtils.h"

namespace PaceUtils{

   long long evaluateSolution(VVI & V, VI & partition ){
        long long res = 0;

        int N = V.size();
        VI clusterSizes(N,0);
        for( int i=0; i<N; i++ ){
            clusterSizes[ partition[i] ]++;
        }

        for( int i=0; i<N; i++ ){
            long long s = clusterSizes[i];
            res += s * (s-1) / 2;
        }

        auto edges = GraphUtils::getGraphEdges(V);
        for( PII e : edges ){
            auto a = e.first, b = e.second;
            if( partition[a] == partition[b] ) res--;
            else res++;
        }

        return res;
   }

    tuple<int, int, int> getEdgeModificationStatistics(VVI &V, VI &partition) {
        int insertions = 0, deletions = 0, modifications = 0;
        int N = V.size();

        VI clusterSizes(N,0);
        for( int i=0; i<N; i++ ){
            clusterSizes[ partition[i] ]++;
        }

        for( int i=0; i<N; i++ ){
            long long s = clusterSizes[i];
            modifications  += s * (s-1) / 2;
            insertions += s * (s-1) / 2;
        }

        auto edges = GraphUtils::getGraphEdges(V);
        for( PII e : edges ){
            auto a = e.first, b = e.second;
            if( partition[a] == partition[b] ){
                modifications--;
                insertions--;
            }
            else{
                modifications++;
                deletions++;
            }
        }

        return {insertions, deletions, modifications};
    }

    VVI partitionToClusters(VI &partition) {
        return StandardUtils::partitionToLayers(partition); // OLD
    }

    VI properlyRemapPartition(VI part){
        int M = *max_element(ALL(part));
        VI mapper(M+1,-1);
        int cnt = 0;
        for( int i=0; i<part.size(); i++ ){
            int p = part[i];
            if(mapper[p] == -1){
                mapper[p] = cnt;
                cnt++;
            }
            part[i] = mapper[p];
        }
        assert(cnt <= part.size());
        return part;
    }

    VI clustersToPartition(VVI &clusters) {
        return StandardUtils::layersToPartition(clusters);
    }


    vector<tuple<int, int, int> > getInducedP2Paths(VVI &V) {
        vector<tuple<int, int, int>> paths;

        int N = V.size();
        VB neigh(N,false);

        for( int a=0; a<N; a++ ){

            for( int b : V[a] ) neigh[b] = true;

            for( int b : V[a] ){

                for( int c : V[b] ){
                    if( a >= c ) continue;
                    if( neigh[c] ) continue; // if abc is a triangle, nothing to add

                    paths.emplace_back( a,b,c );
                }
            }

            for( int b : V[a] ) neigh[b] = false;
        }
        return paths;
    }

    ClusterGraph convertToClusterGraph(VVI &V, VI &inCl) {
        return ClusterGraph( &V, inCl );
    }

    void test() {

        {
            VVI V = {
                    {1, 2},
                    {0, 4},
                    {0, 3, 4},
                    {2},
                    {1, 2, 5},
                    {4, 6},
                    {5}
            };

            int N = V.size();
            VI part(N, 0);
            part[4] = part[5] = part[6] = 1;

            LL res = evaluateSolution(V, part);

            assert(res == 6);

            clog << "Evaluate solution unweighed passed" << endl << endl;
        }

        {
            VVI clusters = {
                    {1,3,5},
                    {2,4},
                    {7}
                    // 0,6 is missing on purpose!
            };

            VI part = { -1,0, 1,0,1,0,-1,2 };

            VI part2 = clustersToPartition(clusters);
            for( int i=0; i<8; i++ ) assert( part[i] == part2[i] );

            VVI layers = partitionToClusters(part);
            for(int i=0; i<3; i++){
                assert( equal( ALL(layers[i]), ALL(clusters[i])  ) );
            }

            clog << "partitionToClusters and clustersToPartition passed" << endl;
        }


        clog << "Pace Utils tests passed" << endl;
    }

    void convertDIMACSTestToVVI() {
        int N,E;
        string s;
        cin >> s >> s >> N >> E;
        VVI V(N);
        for( int i=0; i<E; i++ ){
            int a,b;
            cin >> a >> b;
            V[a].push_back(b);
            V[b].push_back(a);
        }
        for(int i=0; i<N; i++){
            if( i >= 1 ) cout << ",";
            if( i % 7 == 0 ) cout << endl;
            cout << "{";
            for(int j=0; j<V[i].size(); j++){
                if( j >= 1 ) cout << ",";
                cout << V[i][j];
            }
            cout << "}";
        }
        exit(1);
    }

    VI mapClgPartitionToOriginalPartition(ClusterGraph &clg, VI &part) {
        VI res(clg.origV->size(),-1);
        for( int i=0; i<part.size(); i++ ){
            for( int d : clg.clusterNodes[i] ) res[d] = part[i];
        }
        return properlyRemapPartition(res);
    }

    VI mapOriginalPartitionToClgPartition(ClusterGraph &clg, VI &part) {
        VI res(clg.N,-1);
        for( int i=0; i<clg.N; i++ ) res[i] = part[ clg.clusterNodes[i][0] ];
        return properlyRemapPartition(res);
   }

    long long evaluateState(State &st) {
       VI part = mapClgPartitionToOriginalPartition(*st.clg, st.inCl);
        return evaluateSolution( *st.clg->origV, part );
    }



}