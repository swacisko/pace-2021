//
// Created by sylwester on 3/16/21.
//

#include <graphs/GraphUtils.h>
#include <CONTESTS/PACE21/test_graphs.h>
#include "CONTESTS/PACE21/heur/ClusterGraph.h"

ostream& operator<<( ostream& str, ClusterGraph& g ){
    for( int i=0; i<g.N; i++ ){
        str << "Node #" << i << "-->   nw: " << g.nw(i) << ",  V[i]: " << g.V[i] << endl;
    }

    return str;
}



ClusterGraph::ClusterGraph(VVI *V, VI partition) {
    origV = V;
    this->partition = partition;

    createGraph();
}

void ClusterGraph::createGraph() {
    int N = origV->size();

//    partition = origPartition;

    VI part_id_mapper(N, -1); // part_id_mapper[i] is the mapped id of the partition (cluster) i

    int cnt = 0;
    for( int i=0; i<partition.size(); i++ ){
        int p = partition[i];
        if( p < 0 ) continue;

        if(part_id_mapper[p] == -1 ){
            part_id_mapper[p] = cnt;
            cnt++;
            clusterNodes.push_back(VI());
        }

        partition[i] = part_id_mapper[p];
        clusterNodes[ part_id_mapper[p] ].push_back(i);
    }

//    DEBUG(clusterNodes);

    node_weights = VI(cnt,0);

    map<PII,int> edge_weights;
    for( int i=0; i<N; i++ ){
        if( partition[i] < 0 ) continue;

        node_weights[ partition[i] ]++;

        for( int d : (*origV)[i] ){
            if( i < d ){
                int a = partition[i];
                int b = partition[d];
                if(b<0) continue;

                if(a>b) swap(a,b);
                if( a == b ) continue;
                edge_weights[ {a,b} ]++;
            }
        }
    }

    this->N = cnt;
    V = VVPII(cnt);
    for( auto [e,w] : edge_weights ){
        auto [a,b] = e;
        GraphUtils::addEdge( V,a,b,w );
    }
}


VPII ClusterGraph::getClusterSizesCount() {
    unordered_map<int,int,fib_hash> zb;
    for( int d : node_weights ) zb[d]++;
    VPII res(ALL(zb));
    sort( res.rbegin(), res.rend() );
    return res;
}


void ClusterGraph::test(){

    VVI V = CE_test_graphs::cluster_graph_test;
    VI part = CE_test_graphs::cluster_graph_test_partition;

    ClusterGraph g(&V,part);

    VI part_id_mapper_test = { 0,1,1,2,3,3,4,4,4,5, -1, 6 };
    DEBUG(g.partition);

    assert( equal( ALL(g.partition), ALL(part_id_mapper_test) ) );

    DEBUG(g);
    DEBUG(g.clusterNodes);
}


InducedClusterGraph induce(ClusterGraph &clg, VI &nodes) {
    InducedGraphPI ig = GraphInducer::induce( clg.V, nodes );
    InducedClusterGraph g(clg, ig);
    return g;
}

InducedClusterGraph induce(ClusterGraph &clg, VI &nodes, VI &helper) {
    InducedGraphPI ig = GraphInducer::induceNoPerm( clg.V, nodes, helper );
    InducedClusterGraph g(clg, ig);
    return g;
}


InducedClusterGraph::InducedClusterGraph(ClusterGraph & clg, InducedGraphPI ig) {
    swap(par, ig.par);
    swap(nodes, ig.nodes);
    swap(perm, ig.perm);
    swap(edges, ig.edges);
    swap(V, ig.V);
    clPar = &clg;
    node_weights.resize(nodes.size());
    for( int i=0; i<nodes.size(); i++ ) node_weights[i] = clg.nw( nodes[i] );
}

ostream& operator<<( ostream& str, InducedClusterGraph& cl ){
    str
//        << "Par: " << *cl.par << endl
        << "Nodes: " << cl.nodes << endl
        << "Node weights: " << cl.node_weights << endl
        << "Perm: " << cl.perm << endl
        << "V: " << cl.V << endl;
    return str;

}

