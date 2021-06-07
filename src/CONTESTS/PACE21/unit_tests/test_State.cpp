//
// Created by sylwester on 3/22/21.
//

#include "CONTESTS/PACE21/heur/State.h"
#include "CONTESTS/PACE21/test_graphs.h"
#include "gtest/gtest.h"


class StateFixture : public ::testing::Test {
public:

    static void SetUpTestSuite() {
        old_clog_buf = clog.rdbuf();
        clog.rdbuf(cout.rdbuf());

        auto edges = CE_test_graphs::cluster_graph_test2_edges;
        VVI V = CE_test_graphs::createVVIfromVPII(edges);
        VI partition = CE_test_graphs::cluster_graph_test2_partition;

        clg = new ClusterGraph(&V, partition);
//        clog << "Cluster created" << endl;
    }

    static void TearDownTestSuite() {
        delete clg;
        clog.rdbuf(old_clog_buf );
    }

    static ClusterGraph *clg;
    static streambuf* old_clog_buf;
};

ClusterGraph* StateFixture::clg = nullptr;
streambuf* StateFixture::old_clog_buf = nullptr;

TEST_F( StateFixture, initialize_state_single_node ){
    ClusterGraph & clg = (*StateFixture::clg);

    ASSERT_EQ( clg.V[0].size(), 3 );
    ASSERT_EQ( clg.V[1].size(), 3 );
    ASSERT_EQ( clg.V[2].size(), 3 );
    ASSERT_EQ( clg.V[3].size(), 6 );
    ASSERT_EQ( clg.V[4].size(), 3 );
    ASSERT_EQ( clg.V[5].size(), 2 );
    ASSERT_EQ( clg.V[6].size(), 3 );
    ASSERT_EQ( clg.V[7].size(), 2 );
    ASSERT_EQ( clg.V[8].size(), 3 );

//    DEBUG(clg);

    State st(clg, SINGLE_NODES);

    auto equalSets = [](VI a, VI b){
        return equal(ALL(a), ALL(b));
    };

    for( int i=0; i<st.N; i++ ){
        ASSERT_PRED2( equalSets, st.clusters[i].g.nodes, VI(1,i) );
        ASSERT_EQ( st.idInCl[i],0);
        ASSERT_EQ( st.degInCl[i],0 );
        ASSERT_EQ( st.inCl[i],i );
    }
}

TEST_F( StateFixture, initialize_state_random_matching ){
    ClusterGraph & clg = (*StateFixture::clg);

    ASSERT_EQ( clg.V[0].size(), 3 );
    ASSERT_EQ( clg.V[1].size(), 3 );
    ASSERT_EQ( clg.V[2].size(), 3 );
    ASSERT_EQ( clg.V[3].size(), 6 );
    ASSERT_EQ( clg.V[4].size(), 3 );
    ASSERT_EQ( clg.V[5].size(), 2 );
    ASSERT_EQ( clg.V[6].size(), 3 );
    ASSERT_EQ( clg.V[7].size(), 2 );
    ASSERT_EQ( clg.V[8].size(), 3 );

    State st(clg, RANDOM_MATCHING);

    int cnt_cl_size_2 = 0;
    for(auto & cl : st.clusters ){
        ASSERT_LE(cl.size(), 2);
        if( cl.size() == 2 ){
            cnt_cl_size_2++;
            ASSERT_EQ( st.inCl[cl.g.nodes[0]], st.inCl[cl.g.nodes[1]] );
            int w = 0,a = cl.g.nodes[0], b = cl.g.nodes[1];
            for( auto e : clg.V[a] ) if( e.first == b ) w = e.second;
            ASSERT_LE( st.degInCl[a],w );
            ASSERT_LE( st.degInCl[b],w );
        }
    }

    ASSERT_GE( cnt_cl_size_2,3 );
}

TEST_F( StateFixture, mergeClusters ){
    ClusterGraph & clg = (*StateFixture::clg);
    State st(clg, SINGLE_NODES);

//    DEBUG(st.clusters);

    VVI to_merge = { {0,2,1}, {5,4,3}, {6,7} }; // sets to merge are always sorted, so we can pass them unsorted
    st.mergeClusters(to_merge);

    auto equalSets = [](VI a, VI b){
        return equal(ALL(a), ALL(b));
    };

    // cluster 8 should be kept untouched and therefore occurs as first after merging other clusters
    ASSERT_PRED2( equalSets, st.clusters[0].g.nodes, VI({8}) );

    ASSERT_PRED2( equalSets, st.clusters[1].g.nodes, VI({0,1,2}) );
    ASSERT_PRED2( equalSets, st.clusters[2].g.nodes, VI({3,4,5}) );
    ASSERT_PRED2( equalSets, st.clusters[3].g.nodes, VI({6,7}) );

    ASSERT_EQ(st.clusters.size(), 5); // including empty cluster

    ASSERT_EQ( st.inCl[0], st.inCl[1] );
    ASSERT_EQ( st.inCl[0], st.inCl[2] );
    ASSERT_EQ( st.inCl[0], 1 );

    ASSERT_EQ( st.inCl[3], st.inCl[4]);
    ASSERT_EQ( st.inCl[3], st.inCl[5]);
    ASSERT_EQ( st.inCl[3], 2);

    ASSERT_EQ( st.inCl[6], st.inCl[7]);
    ASSERT_EQ( st.inCl[6], 3);

    ASSERT_EQ( st.inCl[8], 0);
}

TEST_F( StateFixture, cl_neigh_graph ) {
    auto edges = CE_test_graphs::swpcndedge_test2_edges;
    VVI V = CE_test_graphs::createVVIfromVPII(edges);
    VI partition = CE_test_graphs::swpcndedge_test2_partition;
    ClusterGraph clg(&V, partition);
    State st(clg, SINGLE_NODES);

    VVI to_merge = {{0, 2, 1}, {4, 3}, {6, 5}, {7,8}}; // sets to merge are always sorted, so we can pass them unsorted
    st.mergeClusters(to_merge);
    st.createClNeighGraph(); // needs to be called after merging clusters.

//    DEBUG(st.clusters);
//    DEBUG(st.cl_neigh_graph);

    // clusters:  {0,1,2}, {3,4}, {5,6}, {7,8}

    ASSERT_EQ( st.cl_neigh_graph[0], VPII( {PII(1,5), PII(2,3) }) );
    ASSERT_EQ( st.cl_neigh_graph[1], VPII( {PII(3,9), PII(1,1) }) );
    ASSERT_EQ( st.cl_neigh_graph[2], VPII( {PII(3,11)} ) );
    ASSERT_EQ( st.cl_neigh_graph[3], VPII( {PII(0,3), PII(2,2)} ) );
    ASSERT_EQ( st.cl_neigh_graph[4], VPII( {PII(2,4), PII(0,3)} ) );
    ASSERT_EQ( st.cl_neigh_graph[5], VPII( {PII(1,6), PII(0,3)} ) );
    ASSERT_EQ( st.cl_neigh_graph[6], VPII() );
    ASSERT_EQ( st.cl_neigh_graph[7], VPII( {PII(0,9)} ) );
    ASSERT_EQ( st.cl_neigh_graph[8], VPII( {PII(0,11)} ) );
}

TEST_F( StateFixture, apply_swap ) {
//    DEBUG(*clg);

    {
        State st(*clg, SINGLE_NODES);
        VVI to_merge = { {0,1,2,3}, {4,5}, {6,8} }; // sets to merge are always sorted, so we can pass them unsorted
        st.mergeClusters(to_merge);

//        DEBUG(st);
        VPII to_swap = {{0, st.inCl[4]},
                        {2, st.getIdOfEmptyCluster()},
                        {3, st.getIdOfEmptyCluster() + 1},
                        {5, st.getIdOfEmptyCluster()},
                        {6, st.inCl[2]}, // after swap 6 should be with 1 in cluster, since 2 was removed to empty cluster
                        {7, st.inCl[8]}
        };

        st.applySwap(to_swap);

//        DEBUG(st);

        ASSERT_EQ(st.inCl[0], st.inCl[4]);
        ASSERT_EQ(st.inCl[1], st.inCl[6]);
        ASSERT_EQ(st.inCl[7], st.inCl[8]);
        ASSERT_EQ(st.inCl[2], st.inCl[5]);
        for (int x : VI({0, 1, 2, 4, 5, 6, 7, 8})) ASSERT_NE(st.inCl[3], st.inCl[x]);
    }

    { // almost the same as above, but {3} and {2,5} are moved to cluster with higher ids
        State st(*clg, SINGLE_NODES);
        VVI to_merge = { {0,1,2,3}, {4,5}, {6,8} }; // sets to merge are always sorted, so we can pass them unsorted
        st.mergeClusters(to_merge);

//        DEBUG(st);
        VPII to_swap = {{0, st.inCl[4]},
                        {2, st.getIdOfEmptyCluster()+4},
                        {3, st.getIdOfEmptyCluster()+2},
                        {5, st.getIdOfEmptyCluster()+4},
                        {6, st.inCl[2]}, // after swap 6 should be with 1 in cluster, since 2 was removed to empty cluster
                        {7, st.inCl[8]}
        };

        st.applySwap(to_swap);

//        DEBUG(st);

        ASSERT_EQ(st.inCl[0], st.inCl[4]);
        ASSERT_EQ(st.inCl[1], st.inCl[6]);
        ASSERT_EQ(st.inCl[7], st.inCl[8]);
        ASSERT_EQ(st.inCl[2], st.inCl[5]);
        for (int x : VI({0, 1, 2, 4, 5, 6, 7, 8})) ASSERT_NE(st.inCl[3], st.inCl[x]);
    }
}

TEST_F(StateFixture, calculate_state_result){
    State st(*clg, SINGLE_NODES);
    {
        VVI to_merge = {{0, 1, 2, 3}, {4, 5}, {6, 8}}; // sets to merge are always sorted, so we can pass them unsorted
        st.mergeClusters(to_merge);


        ASSERT_EQ(st.calculateResultForState(),
                  (8 + 3 + 5 + 8 + 3) // add edges in cluster {0,1,2,3}
                  + ( 0 ) // add edges in cluster {4,5}
                  + ( 0 ) // add edges in cluster {6,8}
                  + ( 0 ) // add edges in cluster {7}
                  + ( 2+2+2+3+1+1+1+2 ) // add edges between clusters
                  );
    }

    {
        st.initializeStateData(SINGLE_NODES); // creating new state
        VVI to_merge = {{0, 1, 2}, {4, 5,7}, {3,6, 8}}; // sets to merge are always sorted, so we can pass them unsorted
        st.mergeClusters(to_merge);


        ASSERT_EQ(st.calculateResultForState(),
                  (3+5) // add edges in cluster {0,1,2}
                  + ( 2+2 ) // add edges in cluster {4,5,7}
                  + ( 3+4 ) // add edges in cluster {3,6,8}
                  + ( 2+2+2+3+4+1+1+1+2 ) // add edges between clusters
        );
    }
}