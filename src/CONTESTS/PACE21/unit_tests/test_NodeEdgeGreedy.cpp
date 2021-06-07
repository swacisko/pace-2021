//
// Created by sylwester on 4/22/21.
//

#include <graphs/GraphUtils.h>
#include <CONTESTS/PACE21/heur/PaceUtils.h>
#include <CONTESTS/PACE21/heur/SwapCandidates/SwpCndNode.h>
#include <CONTESTS/PACE21/heur/StateImprovers/NodeEdgeGreedy.h>
#include "CONTESTS/PACE21/test_graphs.h"
#include "gtest/gtest.h"


class NEGreedyFixture : public ::testing::Test {
public:

    static void SetUpTestSuite() {
        old_clog_buf = clog.rdbuf();
        clog.rdbuf(cout.rdbuf());

        auto edges = CE_test_graphs::swpcndedge_test2_edges;
        V = CE_test_graphs::createVVIfromVPII(edges);
        VI partition = CE_test_graphs::swpcndedge_test2_partition;

        clg = new ClusterGraph(&V, partition);

        st = new State(*clg, SINGLE_NODES);
        VVI to_merge = {{0, 1, 2},
                        {3, 4},
                        {5, 6},
                        {7, 8}}; // the only single nodes remain 1 and 5
        st->mergeClusters(to_merge);
    }

    static void TearDownTestSuite() {
        delete clg; clg = nullptr;
        if(st != nullptr) delete st; st = nullptr;
        clog.rdbuf(old_clog_buf );
    }

    static ClusterGraph *clg;
    static streambuf* old_clog_buf;
    static State* st;
    static VVI V;
};

ClusterGraph* NEGreedyFixture::clg = nullptr;
streambuf* NEGreedyFixture::old_clog_buf = nullptr;
State* NEGreedyFixture::st = nullptr;
VVI NEGreedyFixture::V;


TEST_F( NEGreedyFixture, test_performance ){

    NodeEdgeGreedy neg(*st);
    neg.move_frequency = clg->N;

    neg.improve();

    ASSERT_EQ( neg.best_result, 33 );
}

TEST_F( NEGreedyFixture, test_performance2 ){
    for(int i=0; i<5; i++) {
        {
            auto edges = CE_test_graphs::cluster_graph_test2_edges;
            VVI V = CE_test_graphs::createVVIfromVPII(edges);
            VI partition = CE_test_graphs::cluster_graph_test2_partition;
            ClusterGraph clg(&V, partition);

            State st(clg, SINGLE_NODES);

            NodeEdgeGreedy neg(st);
            neg.move_frequency = clg.N;

            neg.improve();

            ASSERT_EQ(neg.best_result, 20);
        }

        {
            auto edges = CE_test_graphs::cluster_graph_test2_edges;
            VVI V = CE_test_graphs::createVVIfromVPII(edges);
            VI partition = CE_test_graphs::cluster_graph_test2_partition;
            iota(ALL(partition),0);
            ClusterGraph clg(&V, partition);

            State st(clg, SINGLE_NODES);

            NodeEdgeGreedy neg(st);
            neg.move_frequency = clg.N;

            neg.improve();

            ASSERT_LE(neg.best_result, 20);
            ASSERT_GE(neg.best_result, 18);
        }
    }
}

TEST_F( NEGreedyFixture, test_edge_swaps ){

    DEBUG(*st);

    for(int i=0; i<5; i++) {
        {
            auto edges = CE_test_graphs::cluster_graph_test2_edges;
            VVI V = CE_test_graphs::createVVIfromVPII(edges);
            VI partition = CE_test_graphs::cluster_graph_test2_partition;
            iota(ALL(partition),0);
            ClusterGraph clg(&V, partition);

            State st(clg, SINGLE_NODES);

            NodeEdgeGreedy neg(st);
            neg.move_frequency = 1;

            neg.improve();

            ASSERT_EQ(neg.best_result, 18);
        }

        {
            NodeEdgeGreedy neg(*st);
            neg.move_frequency = 2;
            neg.use_node_swaps = true;

            neg.improve();

            EXPECT_EQ(neg.best_result, 33);
        }
    }
}

TEST_F( NEGreedyFixture, test_interchange_value ){

    auto edges = CE_test_graphs::neg_interchange_edges;
    VVI V = CE_test_graphs::createVVIfromVPII(edges);
    VI partition = CE_test_graphs::neg_interchange_partition;
    ClusterGraph clg(&V, partition);

    State st(clg, SINGLE_NODES);
    VVI to_merge = { {0,1,2}, {3,4,5} };
    st.mergeClusters( to_merge );

    DEBUG(st);

    NodeEdgeGreedy neg(st);
    neg.move_frequency = 1;
    neg.use_node_interchanging = true;
    neg.use_edge_swaps = false;
    neg.use_node_swaps = true;

    neg.improve();


    SUCCEED() << "THIS TEST HAS NO ASSERTIONS!!" << endl;
}

TEST_F( NEGreedyFixture, test_chain2_swaps ){

    auto edges = CE_test_graphs::cluster_graph_test2_edges;
    VVI V = CE_test_graphs::createVVIfromVPII(edges);
    VI partition = CE_test_graphs::cluster_graph_test2_partition;
    ClusterGraph clg(&V, partition);

    State st(clg, SINGLE_NODES);
    VVI to_merge = { {0,4}, {1,2,3}, {5,7}, {6,8} };
    st.mergeClusters( to_merge );

    DEBUG(st);

    NodeEdgeGreedy neg(st);
    neg.move_frequency = 1;
    neg.use_node_interchanging = true;
    neg.use_edge_swaps = false;
    neg.use_node_swaps = false;
    neg.use_join_clusters = false;
    neg.use_node_interchanging = false;


    neg.improve();


    SUCCEED() << "THIS TEST HAS NO ASSERTIONS!!" << endl;
}

TEST_F( NEGreedyFixture, test_chain2_swaps2 ){

    auto edges = CE_test_graphs::swpcndedge_test2_edges;
    VVI V = CE_test_graphs::createVVIfromVPII(edges);
    VI partition = CE_test_graphs::swpcndedge_test2_partition;
    ClusterGraph clg(&V, partition);

    State st(clg, SINGLE_NODES);
    VVI to_merge = { {0,1,2}, {3,4}, {5,6}, {7,8} };
    st.mergeClusters( to_merge );

//    DEBUG(st);

    NodeEdgeGreedy neg(st);
    neg.move_frequency = 1;
    neg.use_node_interchanging = true;
    neg.use_edge_swaps = false;
    neg.use_node_swaps = false;
    neg.use_join_clusters = false;
    neg.use_node_interchanging = false;
//    neg.allow_perturbations = false;

    neg.improve();


    SUCCEED() << "THIS TEST HAS NO ASSERTIONS!!" << endl;
}