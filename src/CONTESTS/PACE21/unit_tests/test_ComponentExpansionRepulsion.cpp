//
// Created by sylwester on 4/8/21.
//

#include <graphs/GraphUtils.h>
#include "CONTESTS/PACE21/heur/SwapCandidates/ComponentExpansionRepulsion.h"
#include "CONTESTS/PACE21/test_graphs.h"
#include "gtest/gtest.h"


class ComponentExpansionRepulsionFixture : public ::testing::Test {
public:
    static void SetUpTestSuite() {
        old_clog_buf = clog.rdbuf();
        clog.rdbuf(cout.rdbuf());

        auto edges = CE_test_graphs::comp_exp_rep_edges;
        VVI V = CE_test_graphs::createVVIfromVPII(edges);
        VI partition = CE_test_graphs::comp_exp_rep_partition;

        clg = new ClusterGraph(&V, partition);

        st = new State(*clg, SINGLE_NODES);
        VVI to_merge = {{2, 0, 1},
                        {3, 4},
                        {5, 6}}; // the only single nodes remain 1 and 5
        st->mergeClusters(to_merge);
    }

    static void TearDownTestSuite() {
        delete clg; clg = nullptr;
        delete st; st = nullptr;
        clog.rdbuf(old_clog_buf );
    }

    static ClusterGraph *clg;
    static streambuf* old_clog_buf;
    static State* st;
};

ClusterGraph* ComponentExpansionRepulsionFixture::clg = nullptr;
streambuf* ComponentExpansionRepulsionFixture::old_clog_buf = nullptr;
State* ComponentExpansionRepulsionFixture::st = nullptr;

TEST_F(ComponentExpansionRepulsionFixture, test1){
//    DEBUG( *clg );
//    DEBUG(*st);

    st->createClNeighGraph();

    ComponentExpansionRepulsion cer(*st);

    for(int i=0; i<5; i++) {
        {
            cer.keep_only_nonpositive_candidates = true;
//            DEBUG(st->cl_neigh_graph);

            Cluster* cl = &st->clusters[1];
            auto res = cer.createSwapCandidates(*cl);
            ASSERT_EQ(res.second.size(), 1);
            ASSERT_EQ(res.second[0].getNodesToSwap(), VPII({{1, 3}}));
            ASSERT_EQ(res.second[0].swpVal(), 0);

            auto induced_orders = ExpansionOrderRepulsion::induceOrders( *st, *res.first );
            for( auto & eo : induced_orders ) ASSERT_EQ( eo.cl->id, cl->id );

            ASSERT_EQ( induced_orders.size(),3 );

            ASSERT_EQ( induced_orders[0].ord, VI({0}) );

            ASSERT_EQ( induced_orders[1].ord, VI({1}) );

            ASSERT_EQ( induced_orders[2].ord, VI({2}) );

            delete res.first;

            int N = cer.N;
            ASSERT_EQ( cer.was, VB(2*N,false) );
            ASSERT_EQ( cer.inS, VB(N,false) );
            ASSERT_EQ( cer.cluster_weights, VI(cer.cluster_weights.size(),0) );
            ASSERT_EQ( cer.eInS, VI(cer.eInS.size(),0) );

            ASSERT_TRUE(cer.S.empty());
            ASSERT_TRUE(cer.cl_neigh.empty());
            ASSERT_EQ(cer.sumNWinS ,0);

            for( int i=0; i<N; i++ ){
                ASSERT_TRUE(cer.edges_to_cluster[i].empty() );
                ASSERT_TRUE(cer.swpval_core_for_node[i].empty() );
                ASSERT_TRUE(cer.cluster_sets[i].empty() );
            }
        }

        {
            cer.keep_only_nonpositive_candidates = false; // the only difference from test above
//            DEBUG(st->cl_neigh_graph);

            auto res = cer.createSwapCandidates(st->clusters[1]);
            ASSERT_EQ(res.second.size(), 3);
            ASSERT_EQ(res.second[0].getNodesToSwap(), VPII({{1, 3}}));
            ASSERT_EQ(res.second[0].swpVal(), 0);

            ASSERT_EQ(res.second[1].getNodesToSwap(), VPII({{1, 3}, {2, 2}}));
            ASSERT_EQ(res.second[1].swpVal(), 1);

            ASSERT_EQ(res.second[2].getNodesToSwap(), VPII({{1, 3}, {2, 2}, {0, 5}}));
            ASSERT_EQ(res.second[2].swpVal(), 1);

            delete res.first;
        }

        { // other test, starting from different cluster
            cer.keep_only_nonpositive_candidates = false;

            auto res = cer.createSwapCandidates(st->clusters[st->inCl[3]]);
            ASSERT_EQ(res.second.size(), 2);

            ASSERT_EQ(res.second[0].swpVal(), -2);

            ASSERT_EQ(res.second[1].getNodesToSwap(),
                      VPII({{4, st->getIdOfEmptyCluster()+1}, {3, st->getIdOfEmptyCluster()+2}}));
            ASSERT_EQ(res.second[1].swpVal(), -2);

            delete res.first;
        }

        { // other test, starting from different cluster
            cer.keep_only_nonpositive_candidates = false;
            auto res = cer.createSwapCandidates(st->clusters[st->inCl[6]]);
            ASSERT_EQ(res.second.size(), 2);
            delete res.first;


            cer.keep_only_nonpositive_candidates = false;
            res = cer.createSwapCandidates(st->clusters[st->inCl[7]]);
            ASSERT_EQ(res.second.size(), 0);
            // there is created one swap candidate {7,empty_luster}. Since
            // all nodes from cluster are moved to the same empty cluster, we did not return it as a candidate
            delete res.first;
        }
    }
}

TEST_F(ComponentExpansionRepulsionFixture, test2){
    auto edges = CE_test_graphs::comp_exp_rep_edges;
    VVI V = CE_test_graphs::createVVIfromVPII(edges);

    VPII to_remove = { {2,9}, {1,9}, {3,7} };
    GraphUtils::removeEdges( V, to_remove );

    VI partition = CE_test_graphs::comp_exp_rep_partition;
    ClusterGraph clg(&V, partition);
    State st(clg, SINGLE_NODES);
    VVI to_merge = {{0, 1, 2}, {3, 4}, {5, 6}}; // the only single nodes remain 1 and 5
    st.mergeClusters(to_merge);
    st.createClNeighGraph();

    ComponentExpansionRepulsion cer(st);
    for(int i=0; i<5; i++) {
        cer.keep_only_nonpositive_candidates = false;

        auto res = cer.createSwapCandidates(st.clusters[1]);
        ASSERT_EQ(res.second.size(), 3);
        ASSERT_EQ(res.second[0].getNodesToSwap(), VPII({{2, 5}}));
        ASSERT_EQ(res.second[0].swpVal(), 2);
        ASSERT_EQ(res.second[1].getNodesToSwap(), VPII({{2, 5}, {0, 5}}));
        ASSERT_EQ(res.second[1].swpVal(), 2);
        ASSERT_EQ(res.second[2].getNodesToSwap(), VPII({{2, 5}, {0, 5}, {1, 5}}));
        ASSERT_EQ(res.second[2].swpVal(), 0);

        auto induced_orders = ExpansionOrderRepulsion::induceOrders( st, *res.first );

        ASSERT_EQ( induced_orders.size(),1 );

        ASSERT_EQ( induced_orders[0].ord, VI({2,0,1}) );
        ASSERT_EQ( induced_orders[0].cl->id, 1 );

        delete res.first;
    }
}


TEST_F(ComponentExpansionRepulsionFixture, test3){
    auto edges = CE_test_graphs::comp_exp_rep_edges;
    VVI V = CE_test_graphs::createVVIfromVPII(edges);

    VPII to_remove = { {2,9}, {1,9}, {3,7},
                       {0,2}, {0,3}, {1,4}};
    GraphUtils::removeEdges( V, to_remove );

    VI partition = CE_test_graphs::comp_exp_rep_partition;
    ClusterGraph clg(&V, partition);
    State st(clg, SINGLE_NODES);
    VVI to_merge = {{0, 1, 2}, {3, 4}, {5, 6}}; // the only single nodes remain 1 and 5
    st.mergeClusters(to_merge);
    st.createClNeighGraph();


    ComponentExpansionRepulsion cer(st);
    for(int i=0; i<5; i++) {
        cer.keep_only_nonpositive_candidates = false;

        auto res = cer.createSwapCandidates(st.clusters[1]);
        ASSERT_EQ(res.second.size(), 3);
        ASSERT_EQ(res.second[0].getNodesToSwap(), VPII({{1, 5}}));
        ASSERT_EQ(res.second[0].swpVal(), -2);
        ASSERT_EQ(res.second[1].getNodesToSwap(), VPII({{1, 5}, {2, 6}}));
        ASSERT_EQ(res.second[1].swpVal(), -2);
        ASSERT_EQ(res.second[2].getNodesToSwap(), VPII({{1, 5}, {2, 6}, {0, 7}}));
        ASSERT_EQ(res.second[2].swpVal(), -2);

        auto induced_orders = ExpansionOrderRepulsion::induceOrders( st, *res.first );
//            clog << "induced_orders of order " << StandardUtils::zip( res.first->ord, res.first->move_to ) << endl;
//            for( auto & eo : induced_orders ){
//                clog << eo.ord << "  -->  cluster: " << *eo.cl << endl;
//            }

        ASSERT_EQ( induced_orders.size(),3 );

        ASSERT_EQ( induced_orders[0].ord, VI({st.idInCl[0]}) );

        ASSERT_EQ( induced_orders[1].ord, VI({st.idInCl[1]}) );

        ASSERT_EQ( induced_orders[2].ord, VI({st.idInCl[2]}) );

        delete res.first;
    }
}