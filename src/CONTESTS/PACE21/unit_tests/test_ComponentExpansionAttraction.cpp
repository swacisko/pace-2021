//
// Created by sylwester on 4/12/21.
//


#include <graphs/GraphUtils.h>
#include "CONTESTS/PACE21/heur/SwapCandidates/ComponentExpansionAttraction.h"
#include "CONTESTS/PACE21/test_graphs.h"
#include "gtest/gtest.h"


class ComponentExpansionAttractionFixture : public ::testing::Test {
public:
    static void SetUpTestSuite() {
        old_clog_buf = clog.rdbuf();
        clog.rdbuf(cout.rdbuf());

        auto edges = CE_test_graphs::swpcndedge_test2_edges;
        V = CE_test_graphs::createVVIfromVPII(edges);
        VI partition = CE_test_graphs::swpcndedge_test2_partition;

        clg = new ClusterGraph(&V, partition);

        st = new State(*clg, SINGLE_NODES);
        VVI to_merge = {{0,1,2}, {3,4}, {5,6}, {7,8}}; // the only single nodes remain 1 and 5
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
    static VVI V;
};

ClusterGraph* ComponentExpansionAttractionFixture::clg = nullptr;
streambuf* ComponentExpansionAttractionFixture::old_clog_buf = nullptr;
State* ComponentExpansionAttractionFixture::st = nullptr;
VVI ComponentExpansionAttractionFixture::V;


TEST_F(ComponentExpansionAttractionFixture, test1){

    ComponentExpansionAttraction cea(*st);

    for( int i=0; i<5; i++ ) {
        Cluster *cl = &st->clusters[0];

        cea.keep_only_nonpositive_candidates = false;
        auto res = cea.createSwapCandidates(*cl);

        ASSERT_EQ(res.second.size(), 6);

        ASSERT_EQ(res.second[0].getNodesToSwap(), VPII({{8, 0}}));
        ASSERT_EQ(res.second[0].swpVal(), -7);

        ASSERT_EQ(res.second[1].getNodesToSwap(), VPII({{8, 0}, {3, 0}}));
        ASSERT_EQ(res.second[1].swpVal(), -2);

        ASSERT_EQ(res.second[2].getNodesToSwap(), VPII({{8, 0}, {3, 0}, {5, 0}}));
        ASSERT_EQ(res.second[2].swpVal(), 12);

        ASSERT_EQ(res.second[3].getNodesToSwap(), VPII({{8, 0}, {3, 0}, {5, 0}, {4, 0}}));
        ASSERT_EQ(res.second[3].swpVal(), 24);

        ASSERT_EQ(res.second[4].getNodesToSwap(), VPII({{8, 0}, {3, 0}, {5, 0}, {4, 0}, {6, 0}}));
        ASSERT_EQ(res.second[4].swpVal(), 52);

        ASSERT_EQ(res.second[5].getNodesToSwap(), VPII({{8, 0}, {3, 0}, {5, 0}, {4, 0}, {6, 0}, {7, 0}}));
        ASSERT_EQ(res.second[5].swpVal(), 88);

        delete res.first;


        // another test
        cea.keep_only_nonpositive_candidates = true;
        res = cea.createSwapCandidates(*cl);
        ASSERT_EQ( res.second.size(), 2 );

        delete res.first;
    }
}

TEST_F(ComponentExpansionAttractionFixture, test2){
    auto edges = CE_test_graphs::cluster_graph_test2_edges;
    VVI V = CE_test_graphs::createVVIfromVPII(edges);
    VI partition = CE_test_graphs::cluster_graph_test2_partition;
    ClusterGraph clg(&V, partition);

    State st(clg, SINGLE_NODES);
    VVI to_merge = { {0,4}, {5,7}, {1,2,3} };
    st.mergeClusters(to_merge);

    ComponentExpansionAttraction cea(st);

    for( int i=0; i<5; i++ ) {
        DEBUG(i);
        Cluster *cl = &st.clusters[ st.inCl[5] ];

        cea.keep_only_nonpositive_candidates = false;
        auto res = cea.createSwapCandidates(*cl);

        ASSERT_EQ(res.second.size(), 6);

        ASSERT_EQ(res.second[0].getNodesToSwap(), VPII({{3, cl->id}}));
        ASSERT_EQ(res.second[0].swpVal(), -8);

        ASSERT_EQ(res.second[1].getNodesToSwap(), VPII({{3, cl->id}, {4, cl->id}}));
        ASSERT_EQ(res.second[1].swpVal(), -6);

        ASSERT_EQ(res.second[2].getNodesToSwap(), VPII({{3, cl->id}, {4, cl->id}, {8,cl->id}}));
        ASSERT_EQ(res.second[2].swpVal(), -2);

        ASSERT_EQ(res.second[3].getNodesToSwap(), VPII({{3, cl->id}, {4, cl->id}, {8,cl->id}, {2,cl->id}}));

        ASSERT_EQ(res.second[4].getNodesToSwap(),
                  VPII({{3, cl->id}, {4, cl->id}, {8,cl->id}, {2,cl->id}, {1,cl->id}}));

        ASSERT_EQ(res.second[5].getNodesToSwap(),
                  VPII({{3, cl->id}, {4, cl->id}, {8,cl->id}, {2,cl->id}, {1,cl->id}, {0,cl->id}}));

        auto induced_orders = ExpansionOrderAttraction::induceOrders( st, *res.first );
//        DEBUG(induced_orders);
        ASSERT_EQ(induced_orders.size(), 3 );
        ASSERT_EQ( induced_orders[0].ord, VI({st.idInCl[8]}) );
        ASSERT_EQ( induced_orders[1].ord, VI({st.idInCl[3],st.idInCl[2],st.idInCl[1]}) );
        ASSERT_EQ( induced_orders[2].ord, VI({st.idInCl[4],st.idInCl[0]}) );

        delete res.first;


        // another test
        cea.keep_only_nonpositive_candidates = true;
        res = cea.createSwapCandidates(*cl);
        ASSERT_EQ(res.second.size(), 3);

        delete res.first;

        { // another test, starting from different cluster
            cl = &st.clusters[ st.inCl[6] ];
            cea.keep_only_nonpositive_candidates = true;
            res = cea.createSwapCandidates(*cl);
//            DEBUG(res.second);
            ASSERT_EQ(res.second.size(), 3);
            ASSERT_EQ( res.second[0].swpVal(),-8 );
            ASSERT_EQ( res.second[1].swpVal(),-5 );
            ASSERT_EQ( res.second[2].swpVal(),-1 );

            ASSERT_EQ(res.second[2].getNodesToSwap(), VPII({{3, cl->id}, {2, cl->id}, {8,cl->id}}));

            delete res.first;

        }

        for( int i=0; i<st.N; i++ ){ // checking if everything was cleared
            ASSERT_EQ( cea.was[i], false );
            ASSERT_EQ( cea.inS[i], false );
            ASSERT_EQ( cea.eToS[i], 0 );
            ASSERT_EQ( cea.edges_in_cl[i], 0 );
            if( i < cea.cluster_weights.size() ) ASSERT_EQ( cea.cluster_weights[i], 0 );
            if( i < cea.cluster_cores.size() ) ASSERT_TRUE( cea.cluster_cores[i].empty() );
        }

    }
}



