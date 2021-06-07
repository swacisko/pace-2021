//
// Created by sylwester on 4/1/21.
//

#include "CONTESTS/PACE21/heur/SwapCandidates/SwpCndEO.h"
#include "CONTESTS/PACE21/test_graphs.h"
#include "gtest/gtest.h"

class SwpCndEOFixture : public ::testing::Test {
public:

    static void SetUpTestSuite() {
        old_clog_buf = clog.rdbuf();
        clog.rdbuf(cout.rdbuf());

        auto edges = CE_test_graphs::swpcndeo_edges;
        V = CE_test_graphs::createVVIfromVPII(edges);
        VI partition = CE_test_graphs::swpcndeo_partition;

        clg = new ClusterGraph(&V, partition);

        st = new State(*clg, SINGLE_NODES);
        VVI to_merge = {{2,3,4,6,7}}; // the only single nodes remain 1 and 5
        st->mergeClusters(to_merge);
    }

    static void TearDownTestSuite() {
        delete clg;
        clog.rdbuf(old_clog_buf );
    }

    static ClusterGraph *clg;
    static streambuf* old_clog_buf;
    static State* st;
    static VVI V;
};

ClusterGraph* SwpCndEOFixture::clg = nullptr;
streambuf* SwpCndEOFixture::old_clog_buf = nullptr;
State* SwpCndEOFixture::st = nullptr;
VVI SwpCndEOFixture::V;

TEST_F(SwpCndEOFixture, test1){

//    st->createClNeighGraph();
//    DEBUG(st->cl_neigh_graph);
//    DEBUG(*st->clg);
//    DEBUG(*st);


    // clusters {0}, {1}, {5}, {2,3,4,6,7}, {}
    Cluster* cl = &st->clusters[3];

    VI ord = {2, 4, 7, 3, 6};
    for (int &d : ord) d = st->idInCl[d];
    ExpansionOrder eo(ord, cl);

    SwpCndEOCreator cr(*st);
    cr.keep_only_best_cluster_to_move_to = true;
    cr.keep_only_nonpositive_candidates = false;

    auto res = cr.createSwapCandidates(eo);

    DEBUG(res);
    ASSERT_EQ(res.size(), ord.size() - 1); // removed moving whole cluster to empty cluster
}


TEST_F(SwpCndEOFixture, test2){

    // clusters {0}, {1}, {5}, {2,3,4,6,7}, {}
    Cluster* cl = &st->clusters[3];

    VI ord = {2, 4, 7, 3, 6};
    for (int &d : ord) d = st->idInCl[d];
    ExpansionOrder eo(ord, cl);

    SwpCndEOCreator cr(*st);
    cr.keep_only_best_cluster_to_move_to = false;
    cr.keep_only_nonpositive_candidates = false;

    auto res = cr.createSwapCandidates(eo);

    ASSERT_EQ(res.size(), 18 - 1);// removed moving whole cluster to empty cluster

    for( auto & cnd : res ){
        if( cnd.getNodesToSwap() == VPII( { {2,4}, {4,4}, {7,4} } ) ){
            ASSERT_EQ(cnd.swpVal(),-8);
        }else if( cnd.getNodesToSwap() == VPII( { {2,0}, {4,0}, {7,0}, {3,0} } ) ){
            ASSERT_EQ(cnd.swpVal(),8);
        }
    }
}


TEST_F(SwpCndEOFixture, test3){

    // clusters {0}, {1}, {5}, {2,3,4,6,7}, {}
    Cluster* cl = &st->clusters[3];

    VI ord = {2, 4, 7, 3, 6};
    for (int &d : ord) d = st->idInCl[d];
    ExpansionOrder eo(ord, cl);

    SwpCndEOCreator cr(*st);
    cr.keep_only_best_cluster_to_move_to = false;
    cr.keep_only_nonpositive_candidates = true;

    auto res = cr.createSwapCandidates(eo);

    ASSERT_EQ(res.size(), 10 - 1); // removed moving whole cluster to empty cluster

    for( auto & cnd : res ){
        if( cnd.getNodesToSwap() == VPII( { {2,4}, {4,4}, {7,4} } ) ){
            ASSERT_EQ(cnd.swpVal(),-8);
        }else if( cnd.getNodesToSwap() == VPII( { {2,0}, {4,0}, {7,0}, {3,0} } ) ){
            ASSERT_EQ(cnd.swpVal(),8);
        }else if( cnd.getNodesToSwap() == VPII( { {2,1} } ) ){
            ASSERT_EQ(cnd.swpVal(),0);
        }
    }
}

TEST_F(SwpCndEOFixture, test4){

    // clusters {0}, {1}, {5}, {2,3,4,6,7}, {}
    Cluster* cl = &st->clusters[3];

    VI ord1 = {2, 4, 7, 3, 6};
    VI ord2 = {4,7,3,2,6};
    VI ord3 = { 2,3,4,6,7 };

    for (int &d : ord1) d = st->idInCl[d];
    for (int &d : ord2) d = st->idInCl[d];
    for (int &d : ord3) d = st->idInCl[d];
    ExpansionOrder eo1(ord1, cl);
    ExpansionOrder eo2(ord2, cl);
    ExpansionOrder eo3(ord3, cl);


    SwpCndEOCreator cr(*st);
    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        {
            cr.keep_only_best_cluster_to_move_to = true;
            cr.keep_only_nonpositive_candidates = false;

            vector<ExpansionOrder> orders = {eo1, eo2, eo3};
            vector<ExpansionOrder*> order_ptrs; for( auto & o : orders ) order_ptrs.push_back(&o);
            auto res = cr.createSwapCandidates(order_ptrs);
//        DEBUG(res);
            // {2}, {2,4}, {2,4,7}, {2,4,7,3}, {4}, {4,7}, {4,7,3}, {2,3}, {2,3,4}, {2,3,4,6}
            ASSERT_EQ(res.size(), 10);

            cr.keep_only_nonpositive_candidates = true;
            res = cr.createSwapCandidates(order_ptrs);

//        DEBUG(res);
            ASSERT_EQ(res.size(), 9 - 1); // removed moving whole cluster to empty cluster
        }
    }
}