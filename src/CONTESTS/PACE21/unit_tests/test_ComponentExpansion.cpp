//
// Created by sylwester on 3/21/21.
//

#include "CONTESTS/PACE21/test_graphs.h"
#include "CONTESTS/PACE21/heur/EOCreators/ComponentExpansion.h"
#include "graphs/GraphUtils.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

using namespace ::testing;
using CE=ComponentExpansion;


class EOFixture : public ::testing::Test {
public:

    static void SetUpTestSuite() {
        old_clog_buf = clog.rdbuf();
        clog.rdbuf(cout.rdbuf());

        VI nw = CE_test_graphs::EOnw;
        vector<tuple<int, int, int>> edges = CE_test_graphs::EOgraphEdges;
//        DEBUG(nw);
        int N = 0;
        for (auto &[a, b, c] : edges) N = max(N, max(a, b));
        auto g = VVPII(N + 1);
        for (auto[a, b, c] : edges) {
            g[a].emplace_back(b, c);
            g[b].emplace_back(a, c);
        }

//        clog << "g created: " << g << endl;
        clg = new ClusterGraph();
        clg->node_weights = nw;
        clg->V = g;
        clg->N = clg->V.size();
//        clog << "Cluster graph created" << endl;

        VI nodes(clg->N);
        iota(ALL(nodes), 0);
        cl = new Cluster(*clg, nodes, 0);
//        clog << "Cluster created" << endl;
    }

    static void TearDownTestSuite() {
        delete clg;
        delete cl;

        clog.rdbuf(old_clog_buf );
    }

    static ClusterGraph *clg;
    static Cluster *cl;
    static streambuf* old_clog_buf;
};

ClusterGraph* EOFixture::clg = nullptr;
Cluster* EOFixture::cl = nullptr;
streambuf* EOFixture::old_clog_buf = nullptr;



TEST_F(EOFixture, rule_12 ){

    CE ce(*cl);

    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        for (int use_heap = 0; use_heap <= 1; use_heap++) {
            ce.setCmpRules({1, 2});
            auto res = ce.getExpansionOrder({0, 1, 2}, use_heap);

            VI expected = {0, 1, 2, 3, 6, 12, 4, 10, 5, 11, 8, 9, 7};
            EXPECT_TRUE(equal(ALL(res.ord), ALL(expected)));

            { // another test
                ce.setCmpRules({3});
                auto res = ce.getExpansionOrder({0, 1, 2}, use_heap);

                VI expected = {0, 1, 2, 3, 8, 7, 6, 12, 4, 10, 9, 5, 11};
                EXPECT_TRUE(equal(ALL(res.ord), ALL(expected)));
            }
        }
    }
}

TEST_F(EOFixture, rule_21 ){

    CE ce(*cl);

    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        for (int use_heap = 0; use_heap <= 1; use_heap++) {
            ce.setCmpRules({2, 1});
            auto res = ce.getExpansionOrder({0, 1, 2}, use_heap);

            VI expected = {0, 1, 2, 3, 8, 7, 4, 10, 9, 5, 11, 6, 12};
            EXPECT_TRUE(equal(ALL(res.ord), ALL(expected)));

//            ENDL(5);
        }
    }
}

TEST_F(EOFixture, rule_3 ){

    CE ce(*cl);

    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        for (int use_heap = 0; use_heap <= 1; use_heap++) {
            ce.setCmpRules({3});
            auto res = ce.getExpansionOrder({0, 1, 2}, use_heap);

            VI expected = {0, 1, 2, 3, 8, 7, 6, 12, 4, 10, 9, 5, 11};
            EXPECT_TRUE(equal(ALL(res.ord), ALL(expected)));

//            ENDL(5);
        }
    }
}

TEST_F(EOFixture, rule_4 ){
    CE ce(*cl);

    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        ce.setCmpRules({4});
        auto res = ce.getExpansionOrder({0, 1, 2}, false);

        VI expected = {0, 1, 2, 3, 4, 10, 6, 12, 5, 11, 8, 9, 7};
        EXPECT_TRUE(equal(ALL(res.ord), ALL(expected)));

        { // another test
            ce.setCmpRules({1, 2});
            auto res = ce.getExpansionOrder({0, 1, 2}, false);

            VI expected = {0, 1, 2, 3, 6, 12, 4, 10, 5, 11, 8, 9, 7};
            EXPECT_TRUE(equal(ALL(res.ord), ALL(expected)));
        }
    }

//    ENDL(5);
}

TEST_F(EOFixture, rule_5 ){
    CE ce(*cl);
    ce.setCmpRules({5});
    auto res = ce.getExpansionOrder({0, 1, 2}, false);

    VI expected = {0, 1, 2, 6, 12, 5, 11,4, 9, 10, 3, 7, 8 };
    EXPECT_TRUE(equal(ALL(res.ord), ALL(expected)));

//    ENDL(5);
}


TEST_F(EOFixture, rule_10 ){
    CE ce(*cl);

    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        ce.setCmpRules({10});
        auto res = ce.getExpansionOrder({0, 1, 2}, false);

        VI expected = {0, 1, 2, 3, 8, 7, 6, 12, 4, 10, 9, 5, 11};
        EXPECT_TRUE(equal(ALL(res.ord), ALL(expected)));
    }

//    ENDL(5);
}

TEST_F(EOFixture, rule_11 ){
    CE ce(*cl);

    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        ce.setCmpRules({11});
        auto res = ce.getExpansionOrder({0, 1, 2}, false);

        VI expected = {0, 1, 2, 6, 12, 3, 8, 7, 4, 10, 9, 5, 11};
        EXPECT_TRUE(equal(ALL(res.ord), ALL(expected)));
    }

//    ENDL(5);
}