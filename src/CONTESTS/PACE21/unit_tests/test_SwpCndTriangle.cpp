//
// Created by sylwester on 3/29/21.
//


#include "CONTESTS/PACE21/heur/SwapCandidates/SwpCndTriangle.h"
#include "CONTESTS/PACE21/test_graphs.h"
#include "gtest/gtest.h"


class SwpCndTriangleFixture : public ::testing::Test {
public:

    static void SetUpTestSuite() {
        old_clog_buf = clog.rdbuf();
        clog.rdbuf(cout.rdbuf());

        auto edges = CE_test_graphs::swpcndedge_test2_edges;
        VVI V = CE_test_graphs::createVVIfromVPII(edges);
        VI partition = CE_test_graphs::swpcndedge_test2_partition;

        clg = new ClusterGraph(&V, partition);

        st = new State(*clg, SINGLE_NODES);
        VVI to_merge = {{0, 1, 2},
                        {3, 4},
                        {5, 6},
                        {7, 8}}; // the only single nodes remain 1 and 5
        st->mergeClusters(to_merge);

//        clog << "Cluster created" << endl;
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

ClusterGraph* SwpCndTriangleFixture::clg = nullptr;
streambuf* SwpCndTriangleFixture::old_clog_buf = nullptr;
State* SwpCndTriangleFixture::st = nullptr;


TEST_F(SwpCndTriangleFixture, test1){
    ClusterGraph &clg = (*SwpCndTriangleFixture::clg);
    State &st = (*SwpCndTriangleFixture::st);

    SwpCndTriangleCreator cr(st);

    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        cr.keep_only_nonpositive_candidates = false;
        cr.keep_only_best_cluster_to_move_to = false;
        auto candidates = cr.create_MoveTo_SwapCandidates(true);

    //    clog << "Candidates:" << endl;
    //    for(auto &cnd : candidates ) clog << cnd << endl;


        ASSERT_EQ(candidates.size(), 7);
        for (auto &cnd : candidates) {
            VI nodes = cnd.getNodes();
            sort(ALL(nodes));
            if (nodes == VI({0, 1, 3})) ASSERT_EQ(cnd.swpVal(), 0);
//            if (nodes == VI({0, 1, 2})) ASSERT_EQ(cnd.swpVal(), 0);
            if (nodes == VI({0, 3, 4})) ASSERT_EQ(cnd.swpVal(), -4);
            if (nodes == VI({0, 3, 5})) ASSERT_EQ(cnd.swpVal(), -6);
            if (nodes == VI({0, 4, 5})) ASSERT_EQ(cnd.swpVal(), -8);
            if (nodes == VI({1, 2, 7})) ASSERT_EQ(cnd.swpVal(), -9);
            if (nodes == VI({1, 2, 8})) ASSERT_EQ(cnd.swpVal(), -13);
            if (nodes == VI({3, 4, 5})) ASSERT_EQ(cnd.swpVal(), -6);
        }


        { // this test is the copy-paste test from test2. It is here to check, whether data is correctly cleared.
            cr.keep_only_nonpositive_candidates = false;
            cr.keep_only_best_cluster_to_move_to = false;
            auto candidates = cr.create_MoveTo_SwapCandidates(false);
            ASSERT_EQ(candidates.size(), 16);
            int cnt = 0, cnt2 = 0, cnt3 = 0;
            for (auto &cnd : candidates) {
                VI nodes = cnd.getNodes();
                sort(ALL(nodes));
                if (nodes == VI({0, 1, 2})) cnt++;
                if (nodes == VI({0, 1, 3})) cnt2++;
                if (nodes == VI({3, 4, 5})) cnt3++;
            }

            ASSERT_EQ(cnt, 3);
            ASSERT_EQ(cnt2, 3);
            ASSERT_EQ(cnt3, 2);
        }
    }
}

TEST_F(SwpCndTriangleFixture, test2){
    ClusterGraph &clg = (*SwpCndTriangleFixture::clg);
    State &st = (*SwpCndTriangleFixture::st);

    SwpCndTriangleCreator cr(st);

    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        cr.keep_only_nonpositive_candidates = false;
        cr.keep_only_best_cluster_to_move_to = false;
        auto candidates = cr.create_MoveTo_SwapCandidates(false);

    //    clog << "Candidates:" << endl;
    //    for(auto &cnd : candidates ) clog << cnd << endl;

        ASSERT_EQ(candidates.size(), 16);

        int cnt = 0, cnt2 = 0, cnt3 = 0;
        for (auto &cnd : candidates) {
            VI nodes = cnd.getNodes();
            sort(ALL(nodes));
            if (nodes == VI({0, 1, 2})) cnt++;
            if (nodes == VI({0, 1, 3})) cnt2++;
            if (nodes == VI({3, 4, 5})) cnt3++;
        }

        ASSERT_EQ(cnt, 3);
        ASSERT_EQ(cnt2, 3);
        ASSERT_EQ(cnt3, 2);
    }
}

TEST_F(SwpCndTriangleFixture, test3){
    ClusterGraph &clg = (*SwpCndTriangleFixture::clg);
    State &st = (*SwpCndTriangleFixture::st);

    SwpCndTriangleCreator cr(st);

    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        cr.keep_only_nonpositive_candidates = true;
        cr.keep_only_best_cluster_to_move_to = true;
        auto candidates = cr.create_MoveTo_SwapCandidates(false);

//        clog << "Candidates:" << endl;
//        for(auto &cnd : candidates ) clog << cnd << endl;

        ASSERT_EQ(candidates.size(), 7);
    }
}