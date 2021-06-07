//
// Created by sylwester on 3/25/21.
//

#include "CONTESTS/PACE21/heur/SwapCandidates/SwpCndNode.h"
#include "CONTESTS/PACE21/test_graphs.h"
#include "gtest/gtest.h"


class SwpCndNodeFixture : public ::testing::Test {
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

ClusterGraph* SwpCndNodeFixture::clg = nullptr;
streambuf* SwpCndNodeFixture::old_clog_buf = nullptr;

TEST_F(SwpCndNodeFixture, test_all_onlybest){
    ClusterGraph & clg = (*SwpCndNodeFixture::clg);
    State st(clg, SINGLE_NODES);
    VVI to_merge = { {0,4}, {7,8}, {2,6,3} }; // the only single nodes remain 1 and 5
    st.mergeClusters(to_merge);

    // clusters after merging: {1}, {5}, {0,4}, {2,6,3}, {7,8}

    SwpCndNodeCreator cr(st);

    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        cr.keep_only_best_cluster_to_move_to = false;
        cr.keep_only_nonpositive_candidates = false;

        auto candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[3]]);

        VVI cnt(st.N);
        for (auto cnd : candidates) {
            cnt[cnd.getNodesToSwap()[0].first].push_back(cnd.getMoveNodesTo()[0]);
        }

        for (auto &v : cnt) sort(ALL(v));

        const int empty = st.getIdOfEmptyCluster();
        ASSERT_EQ(cnt[2], VI({0, empty}));
        ASSERT_EQ(cnt[3], VI({1, 2, 4, empty}));
        ASSERT_EQ(cnt[6], VI({4, empty}));
    }
}

TEST_F(SwpCndNodeFixture, test_nonpos1){
    ClusterGraph & clg = (*SwpCndNodeFixture::clg);
    State st(clg, SINGLE_NODES);
    VVI to_merge = { {0,4}, {7,8}, {2,6,3} };
    st.mergeClusters(to_merge);

    // clusters after merging: {1}, {5}, {0,4}, {2,6,3}, {7,8}

    SwpCndNodeCreator cr(st);

    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        {
            cr.keep_only_best_cluster_to_move_to = false;
            cr.keep_only_nonpositive_candidates = false;

            auto candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[8]]);
            ASSERT_EQ(candidates.size(), 5);
        }

        {
//            SwpCndNodeCreator cr(st);
            cr.keep_only_best_cluster_to_move_to = false;
            cr.keep_only_nonpositive_candidates = true;

            auto candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[8]]);
            ASSERT_TRUE(candidates.empty());
        }
    }
}


TEST_F(SwpCndNodeFixture, test_nonpos2){
    ClusterGraph & clg = (*SwpCndNodeFixture::clg);
    State st(clg, SINGLE_NODES);
    VVI to_merge = {  {4,5}, {7,8}, {2,6,3} };
    st.mergeClusters(to_merge);

//    DEBUG(st.clusters);
    // clusters after merging: {0}, {4,5}, {2,6,3}, {7,8}


    SwpCndNodeCreator cr(st);
    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        {
            cr.keep_only_best_cluster_to_move_to = false;
            cr.keep_only_nonpositive_candidates = false;

            auto candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[4]]);

            ASSERT_EQ(candidates.size(), 5);
        }

        {
//            SwpCndNodeCreator cr(st);
            cr.keep_only_best_cluster_to_move_to = false;
            cr.keep_only_nonpositive_candidates = true;

            auto candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[4]]);

            ASSERT_EQ(candidates.size(), 1);
            auto cnd = candidates[0];
            ASSERT_EQ(cnd.getNodesToSwap()[0], PII({4, 0}));
            ASSERT_EQ(cnd.swpVal(), 0);
        }
    }
}

TEST_F(SwpCndNodeFixture, test_best){
    ClusterGraph & clg = (*SwpCndNodeFixture::clg);
    State st(clg, SINGLE_NODES);
    VVI to_merge = {  {4,5}, {7,8}, {2,6,3} };
    st.mergeClusters(to_merge);

//    DEBUG(st.clusters);
    // clusters after merging: {0}, {4,5}, {2,6,3}, {7,8}

    SwpCndNodeCreator cr(st);
    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        {
            cr.keep_only_best_cluster_to_move_to = false;
            cr.keep_only_nonpositive_candidates = false;

            auto candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[4]]);
            ASSERT_EQ(candidates.size(), 5);
        }

        {
//            SwpCndNodeCreator cr(st);
            cr.keep_only_best_cluster_to_move_to = true;
            cr.keep_only_nonpositive_candidates = false;

            auto candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[4]]);
            ASSERT_EQ(candidates.size(), 2);
        }

        {
//            SwpCndNodeCreator cr(st);
            cr.keep_only_best_cluster_to_move_to = true;
            cr.keep_only_nonpositive_candidates = true;

            auto candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[8]]);
            ASSERT_EQ(candidates.size(), 0);
        }
    }
}