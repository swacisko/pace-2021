//
// Created by sylwester on 3/26/21.
//


#include "CONTESTS/PACE21/heur/SwapCandidates/SwpCndEdge.h"
#include "CONTESTS/PACE21/test_graphs.h"
#include "gtest/gtest.h"

class SwpCndEdgeFixture : public ::testing::Test {
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

ClusterGraph* SwpCndEdgeFixture::clg = nullptr;
streambuf* SwpCndEdgeFixture::old_clog_buf = nullptr;
State* SwpCndEdgeFixture::st = nullptr;

TEST_F(SwpCndEdgeFixture, test1) {
    ClusterGraph &clg = (*SwpCndEdgeFixture::clg);
    State &st = (*SwpCndEdgeFixture::st);

    SwpCndEdgeCreator cr(st);
    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        {
            cr.keep_only_nonpositive_candidates = false;
            cr.keep_only_best_cluster_to_move_to = false;
            auto candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[1]]);
//    DEBUG(candidates);
            ASSERT_EQ(candidates.size(), 5);

            bool only_common_neighbors = false;
            candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[1]], only_common_neighbors);
//    DEBUG(candidates);
            ASSERT_EQ(candidates.size(), 11);

            ASSERT_EQ(cr.edges_to_cluster, VI(cr.edges_to_cluster.size(), 0));
        }

        {
            ENDL(10);
            cr.keep_only_nonpositive_candidates = false;
            cr.keep_only_best_cluster_to_move_to = false;
//            DEBUG(st.clusters[st.inCl[3]]);
            auto candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[3]]);
//            DEBUG(candidates);
            ASSERT_EQ(candidates.size(), 2);

            bool only_common_neighbors = false;
            cr.keep_only_nonpositive_candidates = true;
            candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[3]], only_common_neighbors);
//            DEBUG(candidates);
            ASSERT_EQ(candidates.size(),1 );

            ASSERT_EQ(cr.edges_to_cluster, VI(cr.edges_to_cluster.size(), 0));
        }
    }
}

TEST_F(SwpCndEdgeFixture, test2) {
    ClusterGraph &clg = (*SwpCndEdgeFixture::clg);
    State &st = (*SwpCndEdgeFixture::st);

    SwpCndEdgeCreator cr(st);
    cr.keep_only_nonpositive_candidates = true;
    cr.keep_only_best_cluster_to_move_to = false;
    for( int i=0; i<5; i++ ){ // checking for possible uncleared arrays
        auto candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[1]]);
//    DEBUG(candidates);
        ASSERT_EQ(candidates.size(), 3);

        bool only_common_neighbors = false;
        candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[1]], only_common_neighbors);
        ASSERT_EQ(candidates.size(), 3);

        ASSERT_EQ(cr.edges_to_cluster, VI(cr.edges_to_cluster.size(), 0));
    }
}

TEST_F(SwpCndEdgeFixture, test3) {
    ClusterGraph &clg = (*SwpCndEdgeFixture::clg);
    State &st = (*SwpCndEdgeFixture::st);

    SwpCndEdgeCreator cr(st);
    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        cr.keep_only_nonpositive_candidates = true;
        cr.keep_only_best_cluster_to_move_to = true;
        auto candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[1]]);
//    DEBUG(candidates);
        ASSERT_EQ(candidates.size(), 2);

        bool only_common_neighbors = false;
        candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[1]], only_common_neighbors);
        ASSERT_EQ(candidates.size(), 2);

        ASSERT_EQ(cr.edges_to_cluster, VI(cr.edges_to_cluster.size(), 0));
    }
}

TEST_F(SwpCndEdgeFixture, test4) {
    ClusterGraph &clg = (*SwpCndEdgeFixture::clg);
    State &st = (*SwpCndEdgeFixture::st);

    SwpCndEdgeCreator cr(st);
    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        cr.keep_only_nonpositive_candidates = false;
        cr.keep_only_best_cluster_to_move_to = true;
        auto candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[1]]);
//    DEBUG(candidates);
        ASSERT_EQ(candidates.size(), 3);

        bool only_common_neighbors = false;
        candidates = cr.create_MoveTo_SwapCandidatesForCluster(st.clusters[st.inCl[1]], only_common_neighbors);
        ASSERT_EQ(candidates.size(), 3);

        ASSERT_EQ(cr.edges_to_cluster, VI(cr.edges_to_cluster.size(), 0));
    }
}

TEST_F(SwpCndEdgeFixture, test_intercluster_edges) {
    ClusterGraph &clg = (*SwpCndEdgeFixture::clg);
    State &st = (*SwpCndEdgeFixture::st);

    SwpCndEdgeCreator cr(st);
    int empty = st.getIdOfEmptyCluster();

    for( int i=0; i<5; i++ ) { // checking for possible uncleared arrays
        {
            cr.keep_only_nonpositive_candidates = true;
            cr.keep_only_best_cluster_to_move_to = false;


            bool only_common_neighbors = false;
            auto candidates = cr.create_MoveTo_SwapCandidates_DifferentClusters(only_common_neighbors);
//        DEBUG(candidates);
//        DEBUG(candidates.size());
            ASSERT_EQ(candidates.size(), 16);

            int cnt = 0;
            for (auto &cnd : candidates) {
                if (cnd.getNodes() == VI({0, 3}) && cnd.getMoveNodesTo() == VI({1, 1})) {
                    cnt++;
                    ASSERT_EQ(cnd.swpVal(), -3);
                }
                if (cnd.getNodes() == VI({0, 3}) && cnd.getMoveNodesTo() == VI({empty, empty})) {
                    cnt++;
                    ASSERT_EQ(cnd.swpVal(), -2);
                }
                if (cnd.getNodes() == VI({0, 5}) && cnd.getMoveNodesTo() == VI({1, 1})) {
                    cnt++;
                    ASSERT_EQ(cnd.swpVal(), -12);
                }
                if (cnd.getNodes() == VI({2, 8}) && cnd.getMoveNodesTo() == VI({empty, empty})) {
                    cnt++;
                    int before = 22;
                    int after = 11;
                    ASSERT_EQ(cnd.swpVal(), after - before);
                }
            }

            ASSERT_EQ(cnt, 4); // checking just for existence of these two
            ASSERT_EQ(cr.edges_to_cluster, VI(cr.edges_to_cluster.size(), 0));
        }

        {
            cr.keep_only_nonpositive_candidates = true;
            cr.keep_only_best_cluster_to_move_to = true;

            int empty = st.getIdOfEmptyCluster();

            bool only_common_neighbors = false;
            auto candidates = cr.create_MoveTo_SwapCandidates_DifferentClusters(only_common_neighbors);
//        DEBUG(candidates);
//        DEBUG(candidates.size());
            ASSERT_EQ(candidates.size(), 10);

            int cnt = 0;
            for (auto &cnd : candidates) {
                if (cnd.getNodes() == VI({0, 3}) && cnd.getMoveNodesTo() == VI({1, 1})) {
                    cnt++;
                    ASSERT_EQ(cnd.swpVal(), -3);
                }
                if (cnd.getNodes() == VI({0, 3}) && cnd.getMoveNodesTo() == VI({empty, empty})) {
                    cnt++; // this should not happen, since we keepn only best candidates here
                    ASSERT_EQ(cnd.swpVal(), -2);
                }
                if (cnd.getNodes() == VI({2, 8}) && cnd.getMoveNodesTo() == VI({empty, empty})) {
                    cnt++;
                    ASSERT_EQ(cnd.swpVal(), -11);
                }
                if (cnd.getNodes() == VI({1, 8}) && cnd.getMoveNodesTo() == VI({empty, empty})) {
                    cnt++;
                    ASSERT_EQ(cnd.swpVal(), -13);
                }
            }

            ASSERT_EQ(cnt, 3); // checking just for existence of these two
            ASSERT_EQ(cr.edges_to_cluster, VI(cr.edges_to_cluster.size(), 0));
        }
    }
}
