//
// Created by sylwester on 3/25/21.
//

#include "CONTESTS/PACE21/heur/SwapCandidates/SwapCandidate.h"
#include "CONTESTS/PACE21/test_graphs.h"
#include "gtest/gtest.h"

class SwapCandidateFixture : public ::testing::Test {
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

ClusterGraph* SwapCandidateFixture::clg = nullptr;
streambuf* SwapCandidateFixture::old_clog_buf = nullptr;

TEST_F(SwapCandidateFixture, update_edges){
    ClusterGraph & clg = (*SwapCandidateFixture::clg);
    State st(clg, SINGLE_NODES);
    VVI to_merge = { {0,1,2},{4,5},{6,7,8} };
    st.mergeClusters(to_merge);

//    DEBUG(st.clusters);

    SwapCandidateCreatorAdapter swpcnd_ad(st);
    swpcnd_ad.keep_only_best_cluster_to_move_to = false;

    VI neigh;
    swpcnd_ad.updateEdgesInNeighboringClustersForNode( 3, st.inCl[3], neigh );
    sort(ALL(neigh));

    ASSERT_EQ( neigh, VI({1,2,3}) );
    ASSERT_EQ( swpcnd_ad.edges_to_cluster[0], 0 );
    ASSERT_EQ( swpcnd_ad.edges_to_cluster[1], 5 );
    ASSERT_EQ( swpcnd_ad.edges_to_cluster[2], 4 );
    ASSERT_EQ( swpcnd_ad.edges_to_cluster[3], 4 );
}

TEST_F(SwapCandidateFixture, swapval1){
    ClusterGraph & clg = (*SwapCandidateFixture::clg);
    State st(clg, SINGLE_NODES);
    VVI to_merge = { {0,1,2},{4,5},{6,7,8} };
    st.mergeClusters(to_merge);

//    DEBUG(st);

    { // swap value for moving node 3 to cluster {0,1,2}
        SwapCandidateCreatorAdapter swpcnd_ad(st);
        swpcnd_ad.keep_only_best_cluster_to_move_to = false;

        VI neigh;
        swpcnd_ad.updateEdgesInNeighboringClustersForNode( 3, st.inCl[3], neigh );

        int swp_val1 = 0;
        int swpcnd_w = clg.node_weights[3];
        int swp_trg_cl = st.inCl[0];
        int deg_in_cl = 0;
        int swap_value = swpcnd_ad.getSwpValForMove(swp_val1, swpcnd_w, deg_in_cl, swp_trg_cl);

        int before = 0 + 5; // edges to add/remove to make {3} a cluster (considering only clusters {3} and {0,1,2})
        int after = 4 * (3+2+1) - 5; // edges to add/remove to make 3 stay in cluster {0,1,2}: add (25-5), remove 0
        ASSERT_EQ( swap_value, after - before );
    }
}


TEST_F(SwapCandidateFixture, swapval2){
    ClusterGraph & clg = (*SwapCandidateFixture::clg);
    State st(clg, SINGLE_NODES);
    VVI to_merge = { {0,1,2},{4,5},{6,7,8} };
    st.mergeClusters(to_merge);

//    DEBUG(st);

    { // swap value for moving node 1 to cluster {6,7,8}
        SwapCandidateCreatorAdapter swpcnd_ad(st);
        swpcnd_ad.keep_only_best_cluster_to_move_to = false;

        VI neigh;
        swpcnd_ad.updateEdgesInNeighboringClustersForNode( 1, st.inCl[1], neigh );

        int swpcnd_w = clg.node_weights[1];
        int swp_val1 = swpcnd_w * (3+1) - (2+1);
        int swp_trg_cl = st.inCl[6];
        int deg_in_cl = 3;
        int swap_value = swpcnd_ad.getSwpValForMove(swp_val1, swpcnd_w, deg_in_cl, swp_trg_cl);

        int before = swp_val1 + 1; // edges to add/remove to move 1 to cluster {6,7,8}
        int after = (swpcnd_w * ( 1+2+1 ) - 1) + 3; // edges to add/remove to make 1 stay in cluster {6,7,8}
        ASSERT_EQ( swap_value, after - before );
    }
}

TEST_F(SwapCandidateFixture, swapval3){
    ClusterGraph & clg = (*SwapCandidateFixture::clg);
    State st(clg, SINGLE_NODES);
    VVI to_merge = { {0,3,4,5},{1,2,6,7,8} };
    st.mergeClusters(to_merge);

//    DEBUG(st);

    { // swap value for moving node 3 from cluster {0,3,4,5} to {1,2,6,7,8}
        SwapCandidateCreatorAdapter swpcnd_ad(st);

        VI neigh;
        swpcnd_ad.updateEdgesInNeighboringClustersForNode( 3, st.inCl[3], neigh );

        int swpcnd_w = 4;
        int swp_val1 = 4 * ( 3+1+1 ) - (4+2+2);
        int swp_trg_cl = st.inCl[6];
        int deg_in_cl = (4+2+2);
        int swap_value = swpcnd_ad.getSwpValForMove(swp_val1, swpcnd_w, deg_in_cl, swp_trg_cl);

        int before = swp_val1 + (3+1+1); // edges to add/remove to move 3 to cluster {1,2,6,7,8}
        int after = ( 4 * ( 2+1+1+2+1 ) - (3+1+1) ) + (4+2+2); // edges to add/remove to make 3 stay in cluster {1,2,6,7,8}
        ASSERT_EQ( swap_value, after - before );
    }
}

TEST_F(SwapCandidateFixture, swapval4){
    ClusterGraph & clg = (*SwapCandidateFixture::clg);
    State st(clg, SINGLE_NODES);
    VVI to_merge = { {0,3,4,5},{1,2,6,7,8} };
    st.mergeClusters(to_merge);

//    DEBUG(st);

    { // swap value for moving node 4 from cluster {0,3,4,5} to {1,2,6,7,8} - no edges between
        SwapCandidateCreatorAdapter swpcnd_ad(st);

        VI neigh;
        swpcnd_ad.updateEdgesInNeighboringClustersForNode( 4, st.inCl[4], neigh );

        int swpcnd_w = 1;
        int swp_val1 = 1*(3+4+1) - (2+2+1);
        int swp_trg_cl = st.inCl[6];
        int deg_in_cl = (2+2+1);
        int swap_value = swpcnd_ad.getSwpValForMove(swp_val1, swpcnd_w, deg_in_cl, swp_trg_cl);

        int before = swp_val1; // edges to add/remove to move 4 to cluster {1,2,6,7,8}
        int after = ( 1*( 2+1+1+2+1 ) ) + (2+2+1); // edges to add/remove to make 4 stay in cluster {1,2,6,7,8}
        ASSERT_EQ( swap_value, after - before );
    }
}

TEST_F(SwapCandidateFixture, swapval5){
    ClusterGraph & clg = (*SwapCandidateFixture::clg);
    State st(clg, SINGLE_NODES);
    VVI to_merge = { {0,3,4,5},{1,2,6,7,8} };
    st.mergeClusters(to_merge);


    { // swap value for moving node 0 to empty cluster
        SwapCandidateCreatorAdapter swpcnd_ad(st);

        VI neigh;
        swpcnd_ad.updateEdgesInNeighboringClustersForNode( 0, st.inCl[0], neigh );

        int swpcnd_w = 3;
        int swp_val1 = 3*(1+4+1) - (2+4);
        int swp_trg_cl = st.getIdOfEmptyCluster();
        int deg_in_cl = (2+4);
        int swap_value = swpcnd_ad.getSwpValForMove(swp_val1, swpcnd_w, deg_in_cl, swp_trg_cl);

        int before = swp_val1; // edges to add/remove to move 4 to cluster {1,2,6,7,8}
        int after = (2+4); // edges to add/remove to make 4 stay in cluster {1,2,6,7,8}
        ASSERT_EQ( swap_value, after - before );
    }
}