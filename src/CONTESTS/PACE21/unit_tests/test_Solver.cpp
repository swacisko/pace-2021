//
// Created by sylwester on 4/15/21.
//

#include <graphs/GraphUtils.h>
#include <CONTESTS/PACE21/heur/PaceUtils.h>
#include <CONTESTS/PACE21/heur/SwapCandidates/SwpCndNode.h>
#include <CONTESTS/PACE21/heur/SwapCandidates/SwpCndTriangle.h>
#include "CONTESTS/PACE21/heur/Solver.h"
#include "CONTESTS/PACE21/test_graphs.h"
#include "gtest/gtest.h"

class SolverFixture : public ::testing::Test {
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

//        clog << "Cluster created" << endl;
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

ClusterGraph* SolverFixture::clg = nullptr;
streambuf* SolverFixture::old_clog_buf = nullptr;
State* SolverFixture::st = nullptr;
VVI SolverFixture::V;


TEST_F(SolverFixture, createPartitionForGivenState){

    auto edges = CE_test_graphs::swpcndedge_test2_edges;
    VVI V = CE_test_graphs::createVVIfromVPII(edges);
    VI partition = CE_test_graphs::swpcndedge_test2_partition;

    Config cnf;
    Solver solver( V, partition, cnf );

    VI part = solver.createPartitionsForGivenState(*st).first;

//    DEBUG(part);

    for( int i=0; i<=7; i++ ) ASSERT_EQ( part[i], part[0] );
    for( int i=8; i<=10; i++ ) ASSERT_EQ( part[i], part[8] );
    for( int i=11; i<=14; i++ ) ASSERT_EQ( part[i], part[11] );
    for( int i=15; i<=20; i++ ) ASSERT_EQ( part[i], part[15] );
}


TEST_F(SolverFixture, createGranularity1){

    auto edges = CE_test_graphs::swpcndedge_test2_edges;
    VVI V = CE_test_graphs::createVVIfromVPII(edges);
    VI partition = CE_test_graphs::swpcndedge_test2_partition;

    Config cnf;
    cnf.coarsen_mode = (contract_all | remove_edges);

    VI initial_partition(V.size(),0);
    iota(ALL(initial_partition),0);
    Solver solver( V, initial_partition, cnf );

    solver.createClusterGraph();

    ASSERT_EQ(solver.V, V);

    partition = solver.createPartitionsForGivenState(*st).first;
//    DEBUG(partition);

    solver.known_solutions.push_back(partition);
    // now clusters are {0,1,2}, {3,4}, {5,6}, {7,8}

    {
        VI part = partition;
        for( int i=11; i<=14; i++ ) part[i] = part[0]; // moving nodes {5,6} of clg to cluster {0,1,2}
        solver.known_solutions.push_back( part );
        // now clusters are {0,1,2,5,6}, {3,4}, {7,8}
    }

    {
        VI part = partition;
        for( int i=15; i<=20; i++ ) part[i] = part[0]; // moving nodes {7,8} of clg to cluster {0,1,2}
        solver.known_solutions.push_back( part );
        // now clusters are {0,1,2,7,8}, {3,4}, {5,6}
    }

    {
        VI part = partition;
        for( int i=11; i<=14; i++ ) part[i] = part[8]; // moving nodes {5,6} of clg to cluster {3,4}
        solver.known_solutions.push_back( part );
        // now clusters are {0,1,2}, {3,4,5,6}, {7,8}
    }

    {
        VI part = partition;
        for( int i=13; i<=14; i++ ) part[i] = part[15]; // moving node {6} of clg to cluster {7,8}
        for( int i=18; i<=20; i++ ) part[i] = part[0]; // moving nodes {8} of clg to cluster {0,1,2}
        // now clusters are {0,1,2,8}, {3,4}, {5}, {6,7}
        solver.known_solutions.push_back( part );
    }

    // now, nodes {3,4} of clg, that is {8,9,10} of V are always in the same cluster.
    // nodes {0,1,2} are always in the same cluster
    // edges between nodes {8,9,10} and nodes from {0,...,7} have always both ends in different clusters, hence they
    // should be removed


//    clog << "known solutions:" << endl;
//    for(VI & v : solver.known_solutions) clog << v << endl;

    solver.granulateSolution();

//    DEBUG(solver.partition);

    for( int i=0; i<=7; i++ ) ASSERT_EQ(solver.partition[i], solver.partition[0]);
    for( int i=8; i<=10; i++ ) ASSERT_EQ(solver.partition[i], solver.partition[8]);
    for( int i=11; i<=12; i++ ) ASSERT_EQ(solver.partition[i], solver.partition[11]);
    for( int i=13; i<=14; i++ ) ASSERT_EQ(solver.partition[i], solver.partition[13]);
    for( int i=15; i<=17; i++ ) ASSERT_EQ(solver.partition[i], solver.partition[15]);
    for( int i=18; i<=20; i++ ) ASSERT_EQ(solver.partition[i], solver.partition[18]);

    ASSERT_NE( solver.partition[11], solver.partition[13] );
    ASSERT_NE( solver.partition[15], solver.partition[18] );
    ASSERT_NE( solver.partition[0], solver.partition[18] );

    ASSERT_EQ( GraphUtils::countEdges(V), GraphUtils::countEdges(solver.V) + 6 ); // 6 edges are removed
}

TEST_F(SolverFixture, createGranularity2){

    auto edges = CE_test_graphs::swpcndedge_test2_edges;
    VVI V = CE_test_graphs::createVVIfromVPII(edges);
    VI partition = CE_test_graphs::swpcndedge_test2_partition;

    Config cnf;
    cnf.coarsen_mode = (contract_all | remove_edges);

    VI initial_partition(V.size(),0);
    iota(ALL(initial_partition),0);
    Solver solver( V, initial_partition, cnf );
    solver.createClusterGraph();

    ASSERT_EQ(solver.V, V);

    partition = solver.createPartitionsForGivenState(*st).first;
//    DEBUG(partition);

    solver.known_solutions.push_back(partition);
    // now clusters are {0,1,2}, {3,4}, {5,6}, {7,8}

    {
        VI part = partition;
        for( int i=11; i<=14; i++ ) part[i] = part[0]; // moving nodes {5,6} of clg to cluster {0,1,2}
        solver.known_solutions.push_back( part );
        // now clusters are {0,1,2,5,6}, {3,4}, {7,8}
    }

    {
        VI part = partition;
        for( int i=11; i<=14; i++ ) part[i] = part[8]; // moving nodes {5,6} of clg to cluster {3,4}
        solver.known_solutions.push_back( part );
        // now clusters are {0,1,2}, {3,4,5,6}, {7,8}
    }


    // now, nodes {3,4} of clg, that is {8,9,10} of V are always in the same cluster.
    // nodes {0,1,2} of clg are always in the same cluster
    // nodes {5,6} of clg are always in the same cluster
    // nodes {7,8} of clg are always in the same cluster
    // edges between nodes {8,9,10} and nodes from {0,...,7} have always both ends in different clusters, hence they
    // should be removed. The same is true for edges between {15-20} and {0-7}


//    clog << "known solutions:" << endl;
//    for(VI & v : solver.known_solutions) clog << v << endl;

    solver.granulateSolution();

//    DEBUG(solver.partition);

    for( int i=0; i<=7; i++ ) ASSERT_EQ(solver.partition[i], solver.partition[0]);
    for( int i=8; i<=10; i++ ) ASSERT_EQ(solver.partition[i], solver.partition[8]);
    for( int i=11; i<=14; i++ ) ASSERT_EQ(solver.partition[i], solver.partition[11]);
    for( int i=15; i<=20; i++ ) ASSERT_EQ(solver.partition[i], solver.partition[15]);

    ASSERT_NE( solver.partition[11], solver.partition[18] );
    ASSERT_NE( solver.partition[0], solver.partition[18] );

    ASSERT_EQ( GraphUtils::countEdges(V), GraphUtils::countEdges(solver.V) + 26 ); // 28 edges are removed
}


TEST_F(SolverFixture, testOverall_1){

    Config cnf;
//    cnf.coarsen_mode = (contract_all | remove_edges);
    cnf.coarsen_mode = contract_matching;
    cnf.granularity_frequency = 20;

    VI initial_partition = CE_test_graphs::swpcndedge_test2_partition;
    Solver solver( V, initial_partition, cnf );

//    clog << "Starting solver.run()" << endl;
    solver.run();

    VVI res_cl = {
            {0,1,8,9,10,11,12}, {2,3,4}, {15,16,17}, {5,6,7,18,19,20}
    };


    DEBUG(solver.best_partition);
    DEBUG(PaceUtils::partitionToClusters(solver.best_partition));

    for( int i=0; i<res_cl.size(); i++ ){
        for( int d : res_cl[i] ) ASSERT_EQ( solver.best_partition[d], solver.best_partition[res_cl[i][0]] );
        if(i>0) ASSERT_NE( solver.best_partition[ res_cl[i-1][0] ], solver.best_partition[ res_cl[i][0] ] );
    }

}

TEST_F(SolverFixture, testOverall_2) {
    Config cnf;
    cnf.granularity_frequency = 10;

    VVI V = CE_test_graphs::createVVIfromVPII(CE_test_graphs::cluster_graph_test2_edges);
    VI initial_partition(V.size());
    iota(ALL(initial_partition), 0);
    Solver solver(V, initial_partition, cnf);

    clog << "Starting solver.run()" << endl;
    solver.run();

    VVI best_sol_clusters = StandardUtils::partitionToLayers(solver.best_partition);
    for (VI &v : best_sol_clusters) sort(ALL(v));
    sort(ALL(best_sol_clusters), [](VI &v1, VI &v2) { return v1[0] < v2[0]; });

    clog << "Best solution clusters: " << best_sol_clusters << endl;

    ASSERT_EQ(PaceUtils::evaluateSolution(V, solver.best_partition), 18);

    VVVI res_cl_pos = {
            {{3,4,5}, {0,1,6,7,8}, {9},{10, 11, 13, 14}, {12},{15}, {2}},
            {{2},{12},{0, 1, 6, 7}, {3,  4,  5},{8, 9, 10, 11}, {13, 14, 15}},
            {{{0, 1, 6, 7, 8}, {2}, {3, 4, 5}, {9}, {10, 13, 14, 15}, {11, 12}}},
            {{0, 1, 6, 7, 8}, {2}, {3, 4, 5}, {9}, {10, 11, 13, 14}, {12, 15}},
            {{0, 1, 6, 7, 8}, {2}, {3, 4, 5}, {9, 10}, {11, 12}, {13, 14, 15}},
            {{0, 1, 6, 7}, {2}, {3, 4, 5}, {8, 9, 10}, {11, 12}, {13, 14, 15}},
            {{0, 1, 2, 6, 7}, {3, 4, 5}, {8, 9, 10, 11}, {12}, {13, 14, 15}}
    };


    for (VVI &x : res_cl_pos) {
        for (VI &v : x) sort(ALL(v));
        sort(ALL(x), [](VI &v1, VI &v2) { return v1[0] < v2[0]; });
    }

    bool solution_in_set_of_optimals = false;
    for (VVI &x : res_cl_pos) if (x == best_sol_clusters) solution_in_set_of_optimals = true;

    ASSERT_TRUE(solution_in_set_of_optimals);
}

TEST_F(SolverFixture, testOverall_3){

    Config cnf;
    cnf.granularity_frequency = 10;

    VVI V = CE_test_graphs::createVVIfromVPII(CE_test_graphs::kern_heur_7_edges);
    VI initial_partition(V.size());
    iota(ALL(initial_partition),0);
    Solver solver( V, initial_partition, cnf );

    clog << "Starting solver.run()" << endl;
    solver.run();

    VVI best_sol_clusters = StandardUtils::partitionToLayers(solver.best_partition);
    for( VI & v : best_sol_clusters ) sort(ALL(v));
    sort(ALL(best_sol_clusters), []( VI& v1, VI& v2 ){ return v1[0] < v2[0]; });

    clog << "Best solution clusters: " << best_sol_clusters << endl;

    ASSERT_EQ( PaceUtils::evaluateSolution( V, solver.best_partition ), 13 );

    VVVI res_cl_pos = {
            { {0,1,3}, {4,5,6,7,8}, {10,9,14}, {11,12,13,15}, {16,19,18},{2,17} },
            { {0,1,3}, {4,5,6,7,8}, {10,9,14,11}, {12,13,15}, {16,19,18},{2,17} },

            {{0,1,2,3}, {4,5,6,7,8}, {10,9,14}, {11,12,13,15}, {16,19,18},{17}},
            {{0,1,2,3}, {4,5,6,7,8}, {10,9,14,11}, {12,13,15}, {16,19,18},{17}},
    };

    for( VVI & x : res_cl_pos ){
        for( VI & v : x ) sort(ALL(v));
        sort(ALL(x), []( VI& v1, VI& v2 ){ return v1[0] < v2[0]; });
    }

    bool solution_in_set_of_optimals = false;
    for( VVI & x : res_cl_pos ) if( x == best_sol_clusters ) solution_in_set_of_optimals = true;

    ASSERT_TRUE(solution_in_set_of_optimals);
}

TEST_F(SolverFixture, testApplySwap){

    Config cnf;
    cnf.granularity_frequency = 10;

    VI initial_partition = CE_test_graphs::swpcndedge_test2_partition;

    Solver solver( V, initial_partition, cnf );

    solver.st = st;
    solver.clg = *clg;

    vector<SwapCandidate*> candidates;
    SwpCndNode cnd1( -1, 0, st->getIdOfEmptyCluster() );
    SwpCndNode cnd2( -2, 3, st->getIdOfEmptyCluster() );

    class TempCnd : public SwapCandidateAdapter{
    public:
        VI v = {5,7};
        VI w = {4,4};
        virtual VI getNodes()override{ return v; }
        virtual VI getMoveNodesTo() override{ return w; }
        virtual LL swpVal()override{ return -3; }
    } cnd3;

    candidates.push_back( &cnd1 );
    candidates.push_back( &cnd2 );
    candidates.push_back( &cnd3 );

//    DEBUG(*solver.st);

    solver.apply_swap_for_candidates(candidates);

//    clog << "After" << endl;
//    DEBUG(*solver.st);
    ASSERT_NE( solver.st->inCl[0], solver.st->inCl[3] );
    ASSERT_EQ( solver.st->inCl[5], solver.st->inCl[7] );

    solver.st = nullptr; // this will be deleted in Fixture::TearDown
}