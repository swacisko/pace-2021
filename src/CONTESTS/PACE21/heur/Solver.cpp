//
// Created by sylwester on 4/14/21.
//

#include <graphs/GraphUtils.h>
#include <datastructures/FAU.h>
#include <utils/RandomNumberGenerators.h>
#include <CONTESTS/PACE21/heur/PaceUtils.h>
#include <CONTESTS/PACE21/heur/SwapCandidates/SwapCandidate.h>
#include <CONTESTS/PACE21/heur/ExpansionOrder.h>
#include <CONTESTS/PACE21/heur/SwapCandidates/ComponentExpansionAttraction.h>
#include <CONTESTS/PACE21/heur/SwapCandidates/ComponentExpansionRepulsion.h>
#include <CONTESTS/PACE21/heur/SwapCandidates/SwpCndNode.h>
#include <CONTESTS/PACE21/heur/SwapCandidates/SwpCndEdge.h>
#include <CONTESTS/PACE21/heur/SwapCandidates/SwpCndTriangle.h>
#include <CONTESTS/PACE21/heur/SwapCandidates/SwpCndEO.h>
#include <CONTESTS/PACE21/heur/Global.h>
#include <CONTESTS/PACE21/heur/StateImprovers/NodeEdgeGreedy.h>
#include <CONTESTS/PACE21/kernelization/CEKernelizer.h>
#include <CONTESTS/PACE21/heur/StateImprovers/SparseGraphTrimmer.h>
#include <CONTESTS/PACE21/heur/StateImprovers/NodeEdgeGreedyNomap.h>
#include <CONTESTS/PACE21/heur/StateImprovers/NodeEdgeGreedyW1.h>
#include "CONTESTS/PACE21/heur/Solver.h"

Solver::Solver(VVI & V, VI initial_partition, Config& cnf, int rec_depth){
    this->cnf = &cnf;
    this->origV = &V;
    this->V = V;
    N = V.size();
    partition = initial_partition;
    recurrence_depth = rec_depth;

    hashes = hashes2 = VLL(N);
    int SOME_SEED = 1233101;
    UniformIntGenerator rnd(0, 1ll * 1'000'000'000 * 1'000'000'000, SOME_SEED);
    for (int i = 0; i < N; i++){
        hashes[i] = rnd.rand();
        hashes2[i] = rnd.rand();
    }

    best_partition = initial_partition;
    best_result = PaceUtils::evaluateSolution(V,best_partition);
}

Solver::~Solver(){
    if(st != nullptr){
        delete st;
        st = nullptr;
    }
}

void Solver::run(int iters) {
    if(Global::checkTle()) return;

    const bool debug = ( recurrence_depth == 0 );
    const bool debug_all = false;

    if(recurrence_depth == 0 && cnf->use_kernelization){
        CEKernelizer kern( *origV );
        kern.fullKernelization(cnf->use_heuristic_kernelization,0);
        partition = kern.inCl;
        clog << "Kernelization done!" << endl;
    }

    int large_iter_cnt = 0;

    while( !Global::checkTle() && large_iter_cnt < iters ){ // CAUTION!! 'Infinite' loop now!!
        if(Global::checkTle()) return;

        large_iter_cnt++;
        if(debug) clog << endl << endl << "***************************  Large iteration #" << large_iter_cnt << endl;

        pair<VI,VI> new_results = largeIteration();

        known_solutions.push_back(new_results.first);
        known_clg_partitions.push_back(new_results.second);

        if(debug_all){
            clog << "State after large iteration represents clusters: " << endl;
            for( auto & cl : st->clusters){
                set<int> v;
                for( int x : cl.g.nodes ){
                    v += clg.clusterNodes[x];
                }
                clog << v << endl;
            }
        }

        compareToBestSolutionAndUpdate(new_results.first); // updating (hopefully) best solution
        if(Global::checkTle()) return;

        if( (large_iter_cnt % cnf->granularity_frequency) == 0 ){
            granulateSolution();
        }

        if(debug){
            clog << endl << "Current iteration result: " << PaceUtils::evaluateSolution(*origV, new_results.first) << endl;
            clog << "Best result so far: " << best_result << endl;
        }
    }

    if(!Global::disable_all_logs){
        clog << "Creators: (calls,improvements):" << endl;
        for( auto & [s,p] : local_search_creator_calls ){
            clog << s << " --> " << p << endl;
        }
    }

}

void Solver::run_recursive() {
    const bool debug = ( recurrence_depth <= 5 );

    if(recurrence_depth == 0 && cnf->use_kernelization){
        CEKernelizer kern( *origV );
//        kern.fullKernelization(cnf->use_heuristic_kernelization,0);

        if(cnf->use_only_fast_exact_kernelization){ // disable all rules except very fast rules
            for( int i=2; i<kern.MAX_RULES; i++ ) kern.setDisableRule(i,true);
            kern.setDisableRule(15,false);
            kern.setDisableRule(16,false);
            kern.setDisableRule(2,false);
        }

        kern.fullKernelization(false,0);
        if(cnf->use_heuristic_kernelization) kern.improveKernelizationUsingHeuristicRules();
        partition = kern.inCl;
    }

    createClusterGraph();

    if(debug){
        clog << "***********   depth: " << recurrence_depth << ", origV - nodes: " << origV->size() << ", edges: " << GraphUtils::countEdges(*origV) <<
             "   |   clg nodes: " << clg.V.size() << ", edges: " << GraphUtils::countEdges( clg.V ) << endl;
    }


    bool use_run_fast = true;
    int N = origV->size();
    int E = GraphUtils::countEdges(*origV);
//    if( N < 10'000 || E < 200'000 ) use_run_fast = false; // original
    if( N < 15'000 && E < 200'000 ) use_run_fast = false;
//    if( N < 15'000 ) use_run_fast = false;

    double density = 2.0*E / N;
    if(N > 100'000 && density < 4.0) use_run_fast = false; // on large sparse graphs we do NOT want run_fast()

    if(use_run_fast) createKnownSolutionsUsingRunFast();
    else createKnownSolutions(); // original

    if( !Global::checkTle() ) {
        granulateSolution();

        if (recurrence_depth < cnf->max_recursion_depth) {
            Solver solver(*origV, partition, *cnf, recurrence_depth + 1);
            solver.run_recursive();

            partition = solver.best_partition;
            compareToBestSolutionAndUpdate(partition);

            assert(solver.best_result == PaceUtils::evaluateSolution(*origV, partition)); // FIXME:optimization
        } else {
            partition = best_partition;
        }

        /**
         * CAUTION!!
         * If we do not use localSearch to create known partitions, and the result of the coarsened graph is
         * worse than the result found by this solver, then the best result found by this solver will not be subjected
         * to SwpCndCreators in localSearch!!!
         *
         * We should either improve the best result found, not the best result found for coarsened graph, or improve both
         * found results.
         */

//        clog << "TEST 1 !! ";
//        DEBUG(PaceUtils::evaluateSolution(*origV, partition));
//        DEBUG(PaceUtils::evaluateSolution(*origV, best_partition));

        refineAfterCoarsening();

        if(!Global::disable_all_logs) {
            clog << "Creators: (calls,improvements):" << endl;
            for (auto &[s, p] : local_search_creator_calls) {
                clog << s << " --> " << p << endl;
            }
        }
    }
}

NEG *Solver::createNegForState(State *st) {
    if( cnf->use_neg_map_version ) return new NodeEdgeGreedy(*st);

    NEG * neg;
    bool W1 = true;
    for( int i=0; i<N; i++ ){
        if( clg.node_weights[i] != 1 ){
            W1 = false;
            break;
        }
    }

    if( W1 ) neg = new NodeEdgeGreedyW1(*st);
    else neg = new NodeEdgeGreedyNomap(*st);

    return neg;
}


void Solver::run_fast() {
    if( Global::checkTle() ) return;
//    if( recurrence_depth == max_rec_depth ) return; // original position here

    if(recurrence_depth == 0 && cnf->use_kernelization){
        CEKernelizer kern( *origV );

        if(cnf->use_only_fast_exact_kernelization){ // disable all rules except very fast rules
            for( int i=2; i<kern.MAX_RULES; i++ ) kern.setDisableRule(i,true);
            kern.setDisableRule(15,false);
            kern.setDisableRule(16,false);
            kern.setDisableRule(2,false);
        }

        kern.fullKernelization(false,0);
        if(cnf->use_heuristic_kernelization) kern.improveKernelizationUsingHeuristicRules();

        assert( kern.inCl.size() == N );
        partition = kern.inCl;

        partition = PaceUtils::properlyRemapPartition(partition);
    }

    int MAX_ITERS_BEFORE = 50 / cnf->neg_perm_fraction; // #TEST
    int MAX_ITERS_AFTER = 130 / cnf->neg_perm_fraction; // #TEST

    createClusterGraph();
    known_solutions.clear(); known_clg_partitions.clear();

    auto setNegConfigs = [&]( NEG & neg, double factor = 1.0 ){
        neg.setConfigurations(*cnf);

        neg.use_edge_swaps = true;
        neg.edge_swaps_frequency = 6;

        neg.use_triangle_swaps = false; //neg.triangle_swaps_frequency = 23;
        neg.use_chain2_swaps = false; // neg.chain2_swaps_frequency = 18;
        neg.use_node_interchanging = false; // neg.node_interchanging_frequency = 18;
        neg.use_join_clusters = false;

        neg.max_iterations_to_do = MAX_ITERS_BEFORE;
        neg.move_frequency = cnf->neg_move_frequency;

        neg.allow_perturbations = false; // original
    };

    /*
     * Number of solutions created before granulation
     */
    int REPS = 4; // original 3

    int E = GraphUtils::countEdges( *origV );
    {

        if( E > 500'000 ) REPS = min(REPS,3); // #TEST

        const bool USE_LINEAR_SCALE = false;
        if(USE_LINEAR_SCALE){ // linearly scaling MAX_ITERS_BEFORE
            double A = 500'000;
            double B = 1'300'000;
            double a = 1 / (B - A);
            double b = -a * A;
            if (E >= A && E <= B) {
                double x = a * E + b;
                MAX_ITERS_BEFORE *= (2 - x);
                DEBUG(MAX_ITERS_BEFORE);
            }
        }
    }

    if(recurrence_depth == cnf->solver_max_rec_depth_run_fast) REPS = 1; // #TEST - on last recursion level, getting only one solution - check if to increase recursion depth!

    for(int r=0; r<REPS; r++){ // creating initial solution

        delete st; st = new State(clg, RANDOM_STATE_PERM); // #TEST - creating a new state in each repetition

        if( cnf->solver_run_fast_induce_first_solution_from_lower_levels && r == 0 ){
            if( inducePartitionFromLowerLevelPartition() ){
                st->applyPartition(partition);
                if(!Global::disable_all_logs){
                    clog << "Applying induced partition to state!" << endl;
                    DEBUG(PaceUtils::evaluateState(*st));
                    DEBUG(PaceUtils::evaluateSolution(*origV, lower_level_best_partition_to_induce));
                }
            }
        }

        NEG * neg = createNegForState(st);

        setNegConfigs(*neg);
        neg->improve();

        VI sol = neg->best_partition;
        known_clg_partitions.push_back(sol);

        if(!Global::disable_all_logs) {
            clog << endl;
            logSpacing(recurrence_depth);
            clog << "Quickly found result for clg->V.size(): " << clg.N << "  -  " << neg->best_result << endl;
        }

        sol = PaceUtils::mapClgPartitionToOriginalPartition(clg, sol);
        known_solutions.push_back(sol);

        compareToBestSolutionAndUpdate(sol);

        delete neg;

        if( Global::checkTle() ) break;
    }

    if( recurrence_depth >= cnf->solver_max_rec_depth_run_fast ) return;
    if( Global::checkTle() ) return;

    if(!Global::disable_all_logs) {
        logSpacing(recurrence_depth);
        clog << "Granulating! " << endl;
    }

    partition = best_partition;
    granulateSolution();

    { // running recursively
        Solver solver(*origV, partition, *cnf, recurrence_depth + 1);

        solver.lower_level_best_partition_to_induce = best_partition;

        solver.run_fast();
        compareToBestSolutionAndUpdate(solver.best_partition);
        partition = PaceUtils::mapOriginalPartitionToClgPartition(clg, solver.best_partition);

        if(!Global::disable_all_logs) {
            clog << endl;
            logSpacing(recurrence_depth);
            clog << "Result after recursion: " << best_result << endl;
        }
    }

    if(!Global::checkTle()){ // refinement of found solution
        delete st;
        st = new State(clg, SINGLE_NODES);

        // #CAUTION! Shouldn't be here refinement of best_partition?
        // It would make sense to improve best_partition if it is better than partition returned by recursive call
        // to run_fast()
        st->applyPartition(partition);

        const bool use_first_neg = true; // original true

        if(use_first_neg){ // NEG-based refinement - first NEG, then local search without NEG

            NEG* neg = createNegForState(st);
            neg->setConfigurations(*cnf);

//            neg->chain2_swaps_frequency = 35; neg->join_clusters_frequency = 45; // #TEST
            neg->triangle_swaps_frequency = 30; // original 20
            if( recurrence_depth > 0 ) neg->use_triangle_swaps_to_other_clusters = false;

            neg->edge_swaps_frequency = 11; // original value 11
            neg->max_nonnegative_iters = 20; // originally this was not here - this should disallow any perturbations
            neg->node_interchanging_frequency = 10; // originally this was not here

            if( E > 1'000'000 ){ neg->max_nonnegative_iters = 7; /* #TEST #CAUTION - smaller values should be faster*/ }

            if( E > 1'000'000 && (2.0*E / origV->size()) > 50 ) neg->use_triangle_swaps = false;
            else if( E > 500'000 ){ // #TEST #CAUTION #CAUTION2 - should I disable triangle swaps here?
                neg->max_nonnegative_iters = 12;
            }

            neg->max_iterations_to_do = MAX_ITERS_AFTER;
            neg->allow_perturbations = false; // original version

            neg->improve();
            partition = neg->best_partition;
            delete neg;

            auto part = PaceUtils::mapClgPartitionToOriginalPartition(clg, partition);
            compareToBestSolutionAndUpdate(part);

            bool improve_using_ls = true;

            if (improve_using_ls) {
                if (!Global::disable_all_logs){
                    logSpacing(recurrence_depth);
                    clog << endl << "Improving run_fast with local search without NEG" << endl;
                }
                auto old_creators = cnf->swpCndCreatorsToUse;

                auto it = remove(ALL(cnf->swpCndCreatorsToUse), triangle);
                cnf->swpCndCreatorsToUse.resize(it - cnf->swpCndCreatorsToUse.begin());

                it = remove(ALL(cnf->swpCndCreatorsToUse), node);
                cnf->swpCndCreatorsToUse.resize(it - cnf->swpCndCreatorsToUse.begin());

                partition = best_partition;
                auto[part_oV, part_clg] = localSearch();
                compareToBestSolutionAndUpdate(part_oV);

                if (!Global::disable_all_logs) clog << "LS finished" << endl;
                cnf->swpCndCreatorsToUse = old_creators;
            }

            bool improve_again_using_neg_map = false;
            if(improve_again_using_neg_map){
                if(!Global::disable_all_logs){
                    logSpacing(recurrence_depth);
                    clog << "Running NEG_map after LS with chain, node_interchange and joins" << endl;
                }
                NodeEdgeGreedy neg(*st);
                neg.setConfigurations(*cnf);
                neg.max_iterations_to_do = 30;
                neg.use_triangle_swaps = false;
                neg.use_chain2_swaps = true; neg.chain2_swaps_frequency = 13;
                neg.use_node_interchanging = true; neg.node_interchanging_frequency = 11;
                neg.use_edge_swaps = true; neg.edge_swaps_frequency = 9;
                neg.use_join_clusters = true; neg.join_clusters_frequency = 19;
                neg.allow_perturbations = false;
                neg.improve();
                partition = neg.best_partition;
                auto part = PaceUtils::mapClgPartitionToOriginalPartition(clg, partition);
                compareToBestSolutionAndUpdate(part);
            }
        }
        else{ // NEG-based refinement - first local search without NEG, then NEG
            bool improve_using_ls = true;

            if (improve_using_ls) { // #TEST
                if (!Global::disable_all_logs){
                    logSpacing(recurrence_depth);
                    clog << endl << "Improving run_fast with local search without NEG" << endl;
                }
                auto old_creators = cnf->swpCndCreatorsToUse;

                auto it = remove(ALL(cnf->swpCndCreatorsToUse), triangle);
                cnf->swpCndCreatorsToUse.resize(it - cnf->swpCndCreatorsToUse.begin());

                it = remove(ALL(cnf->swpCndCreatorsToUse), node);
                cnf->swpCndCreatorsToUse.resize(it - cnf->swpCndCreatorsToUse.begin());

                partition = best_partition;
                auto[part_oV, part_clg] = localSearch();
                compareToBestSolutionAndUpdate(part_oV);

                if (!Global::disable_all_logs) clog << "LS finished" << endl;
                cnf->swpCndCreatorsToUse = old_creators;
            }

            partition = PaceUtils::mapOriginalPartitionToClgPartition(clg, best_partition);

            delete st;
            st = new State(clg, SINGLE_NODES);
            st->applyPartition(partition);

            clog << "Proceeding to NEG improvement" << endl;

            NEG* neg = createNegForState(st);
            neg->setConfigurations(*cnf);
            neg->triangle_swaps_frequency = 30; // original 20
            neg->edge_swaps_frequency = 11; // original - not changed
            neg->max_nonnegative_iters = 20; // originally this was not here - this should disallow any perturbations
            neg->node_interchanging_frequency = 10; // originally this was not here

            if( E > 1'000'000 ){
                neg->max_nonnegative_iters = 7; // #TEST #CAUTION - smaller values should be faster
                neg->use_triangle_swaps = false;
            }

                // here should be else if(...) statement
            else if( E > 500'000 ){ // #TEST #CAUTION #CAUTION2 - should I disable triangle swaps here?
                neg->max_nonnegative_iters = 12;
            }

            neg->max_iterations_to_do = MAX_ITERS_AFTER;
            neg->allow_perturbations = false;

            neg->improve();
            partition = neg->best_partition;
            delete neg;

            auto part = PaceUtils::mapClgPartitionToOriginalPartition(clg, partition);
            compareToBestSolutionAndUpdate(part);
        }
    }

    partition = best_partition;

    {
        delete st;
        st = new State(clg, SINGLE_NODES);
        auto clg_part = PaceUtils::mapOriginalPartitionToClgPartition(clg, partition);
        st->applyPartition(clg_part);
    }

    if(!Global::disable_all_logs) {
        clog << endl;
        logSpacing(recurrence_depth);
        clog << "Best result after refinement for clg->V.size():  " << clg.N << "  -  " << best_result << endl;
    }
}



VPII Solver::getModifications() {
    set<PII> mods;
    VVI clusters = PaceUtils::partitionToClusters(best_partition);
//    cout << "clusters: " << clusters << endl;

    for( VI & v : clusters ){
        for( int i=0; i<v.size(); i++ ){
            for(int j=i+1; j<v.size(); j++){
                mods.insert( { min(v[i],v[j]), max(v[i],v[j]) } );
            }
        }
    }

    VPII edges = GraphUtils::getGraphEdges(*origV);
    for( auto & p : edges ) if(p.first > p.second) swap(p.first, p.second);

    for( auto e : edges ){
        int a = e.first, b = e.second;
        if(a>b) swap(a,b);

        if( best_partition[a] == best_partition[b] ) mods.erase(PII(a,b));
        else mods.insert(PII(a,b));
    }

//    cout << "mods: " << mods << endl;
//    cout << "best_partition: " << best_partition << endl;

    return VPII(ALL(mods));
}

pair<VI,VI> Solver::largeIteration(int iter_cnt) {
    if(Global::checkTle()) return {};

    const bool debug = (recurrence_depth == 0);

    createClusterGraph();

    delete st;
    st = new State( clg, cnf->state_init_type ); // create current state, initialized with cnf->state_init_type

    if(debug){
        clog << "origV - nodes: " << origV->size() << ", edges: " << GraphUtils::countEdges(*origV) <<
            "   |   clg nodes: " << clg.V.size() << ", edges: " << GraphUtils::countEdges( clg.V ) << endl;
    }

    /**
     * run algorithms for current state [st]
     */
    pair<VI,VI> new_best_parts = localSearch(iter_cnt);

    return new_best_parts;
}

void Solver::granulateSolution() {
    if(Global::checkTle()) return;

    const bool debug = ( recurrence_depth == 0 );

    if(cnf->coarsen_mode & remove_edges){ // removing edges that have both ends always in different clusters.
        VPII edges = GraphUtils::getGraphEdges(V);

        if( debug ){
            clog << endl << "Entering granulateSolution(), known_solutions.size(): " << known_solutions.size() << endl;

            const bool write_solutions = false;
            if(write_solutions) {
                for (auto v : known_solutions) {
                    VVI sol_cl = StandardUtils::partitionToLayers(v);
                    for (VI &w : sol_cl) sort(ALL(w));
                    sort(ALL(sol_cl), [](VI &v1, VI &v2) {
                        if (v1.empty()) return false;
                        if (v2.empty()) return true;
                        return v1[0] < v2[0];
                    });

                    clog << sol_cl << "  --->  " << PaceUtils::evaluateSolution(*origV, v) << "  |  " << v << endl;
                }
            }
        }

        VPII edges_to_remove;
        for (auto &[a, b] : edges) {
            bool always_different_clusters = true;
            for (int i = 0; i < known_solutions.size(); i++) {
                if (known_solutions[i][a] == known_solutions[i][b]) {
                    always_different_clusters = false;
                    break;
                }
            }
            if (always_different_clusters) edges_to_remove.push_back({a, b});
        }

        if( debug ){
            clog << "There are " << edges_to_remove.size() << " edges to remove from V, out of " << edges.size()
                 << " all edges" << endl;
//            clog << "Those edges are: " << edges_to_remove << endl;
        }
        GraphUtils::removeEdges( V, edges_to_remove );
    }

    if( cnf->coarsen_mode & contract_all ){ // updating partition

        VVI inPartition(N);
        for (int i = 0; i < N; i++) inPartition[partition[i]].push_back(i);

        unordered_map<LL,VI> zb;
        for( int i=0; i<N; i++ ){
            LL hash = 0;
            for( int j=0; j<known_solutions.size(); j++ ){
                LL cluster_hash = hashes[ known_solutions[j][i] % N ];
                LL position_hash = hashes2[j%N];
                hash ^= ( cluster_hash + position_hash );
            }
            zb[hash].push_back(i);
        }

        FAU fau(N);
        for( auto & [hash,vec] : zb ){
//            if(debug) DEBUG(vec);
            for( int i=1; i<vec.size(); i++ ){
                fau.Union( vec[i], vec[0] );
            }
        }

        int cnt = 0;
        unordered_map<int,int> mapper;
        for( int i=0; i<N; i++ ){
            int p = fau.Find(i);
            if( mapper.count(p) == 0 ){
                mapper[p] = cnt;
                p = cnt;
                cnt++;
            }else p = mapper[p];

            partition[i] = p;
        }
    }
    else if( cnf->coarsen_mode & contract_matching ){
//        VVI inPartition(N);
//        for (int i = 0; i < N; i++) inPartition[partition[i]].push_back(i);

        unordered_map<LL,VI> zb;
        for( int i=0; i<clg.N; i++ ){
            LL hash = 0;
            for( int j=0; j<known_clg_partitions.size(); j++ ){
                LL cluster_hash = hashes[ known_clg_partitions[j][i] % N ];
                LL position_hash = hashes2[j%N];
                hash ^= ( cluster_hash + position_hash );
            }
            zb[hash].push_back(i);
        }

        VVI clusters;
        for( auto & p : zb ) clusters.push_back(p.second);

        vector<tuple<int,int,int>> edges; // (a,b,w)
        VB was(N,false);

//        const int MIN_CL_SIZE_TO_CONTRACT = 3; // #TEST #CAUTION original value 0

        for( VI & cl : clusters ){
//            if( cl.size() < MIN_CL_SIZE_TO_CONTRACT ) continue;

            for(int d : cl) was[d] = true;

            for(int d : cl) {
                for (auto &[p, w] : clg.V[d]) {
                    if (d < p && was[p]) {
                        edges.emplace_back( d,p,w );
//                        clog << "Adding edge (" << d << "," << p << "," << w << ")" << endl;
                    }
                }
            }

            for(int d : cl) was[d] = false;
        }

        sort(ALL(edges), [&]( auto& e1, auto& e2 ){
           int a1 = get<0>(e1), b1 = get<1>(e1), w1 = get<2>(e1);
           int a2 = get<0>(e2), b2 = get<1>(e2), w2 = get<2>(e2);

           const int mode = 2; // #TEST - granulation mode - joins smallest - originally mode = 1
            if(mode == 0){
                // we prefer to join those nodes that most probably should be together - swap value
                int score1 = clg.node_weights[a1] * clg.node_weights[b1] - (w1<<1);
                int score2 = clg.node_weights[a2] * clg.node_weights[b2] - (w2<<1);
                if(score1 != score2) return score1 < score2;
                else return clg.node_weights[a1] + clg.node_weights[b1] < clg.node_weights[a2] + clg.node_weights[b2];
            }
            else if(mode == 1){
                // we prefer the edge that results in smallest total weight of node
                return clg.node_weights[a1] + clg.node_weights[b1] < clg.node_weights[a2] + clg.node_weights[b2];
            }else if(mode == 2){
                int W1 = clg.node_weights[a1] + clg.node_weights[b1];
                int W2 = clg.node_weights[a2] + clg.node_weights[b2];
                if( W1 != W2 ) return W1 < W2;
                else{
                    int score1 = clg.node_weights[a1] * clg.node_weights[b1] - (w1<<1);
                    int score2 = clg.node_weights[a2] * clg.node_weights[b2] - (w2<<1);
                    return score1 < score2;
                }
            }
        });

//        DEBUG(*st);
//        DEBUG(clg);
//        clog << "edges sorted: ";
//        for( auto e : edges ) clog << get<0>(e) << " " << get<1>(e) << " " << get<2>(e) << endl;
//        exit(1);

        FAU fau(N);
        for( VI & cl : clusters ) {
            for( int d : cl ){
                for( int p : clg.clusterNodes[d] ){
                    fau.Union( p, clg.clusterNodes[d][0] );
                }
            }
        }

        for( auto & [a,b,w] : edges ){
            if(!was[a] && !was[b]){
//                if(debug) clog << "Contracting edge (" << a << "," << b << ")" << endl;
                fau.Union( clg.clusterNodes[a][0], clg.clusterNodes[b][0] );
                was[a] = was[b] = true;
            }
        }

        int cnt = 0;
        unordered_map<int,int> mapper;
        for( int i=0; i<N; i++ ){
            int p = fau.Find(i);
            if( mapper.count(p) == 0 ){
                mapper[p] = cnt;
                p = cnt;
                cnt++;
            }else p = mapper[p];

            partition[i] = p;
        }
    }

//    known_solutions.clear();
//    known_clg_partitions.clear(); // #TEST - do not clear known_solutions in granulateSolution()
}

pair<VI,VI> Solver::createPartitionsForGivenState(State &st) {
    VI part = partition;

    VI part_clg(st.clg->V.size());
    for( auto& cl : st.clusters ){
        for( int i=0; i<cl.g.nodes.size(); i++ ){
            int p = cl.g.nodes[i];

            part_clg[p] = cl.id;

            for( int d : st.clg->clusterNodes[p] ) {
                part[d] = cl.id;
            }
        }
    }

    return {part, part_clg};
}

void Solver::createClusterGraph() {
    clg = ClusterGraph( &V, partition );

    if(!Global::disable_all_logs) {
        clog << endl;
        logSpacing(recurrence_depth);
        clog << "clg.V.size(): " << clg.V.size() << endl;
    }
}

bool Solver::compareToBestSolutionAndUpdate(VI part) {
    const bool debug = !Global::disable_all_logs;

    if( part.empty() ) return false;
    int part_res = PaceUtils::evaluateSolution(*origV, part);

    if( part_res < best_result ){
        if(debug){
            clog << "\tNew, better_found solution was found! old_best = " << best_result
                 << " > " << part_res << " = current_best" << endl;
        }
        best_result = part_res;
        best_partition = part;
        return true;
    }

    return false;
}

bool Solver::smallIteration(int iter_cnt, int nonneg_iter_cnt) {
    if(Global::checkTle()) return false;

    const bool debug = false;
    using Timer = TimeMeasurer;

    /**
     * FIXME: why this happens??!!??!! Brainfucked...
     * CAUTION! We need to do that!! otherwise, for some unknown reason, if we use SwpCndNode, then on graph heur199.gr
     * we grow out of memory!!! No idea why this is so, maybe some kind of very nasty memory fragmentation??
     */
    if( (iter_cnt % 5) == 4 ){
        State old = *st;
        delete st;
        st = new State(old);
    }

    bool improved = false;

    vector<SwapCandidate*> candidates;

    vector<SwpCndNode> res_node, res_node_all;
    vector<SwpCndEdge> res_edge, res_edge_all;
    vector<SwpCndTriangle> res_triangle, res_triangle_all;
    vector<SwpCndEO> res_eo, res_eo_all;
    vector<SwpCndEORepulsion> res_eo_rep, res_eo_rep_all;
    vector<SwpCndEOAttraction> res_eo_attr, res_eo_attr_all;

    // CAUTION! If swp candidates store pointers to those orders from exp_ordres, they will be invalidated when the
    // vector is reallocated. Thus, SwpCndEO cannot be stored
    vector<ExpansionOrder*> exp_orders;

    vector<ExpansionOrderAttraction*> attr_orders;
    vector<ExpansionOrderRepulsion*> rep_orders;

    auto createSwpCndEoForOrders = [&]( vector<ExpansionOrder> & orders,  SwpCndEOCreator & cr ){
        int beg = exp_orders.size();
        for( int j=0; j<orders.size(); j++ ){
            exp_orders.push_back( new ExpansionOrder({-1}, nullptr) ); // adding empty order
            swap(exp_orders.back()->ord, orders[j].ord);
            exp_orders.back()->cl = orders[j].cl;
        }

        vector<SwpCndEO> res = cr.createSwapCandidates(exp_orders, beg, exp_orders.size() );

        bool negative = false;
        for( auto & cnd : res ){
            res_eo.push_back(cnd);
            if(cnd.swpVal() < 0) negative = true;
        }

        return negative; // true if there was a negative swap candidate
    };

    auto clearOrdersAndCandidates = [&](){
        for (auto ptr : attr_orders) delete ptr; for (auto ptr : rep_orders) delete ptr;
        for (auto ptr : exp_orders) delete ptr;
        exp_orders.clear(); attr_orders.clear(); rep_orders.clear();

        vector<SwpCndNode>().swap(res_node);
        vector<SwpCndEdge>().swap(res_edge);

        res_node.clear(); res_edge.clear(); res_triangle.clear();
        res_eo.clear(); res_eo_rep.clear(); res_eo_attr.clear();

        vector<SwapCandidate*>().swap(candidates);
    };

    for( int i=0; i<cnf->swpCndCreatorsToUse.size(); i++ ){

        if( !cnf->keep_all_swap_candidates ){
            // if medium_slow or faster, we clear existing candidates, they will not be used further
            clearOrdersAndCandidates();
        }

        if(Global::checkTle()) return false;

        switch( cnf->swpCndCreatorsToUse[i] ){
            case edge_same_cl:{
                Timer::start("SwpCndEdge_same_cl");
                local_search_creator_calls["SwpCndEdge_same_cl"].first++;

                SwpCndEdgeCreator cr(*st);
                cr.keep_only_nonpositive_candidates = cnf->keep_only_nonpositive_candidates;
                cr.keep_only_best_cluster_to_move_to = cnf->keep_only_best_cluster_to_move_to;
                vector<SwpCndEdge> res;
                for( auto & cl : st->clusters ){
                    auto cnds = cr.create_MoveTo_SwapCandidatesForCluster(cl);
                    for( auto & cnd : cnds ){
                        res.push_back(cnd);
                    }
                }

                bool negative = false;
                for( auto & cnd : res ){
                    res_edge.push_back(cnd);
                    if(cnd.swpVal() < 0) negative = true;
                }

                if(negative) local_search_creator_calls["SwpCndEdge_same_cl"].second++;
                Timer::stop("SwpCndEdge_same_cl");
                break;
            }
            case edge_diff_cl:{
                Timer::start( "SwpCndEdge_diff_cl" );
                local_search_creator_calls["SwpCndEdge_diff_cl"].first++;

                SwpCndEdgeCreator cr(*st);
                cr.keep_only_nonpositive_candidates = cnf->keep_only_nonpositive_candidates;
                cr.keep_only_best_cluster_to_move_to = cnf->keep_only_best_cluster_to_move_to;
                bool only_common_neighbors = cnf->use_only_common_neighbors_in_swp_cnd_edge;
                vector<SwpCndEdge> res = cr.create_MoveTo_SwapCandidates_DifferentClusters( only_common_neighbors );

                bool negative = false;
                for( auto & cnd : res ){
                    res_edge.push_back(cnd);
                    if(cnd.swpVal() < 0) negative = true;
                }

                if(negative) local_search_creator_calls["SwpCndEdge_diff_cl"].second++;
                Timer::stop( "SwpCndEdge_diff_cl" );
                break;
            }
            case edge_all:{
                Timer::start( "SwpCndEdge_all" );
                local_search_creator_calls["SwpCndEdge_all"].first++;

                SwpCndEdgeCreator cr(*st);
                cr.keep_only_nonpositive_candidates = cnf->keep_only_nonpositive_candidates;
                cr.keep_only_best_cluster_to_move_to = cnf->keep_only_best_cluster_to_move_to;
                bool only_common_neighbors = cnf->use_only_common_neighbors_in_swp_cnd_edge;
                vector<SwpCndEdge> res = cr.create_MoveTo_SwapCandidates_DifferentClusters( only_common_neighbors );
                for( auto & cl : st->clusters ){
                    auto cnds = cr.create_MoveTo_SwapCandidatesForCluster(cl);
                    for( auto & cnd : cnds ){
                        res.push_back(cnd);
                    }
                }

                bool negative = false;
                for( auto & cnd : res ){
                    res_edge.push_back(cnd);
                    if(cnd.swpVal() < 0) negative = true;
                }

                if(negative) local_search_creator_calls["SwpCndEdge_all"].second++;
                Timer::stop( "SwpCndEdge_all" );
                break;
            }
            case triangle:{
                Timer::start( "SwpCndTriangle" );
                local_search_creator_calls["SwpCndTriangle"].first++;

                SwpCndTriangleCreator cr(*st);
                cr.keep_only_nonpositive_candidates = cnf->keep_only_nonpositive_candidates;
                cr.keep_only_best_cluster_to_move_to = cnf->keep_only_best_cluster_to_move_to;
                bool only_empty_cluster = cnf->use_only_empty_cluster_in_swp_cnd_triangle;
                vector<SwpCndTriangle> res = cr.create_MoveTo_SwapCandidates( only_empty_cluster );

                bool negative = false;
                for( auto & cnd : res ){
                    res_triangle.push_back(cnd);
                    if(cnd.swpVal() < 0) negative = true;
                }

                if(negative){
                    local_search_creator_calls["SwpCndTriangle"].second++;
                    if(!Global::disable_all_logs) clog << "smallIteration triangle negative swpval" << endl;
                }
                Timer::stop( "SwpCndTriangle" );
                break;
            }
            case node: {
                Timer::start("SwpCndNode");
                local_search_creator_calls["SwpCndNode"].first++;


                const bool USE_NEG = true;

                if (USE_NEG) {
                    if (!Global::disable_all_logs)
                        clog << "Before NEG, result: " << PaceUtils::evaluateState(*st) << endl;

                    NEG *neg = createNegForState(st);
                    neg->setConfigurations(*cnf);

                    int before;
                    if (!Global::CONTEST_MODE) {
                        assert(PaceUtils::evaluateState(*st) == neg->best_result);
                        before = neg->best_result;
                    }

                    neg->move_frequency = cnf->neg_move_frequency;

                    if (!Global::disable_all_logs)
                        clog << "  Before neg.improve(), result: " << neg->best_result << endl;

                    neg->improve();
                    delete st;
                    st = new State(clg, SINGLE_NODES);
                    VVI to_merge = StandardUtils::partitionToLayers(neg->best_partition);
                    st->mergeClusters(to_merge);
                    partition = PaceUtils::mapClgPartitionToOriginalPartition(clg,
                                                                              st->inCl); // #TEST - uncommented seems to work ok

                    if (!Global::disable_all_logs) clog << "  After SwpCndNode, result: " << neg->best_result << endl;

                    if(Global::checkTle()){
                        delete neg;
                        return false;
                    }

                    if (!Global::CONTEST_MODE) {
                        assert(PaceUtils::evaluateState(*st) == neg->best_result);
                        int after = neg->best_result;
                        if (after < before) improved = true; // #TEST #TEST - originally this was not here

                        {
                            VI temp_part = PaceUtils::mapClgPartitionToOriginalPartition(clg, neg->best_partition);
                            assert(neg->best_result == PaceUtils::evaluateSolution(*origV, temp_part));
                        }
                    }

                    delete neg;

                    local_search_creator_calls["SwpCndNode"].second++;
                } else {
                    SwpCndNodeCreator cr(*st);
                    cr.keep_only_nonpositive_candidates = cnf->keep_only_nonpositive_candidates;
                    cr.keep_only_best_cluster_to_move_to = cnf->keep_only_best_cluster_to_move_to;
                    auto res = cr.createSwapCandidates();

                    bool negative = false;
                    for (auto &cnd : res) {
                        res_node.push_back(cnd);
                        if (cnd.swpVal() < 0) negative = true;
                    }
                    if (negative) local_search_creator_calls["SwpCndNode"].second++;
                }

                Timer::stop("SwpCndNode");
                break;
            }
            case exp_ord:{
                Timer::start( "SwpCndEO" );
                local_search_creator_calls["SwpCndEO"].first++;

                SwpCndEOCreator cr(*st);
                cr.keep_only_best_cluster_to_move_to = cnf->keep_only_best_cluster_to_move_to;

                auto orders = cr.createExpansionOrders( cnf->expansionOrderInitialSetProvider );
                if( Global::checkTle() ) return false;

                bool negative = createSwpCndEoForOrders(orders,cr);
                if( Global::checkTle() ) return false;

                if(negative){
                    local_search_creator_calls["SwpCndEO"].second++;
                    if(!Global::disable_all_logs){
                        clog << " --> expOrd negative swpval, time: " << Global::secondsFromStart() << endl;
                    }
                }
                Timer::stop( "SwpCndEO" );
                break;
            }
            case exp_ord_rep:{
                Timer::start( "SwpCndEORep" );
                local_search_creator_calls["SwpCndEORep"].first++;

                ComponentExpansionRepulsion cr(*st);
                cr.min_cluster_size = cnf->min_cluster_size_for_eo_rep;
                cr.keep_only_nonpositive_candidates = cnf->keep_only_nonpositive_candidates; // #TEST - commented to allow more induced orders

                auto [orders, res] = cr.createSwapCandidates();
                rep_orders += orders;
                if( Global::checkTle() ){
                    clearOrdersAndCandidates();
                    return false;
                }

                bool negative = false;
                for( auto & cnd : res ){
                    if( cr.keep_only_nonpositive_candidates ) assert( cnd.swpVal() <= 0 );

                    res_eo_rep.push_back(cnd);
                    if(cnd.swpVal() < 0){
                        negative = true;
//                        clog << "Negative SwpCndEORep for cl.size(): " << cnd.eo->cl->size() << ", k: " << cnd.k << endl;
                    }
                }

                if(!negative){ // inducing orders
                    vector<ExpansionOrder> induced_orders;
                    for( auto * o : orders ){
                        induced_orders += ExpansionOrderRepulsion::induceOrders(*st, *o);
                    }

                    SwpCndEOCreator cr(*st);
                    bool t = createSwpCndEoForOrders(induced_orders, cr);
                    if(t){
//                        if(!Global::disable_all_logs) clog << "   Induced order from SwpCndRep has negative swap value with current negative: "
//                                                           << negative << endl;
                        negative = t;
                    }
                }

                if( Global::checkTle() ){
                    clearOrdersAndCandidates();
                    return false;
                }

                if(negative){
                    local_search_creator_calls["SwpCndEORep"].second++;
                    if(!Global::disable_all_logs){
                        clog << " --> expOrdRep negative swpval, time: " << Global::secondsFromStart() << endl;
                    }
                }
                Timer::stop( "SwpCndEORep" );
                break;
            }
            case exp_ord_attr:{
                Timer::start( "SwpCndEOAttr" );
                local_search_creator_calls["SwpCndEOAttr"].first++;

                ComponentExpansionAttraction cr(*st);
                cr.min_cluster_size = cnf->min_cluster_size_for_eo_attr;
                cr.keep_only_nonpositive_candidates = cnf->keep_only_nonpositive_candidates;  // #TEST - commented to allow more induced orders

                auto [orders, res] = cr.createSwapCandidates();
                attr_orders += orders;
                if( Global::checkTle() ){
                    clearOrdersAndCandidates();
                    return false; // not clearing rep_orders, possible memory leak
                }

                bool negative = false;
                for( auto & cnd : res ){
                    if( cr.keep_only_nonpositive_candidates ) assert( cnd.swpVal() <= 0 );

                    res_eo_attr.push_back(cnd);
                    if(cnd.swpVal() < 0){
                        negative = true;
//                        clog << "Negative SwpCndEOAttr for cl.size(): " << cnd.eo->cl->size() << ", k: " << cnd.k << endl;
                    }
                }

                if(!negative){ // inducing orders
                    vector<ExpansionOrder> induced_orders;
                    for( auto * o : orders ){
                        induced_orders += ExpansionOrderAttraction::induceOrders(*st, *o);
                    }

                    SwpCndEOCreator cr(*st);
                    bool t = createSwpCndEoForOrders(induced_orders, cr);
                    if(t){
                        if(!Global::disable_all_logs){
//                            clog << "   Induced order from SwpCndAttr has negative swap value with current negative: "
//                                 << negative << endl;
                        }
                        negative = t;
                    }
                }

                if( Global::checkTle() ){
                    clearOrdersAndCandidates();
                    return false; // not clearing rep_orders, possible memory leak
                }

                if(negative){
                    local_search_creator_calls["SwpCndEOAttr"].second++;
                    if(!Global::disable_all_logs){
                        clog << " --> expOrdAttr negative swpval, time: " << Global::secondsFromStart() << endl;
                    }
                }
                Timer::stop( "SwpCndEOAttr" );
                break;
            }
            default:{
                clog << "No default swap candidate creator mode" << endl;
            }
        }

        vector<SwapCandidate*>().swap(candidates);

        for( auto & cnd : res_node ) candidates.push_back(&cnd);
        for( auto & cnd : res_edge ) candidates.push_back(&cnd);
        for( auto & cnd : res_triangle ) candidates.push_back(&cnd);
        for( auto & cnd : res_eo ) candidates.push_back(&cnd);
        for( auto & cnd : res_eo_rep ) candidates.push_back(&cnd);
        for( auto & cnd : res_eo_attr ) candidates.push_back(&cnd);

        if(debug){
            DEBUG(*st);

            if( !res_node.empty() ) DEBUG(res_node);
            if( !res_edge.empty() ) DEBUG(res_edge);
            if( !res_triangle.empty() ) DEBUG(res_triangle);
            if( !res_eo.empty() ) DEBUG(res_eo);
            if( !res_eo_attr.empty() ) DEBUG(res_eo_attr);
            if( !res_eo_rep.empty() ) DEBUG(res_eo_rep);

            if(!candidates.empty()){
                clog << "candidates: " << endl;
                for( auto * cnd : candidates ){
                    clog << "[" << cnd->getNodesToSwap() << ", swpval: " << cnd->swpVal() << "]" << endl;
                }
            }
        }


        if( cnf->apply_swap_on_first_negative && !candidates.empty() ){
            bool has_nonpositive = false;
            for( auto * cnd : candidates ){
                if( cnd->swpVal() <= 0 ){ has_nonpositive = true; break; }
            }

            vector<SwapCandidate*> candidates_filtered;
            if( !cnf->keep_only_nonpositive_candidates ) candidates_filtered = candidates;
            else{
                for( auto * cnd : candidates ){
                    if( cnd->swpVal() < 0 ) candidates_filtered.push_back(cnd); // admitting only negative candidates
                }

                if( !candidates_filtered.empty() ){
                    // if there exists a candidate with negative swap value, then we also add all those with 0 swpval
                    // to increase 'perturbations'
                    for( auto * cnd : candidates ){
                        if( cnd->swpVal() == 0 ) candidates_filtered.push_back(cnd); // admitting only negative candidates
                    }
                }
            }

            if(debug && !candidates_filtered.empty()){
                clog << "candidates_filtered: " << endl;
                for( auto * cnd : candidates_filtered ){
                    clog << "[" << cnd->getNodesToSwap() << ", swpval: " << cnd->swpVal() << "]" << endl;
                }
            }

            if(has_nonpositive) {
                if (debug) clog << "Trying to apply changes, there exists a nonpositive swap candidate" << endl;

                /**
                 * ********* CAUTION!!
                 * If we admit application of neutral (swpVal == 0) candidates, then state [st] is updated,
                 * but apply_swap_for_candidates will return false.
                 * We MUST start new smallIteration, because data in stored expansion orders and swap candidates will
                 * most probably be invalidated.
                 *
                 * Solution 1: admit only negative swaps - this may inhibit getting out of local optimums
                 * Solution 2: gather all neutral candidates from all creators, then apply them at the very end and
                 * the update state
                 */
                bool changes = apply_swap_for_candidates(candidates_filtered);
                if (changes) {
                    improved = true;
                    break;
                }
            }
        }
    }

    perturbState(nonneg_iter_cnt);

    { // clearing pointers to avoid memory leak
//        for (auto ptr : exp_orders) delete ptr;
//        for (auto ptr : attr_orders) delete ptr;
//        for (auto ptr : rep_orders) delete ptr;
        clearOrdersAndCandidates();
    }

    return improved;
}

void Solver::perturbState(int nonneg_iter_cnt) {
    return;
}

bool Solver::apply_swap_for_candidates(vector<SwapCandidate *> &candidates) {
    if(candidates.empty()) return false;
    if(Global::checkTle()) return false;

    const bool debug = false;
    bool improved = false;

    std::mt19937_64 drng; // do we need to seed at all? may do not make any shuffle??

    StandardUtils::shuffle(candidates, drng);
    sort(ALL(candidates), []( SwapCandidate* a, SwapCandidate* b ){
        // sorting candidates by their swap value
        return a->swpVal() < b->swpVal(); // old sorting
    });

    if(debug){
        clog << "Candidates sorted: " << endl;
        for( auto * cnd : candidates ){
            clog << "[" << cnd->getNodesToSwap() << ", swpval: " << cnd->swpVal() << "]" << endl;
        }
    }

    if( cnf->swap_application_mode == GREEDY_MAXIMAL_DISJOINT ){
        VB affected_clusters( 2 * st->N,false );
        VPII to_swap;

        int first_empty_for_all_cnds = st->getIdOfEmptyCluster();

        for( auto * cnd : candidates ){
            if(cnd->swpVal() < 0) improved = true;

            VI affected = cnd->getAffectedClusters(*st);
            bool can = true;
            for( int c : affected ){
                if( c < st->getIdOfEmptyCluster() && affected_clusters[c] ){ can = false; break; }
            }

            const bool APPLY_ALL_SWAPS_REGARDLESS_OF_SAME_CLUSTERS = false; // should be false
            if(APPLY_ALL_SWAPS_REGARDLESS_OF_SAME_CLUSTERS) can = true;

            if(can){
                if(debug){
                    clog << "Adding candidate " << cnd->getNodesToSwap() << " with swpval: "
                         << cnd->swpVal() << " to apply swap" << endl;
                }
                for( int c : affected ) affected_clusters[c] = true;

                auto nodes_to_swap = cnd->getNodesToSwap();

                unordered_map<int,int> empty_cl_mapper;
                unordered_map<int,int> marker;
                int first_empty = first_empty_for_all_cnds;
                for( PII & p : nodes_to_swap ){
                    if( p.second < st->getIdOfEmptyCluster() ) marker[p.first] = p.second;
                    else{
                        auto it = empty_cl_mapper.find( p.second );
                        if( it == empty_cl_mapper.end() ){ // cluster p.second is NOT present - creating it
                            empty_cl_mapper[p.second] = first_empty;
                            marker[p.first] = first_empty;
                            first_empty++;
                        }else{ // cluster p.second is PRESENT, with id empty_cl_mapper[p.second]
                            marker[p.first] = empty_cl_mapper[p.second];
                        }
                    }

                    to_swap.push_back( {p.first, marker[p.first]} );
                }

                first_empty_for_all_cnds = first_empty;
            }
        }

        if(debug) clog << "Applying swap: " << to_swap << endl;
        st->applySwap(to_swap);

    }else if( cnf->swap_application_mode == ONLY_BEST_ONE ){
        SwapCandidate * cnd = candidates[0];
        if( cnd->swpVal() <= 0 ){
            if( cnd->swpVal() < 0 ) improved = true;
            auto to_swap = cnd->getNodesToSwap();
            st->applySwap( to_swap );
        }
    }else{
        cerr << "No such option as  cnf->swap_application_mode = " <<  cnf->swap_application_mode << endl;
    }

    if( debug ){
        clog << "After swap, state: " << *st << endl;
    }

    if( candidates.empty() || candidates[0]->swpVal() >= 0 ) return false;
    return improved;
}


pair<VI,VI> Solver::localSearch(int iter_cnt) {
    if(Global::checkTle()) return createPartitionsForGivenState(*st);

    const bool debug = !Global::disable_all_logs; //( recurrence_depth == 0 );


    int iter = 1, iter_nonneg = 0;
    bool changes = true;

    if(debug) clog << "Starting local search, partition result: " << PaceUtils::evaluateSolution(V, partition) << endl;

    {
        delete st;
        st = new State(clg, SINGLE_NODES);
        auto clg_part = PaceUtils::mapOriginalPartitionToClgPartition(clg, partition);
        st->applyPartition(clg_part);
    }

    while( changes || iter_nonneg < cnf->max_nonnegative_iterations ){
        if(debug){
            int cur_res = PaceUtils::evaluateState(*st);
            clog << "\rsmallIteration iteration #" << iter << ", st->clusters.size(): "
                 << st->clusters.size() << ", current_res: " << cur_res << flush;
        }

        int res_before_small_iter;
        if(!Global::disable_all_logs) res_before_small_iter = PaceUtils::evaluateState(*st); // #TEST FIXME:optimize unnecessary

        auto old_neg_do_not_perturb_if_improved = cnf->neg_do_not_perturb_if_improved;
        cnf->neg_do_not_perturb_if_improved = true; // #TEST #CAUTION

        changes = smallIteration(iter, iter_nonneg);

        cnf->neg_do_not_perturb_if_improved = old_neg_do_not_perturb_if_improved;

        if(Global::checkTle()) return createPartitionsForGivenState(*st);

        if(!Global::disable_all_logs){
            int res_after_small_iter = PaceUtils::evaluateState(*st); // #TEST FIXME:optimize unnecessary
            assert(res_after_small_iter <= res_before_small_iter);
        }

        if( changes ) iter_nonneg = 0;
        else iter_nonneg++;

        iter++;
    }

    return createPartitionsForGivenState(*st);
}

void Solver::createKnownSolutions() {
    if(Global::checkTle()) return;

    known_solutions.clear();
    /**
     * modes:
     * 0 - create cnf->granularity_frequency different solutions using NodeEdgeGreedy
     * 1 - create only 1 know solution using localSearch()
     */
    int mode = 0;

    if(mode == 0) {
        if(!Global::disable_all_logs) clog << "Creating known partitions using NEG" << endl;

        for (int r = 0; r < cnf->granularity_frequency; r++) {
            if(Global::checkTle()) return;

            if(!Global::disable_all_logs) clog << "\rCreating known solutions, r: " << r << " / " << cnf->granularity_frequency << ", best_result: "
                                               << best_result << flush;

            const bool ONLY_NEG = cnf->solver_use_only_neg_to_create_known_solutions;
            if(ONLY_NEG) {
                delete st;
//                st = new State(clg, RANDOM_MATCHING); // original
                st = new State(clg, cnf->state_init_type); // original

                NEG* neg = createNegForState(st);
                neg->setConfigurations(*cnf);
                neg->move_frequency = cnf->neg_move_frequency;
                neg->improve();

                VI sol = neg->best_partition;
                known_clg_partitions.push_back(sol);

                if(!Global::disable_all_logs) clog << ", current_result: " << neg->best_result << endl;

                sol = PaceUtils::mapClgPartitionToOriginalPartition(clg, sol);
                known_solutions.push_back(sol);

                compareToBestSolutionAndUpdate(sol);
                delete neg;
            }else{
                delete st;
                st = new State(clg, RANDOM_MATCHING);

                auto [part_oV, part_clg] = localSearch();
                compareToBestSolutionAndUpdate(part_oV);

                if(Global::checkTle()) return;

                known_solutions.push_back(part_oV);
                known_clg_partitions.push_back(part_clg);
            }
        }

        if(Global::checkTle()) return;

        const bool IMPROVE_BEST_RESULT_USING_LOCAL_SEARCH = cnf->solver_improve_best_known_solution_using_local_search;
        if(IMPROVE_BEST_RESULT_USING_LOCAL_SEARCH){
            if(!Global::disable_all_logs) clog << "Calling local search for best result found" << endl;
            partition = best_partition;
            VI clg_part = PaceUtils::mapOriginalPartitionToClgPartition(clg, partition);
            clg_part = PaceUtils::properlyRemapPartition(clg_part);
            st->applyPartition(clg_part);

            auto [part_oV, part_clg] = localSearch();
            if(Global::checkTle()) return;

            compareToBestSolutionAndUpdate(part_oV);

            known_solutions.push_back(part_oV);
            known_clg_partitions.push_back(part_clg);
        }

    }else if(mode == 1){
        if(!Global::disable_all_logs) clog << "Creating single known partition using localSearch" << endl;

        delete st;
        st = new State(clg, RANDOM_MATCHING);

        auto [ part_oV, part_clg ] = localSearch();
        if(Global::checkTle()) return;

        known_clg_partitions.push_back(part_clg);
        known_solutions.push_back(part_oV);
        compareToBestSolutionAndUpdate(part_oV);
    }
}

void Solver::createKnownSolutionsUsingRunFast() {
    if(Global::checkTle()) return;

    known_solutions.clear();

    int REPS = 1;
    int E = GraphUtils::countEdges( *origV );
    if( E < 100'000 ) REPS = 2; // #TEST #CAUTION - parameter change in Solver based on origV sizes
    else if( E < 60'000 ) REPS = 3;
    else if( E < 30'000 ) REPS = 4;

    for(int r=0; r<REPS; r++) {

        { // running run_fast()
            Solver solver(V, partition, *cnf, 0);
            solver.run_fast();

            compareToBestSolutionAndUpdate(solver.best_partition);
            partition = PaceUtils::mapOriginalPartitionToClgPartition(clg, solver.best_partition);

            {// originally uncommented
                known_solutions += solver.known_solutions; // appending known solution
                known_clg_partitions += solver.known_clg_partitions;
            }
        }

        if (!cnf->solver_improve_best_known_solution_using_local_search) { // refinement of a solution found
            delete st;
            st = new State(clg, SINGLE_NODES);
            st->applyPartition(partition);

            NEG* neg = createNegForState(st);
            neg->setConfigurations(*cnf);
            neg->triangle_swaps_frequency = 30;
            neg->edge_swaps_frequency = 11;
            neg->max_nonnegative_iters = 20;
            neg->node_interchanging_frequency = 10;
            neg->chain2_swaps_frequency = 17; // #TEST

            neg->improve();
            partition = neg->best_partition;
            delete neg;

            auto part_oV = PaceUtils::mapClgPartitionToOriginalPartition(clg, partition);
            compareToBestSolutionAndUpdate(part_oV);

            known_solutions.push_back(part_oV);
            known_clg_partitions.push_back(partition);
        } else {
            partition = best_partition;
            VI clg_part = PaceUtils::mapOriginalPartitionToClgPartition(clg, partition);
            clg_part = PaceUtils::properlyRemapPartition(clg_part);
            delete st;
            st = new State(clg, SINGLE_NODES);
            st->applyPartition(clg_part);

            auto[part_oV, part_clg] = localSearch();
            compareToBestSolutionAndUpdate(part_oV);

            if (Global::checkTle()) return;

            known_solutions.push_back(part_oV);
            known_clg_partitions.push_back(part_clg);
        }

    }
}


void Solver::refineAfterCoarsening() {
    if(!Global::disable_all_logs) clog << endl << endl << "Refinement of solution of  clg nodes: " << clg.V.size() << ", edges: "
                                       << GraphUtils::countEdges( clg.V ) << endl;

//    if(!Global::checkTle()){ // //FIXME:optimize
//        int before = PaceUtils::evaluateSolution(*origV, partition);
//        VI clg_part = PaceUtils::mapOriginalPartitionToClgPartition( clg, partition );
//        VI part = PaceUtils::mapClgPartitionToOriginalPartition(clg, clg_part);
//        int after = PaceUtils::evaluateSolution( V, part );
////        DEBUG2(before,after);
//        assert(before == after);
//    }

    VI clg_part = PaceUtils::mapOriginalPartitionToClgPartition( clg, partition );

//    int before;
//    if(!Global::checkTle()){
//        before = PaceUtils::evaluateSolution(*origV, partition);//FIXME:optimize
//        if(!Global::disable_all_logs) clog << "Before refinement, result: " << before << endl;
//    }



    delete st;
    st = new State( clg, SINGLE_NODES );
    VVI to_merge = PaceUtils::partitionToClusters(clg_part);
    st->mergeClusters(to_merge);

//    int after;
//    if(!Global::checkTle()){
//        after = PaceUtils::evaluateState(*st);//FIXME:optimize
//        assert( before == after );
//    }
//
//
//    if(!Global::checkTle()){ // FIXME:optimize
//        VI clg_part = st->inCl;
//        VI part = PaceUtils::mapClgPartitionToOriginalPartition(clg, clg_part);
//        DEBUG(PaceUtils::evaluateSolution( V, part ));
//    }

    auto [part_oV, part_clg] = localSearch();

    compareToBestSolutionAndUpdate( part_oV );

    if(!Global::disable_all_logs) clog << "After refinement, result: " << best_result << endl;
}

void Solver::logSpacing(int depth) {
    for(int i=0; i<depth; i++) clog << "\t";
}

bool Solver::inducePartitionFromLowerLevelPartition() {
    if(lower_level_best_partition_to_induce.empty()) return false;
    partition = PaceUtils::mapOriginalPartitionToClgPartition(clg,lower_level_best_partition_to_induce);
    return true;
}








