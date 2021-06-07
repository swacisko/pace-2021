//
// Created by sylwester on 4/19/21.
//

#include <CONTESTS/PACE21/heur/StateImprovers/NodeEdgeGreedy.h>
#include <CONTESTS/PACE21/heur/StateImprovers/NodeEdgeGreedyNomap.h>
#include <CONTESTS/PACE21/heur/StateImprovers/SparseGraphTrimmer.h>
#include <CONTESTS/PACE21/test_graphs.h>
#include <CONTESTS/PACE21/heur/StateImprovers/NodeEdgeGreedyW1.h>
#include "CONTESTS/PACE21/main_CE.h"

void kernelizationCompare(){

//    ofstream clog_buf("input_kern_stats.txt", ios::out | ios::app);
//    clog.rdbuf(clog_buf.rdbuf());

    int A = 195;
    int B = 197;

    int cnt = 0;

    VI ruleAppliedCnt(1000,0);

    for (int i = A; i <= B; i += 2) {
        string id = "";
        if (i < 100) id += "0";
        if (i < 10) id += "0";
        id += to_string(i);

//        string inst = "exact";
        string inst = "heur";

        string file = inst + "/" + inst + id + ".gr";

        clog << "Processing file: " << file << endl;
        ifstream istr(file);

//        VVI V = GraphReader::readGraphDIMACSWunweighed(istr);
        VVI V = GraphReader::readGraphDIMACSWunweighed(cin);

        GraphUtils::writeBasicGraphStatistics(V);

        const int ITERS = 0;

        for(int iter=0; iter <= ITERS; iter++ ){

            CEKernelizer kern(V,cnt++);

            /* { // testing heur rule 6 and 7
                 for(int t=0; t<kern.MAX_RULES; t++ ) kern.setDisableRule(t,true);
                 for(int t=0; t<kern.MAX_HEUR_RULES; t++ ) kern.setDisableHeurRule(t,true);

                 kern.setDisableRule(16,false); // correct almost clique rule
                 kern.setDisableHeurRule(6,false);
                 kern.setDisableHeurRule(7,false);

                 kern.fullKernelization(false, 0);
                 kern.improveKernelizationUsingHeuristicRules();
                 continue;
             }*/


            kern.fullKernelization(false, 0);

            { // some possible modifications
//                kern.use_heuristic_rules_separately = true;
//            kern.USE_HEUR_RULE_3_WITHOUT_STRICT_SIZE_CONDITION = true;
            }

            kern.improveKernelizationUsingHeuristicRules();

            ENDL(3);

            { // gathering info about kernelizations used and times used
                VPII apps;
                for( int i=1; i<ruleAppliedCnt.size(); i++ ){
                    if( i < kern.ruleAppliedCnt.size() && kern.ruleAppliedCnt[i] > 0 )
                        apps.push_back( {i, kern.ruleAppliedCnt[i]} );
                }

                for( auto [a,b] : apps ) ruleAppliedCnt[a] += b;
            }

        }

        ENDL(7);
    }

    {
        for( int i=0; i<ruleAppliedCnt.size(); i++ ){
            if( ruleAppliedCnt[i] > 0 ){
                clog << "Rule " << i << " was used " << ruleAppliedCnt[i] << endl;
            }
        }
    }
}

void main_CE(){
    std::ios_base::sync_with_stdio(0);
    std::cin.tie(NULL);

    Global::startAlg();

//        clog.rdbuf(nullptr);
//        cerr.rdbuf(nullptr);

    const bool ADD_SIGTERM_CHECK = true;
    if(!ADD_SIGTERM_CHECK){
        Global::max_runtime_in_seconds = 295;
        clog << "Setting maximal time to " << Global::max_runtime_in_seconds << " seconds" << endl;
    }
    else Global::addSigtermCheck();

    TimeMeasurer::start( "Total time" );

    Global::increaseStack();

//    kernelizationCompare();


    VVI V;

    bool read_from_cin = true;
    if(!read_from_cin){
        int test_case_optilio = 35;
        int test_case = 99;
        if(test_case_optilio != -1) test_case = 2*test_case_optilio-1;

        string case_string = "";

        if(test_case < 10) case_string += "0";
        if(test_case < 100) case_string += "0";
        case_string += to_string(test_case);


        clog << "Test case #" << test_case << " -> optil.io#" << (test_case+1) / 2 << endl;
        ifstream istr( "heur/heur" + case_string + ".gr" );
//        ifstream istr( "exact/exact" + case_string + ".gr" );

        V = GraphReader::readGraphDIMACSWunweighed(istr);
    }else V = GraphReader::readGraphDIMACSWunweighed(cin);


    if(!Global::disable_all_logs) GraphUtils::writeBasicGraphStatistics(V);

    Config cnf;
//    cnf.swpCndCreatorsToUse = { exp_ord_attr};
//    cnf.swpCndCreatorsToUse = {node, triangle, exp_ord, exp_ord_attr, exp_ord_rep };
//    cnf.swpCndCreatorsToUse = {node, exp_ord, exp_ord_attr, exp_ord_rep, triangle }; // good order
//    cnf.swpCndCreatorsToUse = {node, exp_ord, exp_ord_rep, exp_ord_attr /*, triangle*/ }; // #TEST
    cnf.swpCndCreatorsToUse = {exp_ord, exp_ord_rep, exp_ord_attr, node /*triangle*/ }; // #TEST - original order
//    cnf.swpCndCreatorsToUse = {node, triangle, exp_ord_rep };

    /**
     * #CAUTION!
     * edge_all seems to cause some bug in smallIteration() - result of state after small iteration is larger than
     * before the iteration
     */

    cnf.granularity_frequency = 10; // original 10
    cnf.max_recursion_depth = 10;
    cnf.keep_only_best_cluster_to_move_to = true;
    cnf.setSpeedMode(medium);
    cnf.state_init_type = RANDOM_MATCHING;
    cnf.coarsen_mode = contract_matching;

    cnf.use_kernelization = false; // #TEST - first main iteration originally without kernelization
    cnf.use_heuristic_kernelization = false; // do not use heuristic in first iterations
    cnf.use_only_fast_exact_kernelization = true; // originally true

    cnf.solver_use_only_neg_to_create_known_solutions = true; //this will be used only if known solutions are created without  run_fast()
    cnf.solver_improve_best_known_solution_using_local_search = true;

    cnf.solver_run_fast_induce_first_solution_from_lower_levels = false;
    cnf.solver_max_rec_depth_run_fast = 3;

    cnf.neg_perm_fraction = 1.0;
    cnf.neg_do_not_perturb_if_improved = true;



    {
        cnf.neg_move_frequency = 1;
        cnf.neg_max_perturb = 10; // #TEST original 10
        cnf.neg_max_nonneg_iters = 10; // original 10
        cnf.neg_use_queue_propagation = false; // original false
        cnf.neg_use_node_interchange = true; // original true
        cnf.neg_use_join_clusters = true; // original true
        cnf.neg_use_two_node_swaps = false; // original false - there is no need to move two nonadjacent nodes FROM a given cluster

        cnf.neg_use_edge_repulsion = true; // using edge repulsion

        cnf.neg_max_best_cl_size_triangle_swaps = 2; // #TEST original 2
        cnf.neg_use_triangle_swaps_to_other_clusters = false;
    }

    { // setting granularity and maximal recursion depth
        int N = V.size();
        int E = GraphUtils::countEdges(V);
        int T = N+E;
        long long T2 = 1ll * N * E;
        double avg_deg = 2.0*E / N;

        if( N > 100'000 ){
            cnf.granularity_frequency = 3;
            cnf.max_recursion_depth = 3;
        }else if( N > 30'000 ){
            cnf.granularity_frequency = 5;
            cnf.max_recursion_depth = 4;
        }else if( N > 10'000 ){
            cnf.granularity_frequency = 5;
            cnf.max_recursion_depth = 5;
        }else if( N > 1'000 && E > 200'000 ){
            cnf.granularity_frequency = 5;
            cnf.max_recursion_depth = 5;
        }else if( N > 1'000 && E < 200'000 ){
            cnf.granularity_frequency = 10;
            cnf.max_recursion_depth = 5;
        }else if( E < 100'000 ){
            cnf.granularity_frequency = 15;
            cnf.max_recursion_depth = 5;
        }

        if( E > 50'000 ){ // #FIXME: this should be in Solver constructor!
            if( avg_deg > 50 ) cnf.neg_use_triangle_swaps = false;
            if( avg_deg > 150 ){
//                cnf.neg_use_edge_swaps = false;
//                cnf.neg_use_node_interchange = false;
            }
        }

//        if( E > 300'000 && avg_deg > 20){ // #TEST - removing triangle check from localSearch if graph is 'dense'
        if( E > 50'000 ){ // #TEST - removing triangle check from localSearch if graph is 'dense'
            auto it = remove(ALL(cnf.swpCndCreatorsToUse), triangle);
            cnf.swpCndCreatorsToUse.resize( it - cnf.swpCndCreatorsToUse.begin() );
        }


        if(avg_deg < 4){ // parameters for graphs with very small average degree
            cnf.neg_use_join_clusters = false;  // original false
            cnf.neg_use_node_interchange = false; cnf.neg_node_interchanging_frequency = 30; // original false
            cnf.neg_use_chain2_swaps = false; cnf.neg_chain2_swaps_frequency = 40; // original false

            cnf.neg_move_frequency = 2;

            cnf.neg_use_edge_swaps = true;
            cnf.neg_edge_swaps_frequency = 15;

            cnf.neg_use_triangle_swaps = true;
            cnf.neg_triangle_swaps_frequency = 70; // originally 70

            cnf.neg_use_queue_propagation = true;
        }
    }

    {
        if(!Global::disable_all_logs) clog << "Checking if NEG_nomap and NEG:W1 can be used" << endl;
        VI part(V.size(),0); iota(ALL(part),0);
        ClusterGraph clg(&V,part);
        State st(clg, RANDOM_STATE_PERM);
        NodeEdgeGreedyNomap nomap(st);
        nomap.use_edge_swaps = nomap.use_triangle_swaps = nomap.use_node_interchanging = nomap.use_join_clusters =
        nomap.use_chain2_swaps = false;
        nomap.max_iterations_to_do = 3;
        nomap.improve();
        int clusters_cnt = nomap.countNonemptyClusters();

        /**
         * If there are few clusters compared to the number of nodes (that is graph is e.g. very dense or locally
         * very dense), then we use 'map' version of NEG, as it is much faster on that instances.
         */
        const double THR = 3.5;
        if( clusters_cnt < clg.V.size() / THR ) cnf.use_neg_map_version = true;
    }



//    cnf.use_neg_map_version = true; // #TEST


    VI init_part(V.size());
    iota(ALL(init_part),0);

    if(!Global::disable_all_logs){
        DEBUG(cnf.granularity_frequency);
        DEBUG(cnf.max_recursion_depth);

        DEBUG(cnf.neg_move_frequency);
        DEBUG(cnf.neg_use_edge_swaps);
        DEBUG(cnf.neg_use_triangle_swaps);
        DEBUG(cnf.neg_use_node_interchange);
        DEBUG(cnf.neg_use_queue_propagation);
        DEBUG(cnf.neg_use_join_clusters);
        DEBUG(cnf.neg_use_chain2_swaps);

        DEBUG(cnf.use_neg_map_version);
    }

    const bool TEST = false;
    if(TEST){

        {
            ClusterGraph clg(&V, init_part);

            for (int i = 0; i <= 0; i++) {
                State st(clg, RANDOM_STATE_PERM);
                DEBUG(st.clusters.size());

                NEG* neg;
//                neg = new NodeEdgeGreedy(st);
                neg = new NodeEdgeGreedyNomap(st);
//                neg = new NodeEdgeGreedyW1(st);

                neg->setConfigurations(cnf);

                neg->do_not_perturb_if_improved = false;
                neg->max_iterations_to_do = 200 * ( 1.0 / cnf.neg_perm_fraction );


                neg->move_frequency = 2;

                neg->use_node_interchanging = true;     neg->node_interchanging_frequency = 20;
                neg->use_join_clusters = true;          neg->join_clusters_frequency = 20;
                neg->use_chain2_swaps = true;           neg->chain2_swaps_frequency = 20;
                neg->use_queue_propagation = false;

                neg->use_triangle_swaps = true;         neg->triangle_swaps_frequency = 30;
                neg->use_edge_swaps = true;             neg->edge_swaps_frequency = 5;
                neg->use_node_swaps = true; // #TEST


                neg->improve();
                DEBUG(neg->best_result);
                ENDL(3);

                delete neg;
            }
        }

    }else {
        double E = GraphUtils::countEdges(V);
        double avg_deg = 2.0 * E / V.size();

        bool use_run_fast = true;
        if(avg_deg < 4) use_run_fast = false;

        VPII best_mods; int best_result = 1e9;

        bool switcher = false;

        int cnt = 0;
        while( !Global::checkTle() ) {
            if( (cnt % 3) == 2 ) cnf.neg_use_triangle_swaps_to_other_clusters = true;
            else cnf.neg_use_triangle_swaps_to_other_clusters = false;

            Solver solver(V, init_part, cnf);

            if(use_run_fast){
                int old_cnf_use_only_fast_exact_kernelization = cnf.use_only_fast_exact_kernelization;
                if( cnf.use_kernelization && E < 50'000 ){
                    if(switcher) cnf.use_only_fast_exact_kernelization = false;
                    else cnf.use_only_fast_exact_kernelization = true;
                    switcher = !switcher;
                }

                solver.run_fast(); // #TEST
                auto [ part_oV, part_clg ] = solver.localSearch();
                solver.compareToBestSolutionAndUpdate(part_oV);

                cnf.use_kernelization = !cnf.use_kernelization; // changing to use / not to use kernelization in next iteration
                //#TEST - do not use kernelization

                cnf.use_only_fast_exact_kernelization = old_cnf_use_only_fast_exact_kernelization;
            }
            else{
//                solver.run_recursive(); // original
                ClusterGraph clg(&V,init_part);
                State st(clg, RANDOM_MATCHING);
                NEG* neg = new NodeEdgeGreedyW1(st);
                neg->setConfigurations(cnf);

                neg->perturb_mode = 0; // cluster joining instead of splitting
                neg->allow_perturbations = true;
                neg->do_not_perturb_if_improved = false;

                {
                    // #TEST
                    neg->prefer_cluster_mode = 1; // prefer moving to smaller clusters
                }

                neg->improve();

                { // CAUTION - setting values to unused solver
                    solver.best_result = neg->best_result;
                    solver.best_partition = neg->best_partition;
                }

                delete neg;
            }

            if( solver.best_result < best_result ){
                best_result = solver.best_result;
                best_mods = solver.getModifications();
            }

            if(!Global::disable_all_logs){
                clog << "Creators: (calls,improvements):" << endl;
                for( auto & [s,p] : solver.local_search_creator_calls ){
                    clog << s << " --> " << p << endl;
                }
                clog << endl << endl << endl << endl << "********************* NEXT MAIN ITERATION, current best: "
                                               << best_result << endl << endl;
            }/*else{
                cerr << endl << endl << endl << endl << "********************* NEXT MAIN ITERATION, current best: "
                     << best_result << endl << endl;
            }*/


            if(cnt++ == 24) break;
        }

//        bool write_mods = Global::CONTEST_MODE;
        bool write_mods = true; // #TEST
        if (write_mods) {
            VPII mods = best_mods;
            for (auto e : mods) cout << e.first+1 << " " << e.second+1 << endl;
        }

        cerr << "Final result: " << best_result << endl;
        cerr << "Total real time: " << Global::secondsFromStart() << endl;
    }

    if(!Global::disable_all_logs) {
        TimeMeasurer::stop("Total time");
        TimeMeasurer::write();
    }

    return;
}
