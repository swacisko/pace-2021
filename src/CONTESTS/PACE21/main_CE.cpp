//
// Created by sylwester on 4/19/21.
//

#include <CONTESTS/PACE21/heur/StateImprovers/NodeEdgeGreedy.h>
#include <CONTESTS/PACE21/heur/StateImprovers/NodeEdgeGreedyNomap.h>
#include <CONTESTS/PACE21/heur/StateImprovers/SparseGraphTrimmer.h>
#include <CONTESTS/PACE21/test_graphs.h>
#include <CONTESTS/PACE21/heur/StateImprovers/NodeEdgeGreedyW1.h>
#include <graphs/components/ConnectedComponents.h>
#include <getopt.h>
#include "CONTESTS/PACE21/main_CE.h"


const string kernelize = "kernelize";
const string lift_solution = "lift-solution";

string kernelization_phase = kernelize;

void initializeParams(int argc, char **argv) {
    string phase = "phase";

    static struct option long_options[] = {
            {phase.c_str(),      optional_argument, 0, 0},
            {0,0,0,0}
    };

    while (optind < argc) {
        int option_index = 0; int c;
        string option, option_name;

        c = getopt_long(argc, argv, "", long_options, &option_index);

        switch (c) {
            case 0: {
                option = string(optarg);
                option_name = string(long_options[option_index].name);

                if (option_name == phase) {
                    string temp = string(optarg);
                    if(temp == lift_solution) kernelization_phase = lift_solution;
                    else if( temp != kernelize ){
                        cerr << "No such option as " << temp << " provided by parameter --" << option_name << endl;
                        cerr << "Proceeding to kernelization phase (by default)" << endl;
                    }
                }
                break;
            }
            default:{
                option = string(optarg);
                option_name = string(long_options[option_index].name);
            }
        }
    }
}

void liftSolution(VVI & V){
    VVI additional_clusters;
    VI mapper;

    auto getEdgeMods = []( VVI & V, VVI G ){
        VPII edgesV = GraphUtils::getGraphEdges(V);
        VPII edgesG = GraphUtils::getGraphEdges(G);
        set<PII> eV(ALL(edgesV)), eG(ALL(edgesG));
        VPII res;
        for( PII p : edgesV ) if( eG.count(p) == 0 ) res.push_back(p);
        for( PII p : edgesG ) if( eV.count(p) == 0 ) res.push_back(p);
        return res;
    };

    auto transformToG = [&]( VVI & newV ){
        VVI res(V.size());
        for( int i=0; i<newV.size(); i++ ){
            for( int d : newV[i] ){
                res[ mapper[i] ].push_back( mapper[d] );
            }
        }
        return res;
    };

    auto liftSolution = [&]( VVI newV ){
//        clog << "Lifting solution! Before, newV: " << endl << newV << endl;
        newV = transformToG(newV);
//        clog << "After, newV: " << endl << newV << endl;
        for( int i=0; i<additional_clusters.size(); i++ ){
            for( int a : additional_clusters[i] ){
                for( int b : additional_clusters[i] ){
                    if(a==b) continue;
                    newV[a].push_back(b);
                }
            }
        }
        return newV;
    };

    auto applyModifications = [&]( VVI & resV, VPII & mods ){
        set<PII> edge_set;
        auto edges = GraphUtils::getGraphEdges(resV);
        for(PII p : edges){
            edge_set.insert(p);
            edge_set.insert({p.second, p.first});
        }

//        DEBUG(edge_set);

        VPII to_add, to_remove;
        for( PII p : mods ){
            if( edge_set.count(p) > 0 ) to_remove.push_back(p);
            else to_add.push_back(p);
        }

//        DEBUG(to_remove);
//        DEBUG(to_add);

        GraphUtils::removeEdges(resV, to_remove);
        for(auto [a,b] : to_add) GraphUtils::addEdge(resV,a,b);
    };

    //***************************************************************

    V = GraphReader::readGraphDIMACSWunweighed(cin);
    VVI resV;
    VPII modifications;

    mapper = VI(V.size(),-1);

    if (kernelization_phase == lift_solution) {
        string s;
        int edges_read = 0, lineNumber = 1;
        int N,M;
        int resN, resM;

        int section = 0;
        int section1_cnt = 0;

        while (!cin.eof()) {
            getline(cin, s);
            if(cin.eof()) break;

            if (s[0] == 'c') {
                stringstream str(s);
                if( s.size() >= string("c mapper").size() && s.substr( 2, 6 ) == "mapper" ) {
                    char e;
                    string temp = "";
                    str >> e >> temp;
                    int a, cnt = 0;

                    str >> cnt;
                    for(int i=0; i<cnt; i++) str >> mapper[i];

//                    if(!Global::disable_all_logs)
//                        DEBUG(mapper);

                }else if( s.size() >= string("c additional_clusters").size() && s.substr( 2,19 ) == "additional_clusters" ){
                    char e; string temp; int a;
                    str >> e >> temp >> a;
                    for(int i=0; i<a; i++){
                        int b,x;
                        cin >> e >> x;
                        additional_clusters.push_back(VI());
                        for( int j=0; j<x; j++ ){
                            cin >> b;
                            additional_clusters.back().push_back(b);
                        }
                    }
                    cin.ignore();

//                    if(!Global::disable_all_logs)
//                        DEBUG(additional_clusters);
                }

            } else if (s[0] == 'p') {
                if(section == 0){
                    stringstream str(s);
                    string nothingBox;
                    str >> nothingBox >> nothingBox >> N >> M;
                    V = VVI(N);
                }else if(section == 1){
                    stringstream str(s);
                    string nothingBox;
                    str >> nothingBox >> nothingBox >> resN >> resM;
                    resV = VVI(resN);
                }

            } else if( s[0] == '#' ){
                section++;
//                if(!Global::disable_all_logs){
//                    clog << "Changing section, now section: " << section << endl;
//                    DEBUG(V);
//                    ENDL(10);
//                    DEBUG(resV);
//                }
            }
            else {
                if(section == 0) {
                    stringstream str(s);
                    int a, b; char e;

                    str >> a >> b;
                    edges_read++;

                    a--; b--;

                    V[a].push_back(b);
                    V[b].push_back(a);
                }else if(section == 1){
                    stringstream str(s);
                    int a, b; char e;

                    if(section1_cnt == 0){
                        str >> a; // reading number 'd' - this is useless in solution lifting
//                        DEBUG(a);
                        section1_cnt++;
                        continue;
                    }

//                    DEBUG(s);
                    str >> a >> b;
                    edges_read++;

//                    DEBUG2(a,b);

                    a--; b--;

//                    clog << "Adding edge " << a << " -> " << b << " to resV" << endl;

                    resV[a].push_back(b);
                    resV[b].push_back(a);
                }else if(section == 2){
//                    cerr << "Reading modifications" << endl;
                    stringstream str(s);
                    int a, b;
                    str >> a >> b;
                    a--; b--;
                    modifications.push_back( {a,b} );
                }
            }
        }

//        if(!Global::disable_all_logs)
//            DEBUG(modifications);


        applyModifications(resV, modifications);
//        cerr << "Graph after modifications: " << resV << endl;

        auto newV = liftSolution(resV);

//        DEBUG(newV);

        auto mods = getEdgeMods(V, newV);
        cout << mods.size() << endl;
        for(PII p : mods) cout << p.first+1 << " " << p.second+1 << endl;

        cerr << "mods.size(): " << mods.size() << endl;
    }
}

void main_CE(int argc, char **argv){
    std::ios_base::sync_with_stdio(0);
    std::cin.tie(NULL);

//    {
//        auto edges = CE_test_graphs::cluster_graph_test2_edges;
//        auto N = -1;
//        for(PII p : edges) N = max(N, max(p.first+1, p.second+1) );
//        cout << "p cep " << N << " " << edges.size() << endl;
//        for(PII p : edges) cout << p.first+1 << " " << p.second+1 << endl;
//        cerr << "Test written to file" << endl;
//        exit(2);
//    }

    initializeParams(argc, argv);

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
    }else {

        if( kernelization_phase == lift_solution ){
            liftSolution(V);
            return;
        }else{
            V = GraphReader::readGraphDIMACSWunweighed(cin);
        }
    }


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

                neg->use_node_interchanging = false;     neg->node_interchanging_frequency = 10;
                neg->use_join_clusters = false;          neg->join_clusters_frequency = 30;
                neg->use_chain2_swaps = false;           neg->chain2_swaps_frequency = 30;
                neg->use_queue_propagation = false;

                neg->use_triangle_swaps = true;         neg->triangle_swaps_frequency = 7;
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
//        if(avg_deg < 4) use_run_fast = false;

        VVI results;

        bool switcher = false;

        int cnt = 0;
        while( !Global::checkTle() ) {
            Solver solver(V, init_part, cnf);

            if(use_run_fast){
                int old_cnf_use_only_fast_exact_kernelization = cnf.use_only_fast_exact_kernelization;
                if( cnf.use_kernelization && E < 300'000 ){
                    if(switcher) cnf.use_only_fast_exact_kernelization = false;
                    else cnf.use_only_fast_exact_kernelization = true;
                    switcher = !switcher;
                }

                solver.run_fast(); // #TEST
                auto [ part_oV, part_clg ] = solver.localSearch();
                solver.compareToBestSolutionAndUpdate(part_oV);

                cnf.use_kernelization = !cnf.use_kernelization; // changing to use / not to use kernelization in next iteration

                cnf.use_only_fast_exact_kernelization = old_cnf_use_only_fast_exact_kernelization;
            }

            results.push_back( solver.best_partition );

            if(cnt++ == 20) break;
        }


        bool write_mods = true;
        if (write_mods) {
            VPII edges_to_remove;
            VPII edges = GraphUtils::getGraphEdges(V);
            {
                for (auto &[a, b] : edges) {
                    bool always_different_clusters = true;
                    for (int i = 0; i < results.size(); i++) {
                        if (results[i][a] == results[i][b]) {
                            always_different_clusters = false;
                            break;
                        }
                    }
                    if (always_different_clusters) edges_to_remove.push_back({a, b});
                }
            }
            for( auto & [a,b] : edges_to_remove ) if(a>b) swap(a,b);

            VPII edges_to_insert;
            {
                Solver solver(V, init_part, cnf);
                solver.known_solutions = results;
                cnf.coarsen_mode = contract_all;
                solver.granulateSolution();
                VI part = solver.partition;
                set<PII> isEdge(ALL(edges));

                VVI clusters = PaceUtils::partitionToClusters(part);
                for( int i=0; i<clusters.size(); i++ ){
                    for( int j=0; j<clusters[i].size(); j++ ){
                        int a = clusters[i][j];
                        for( int k=j+1; k<clusters[i].size(); k++ ){
                            int b = clusters[i][k];
                            if( isEdge.count( {a,b} ) == 0 ) edges_to_insert.push_back( {a,b} );
                        }
                    }
                }
            }
            for( auto & [a,b] : edges_to_insert ) if(a>b) swap(a,b);


            cerr << "There are " << edges_to_remove.size() << " edges removed and " << edges_to_insert.size()
                 << " edges added" << endl;

            VVI resV = V;
            GraphUtils::removeEdges(resV, edges_to_remove);
            for( auto & [a,b] : edges_to_insert) GraphUtils::addEdge(resV,a,b,false);


            VI to_induce;
            VVI components = ConnectedComponents::getConnectedComponents(resV);
            VVI additional_clusters;

            for( VI & comp : components ){
                bool isCluster = true;
                for( int d : comp ) if( resV[d].size() != (comp.size()-1) ) isCluster = false;
                if(isCluster){
                    additional_clusters.push_back(comp);
                }
                else to_induce += comp;
            }

            InducedGraph indResV = GraphInducer::induce(resV,to_induce);

            /**
             * returns a vector of pairs denoting edge modifications needed to transform V into G (or vice versa)
             */
            auto getEdgeMods = []( VVI & V, VVI G ){
                VPII edgesV = GraphUtils::getGraphEdges(V);
                VPII edgesG = GraphUtils::getGraphEdges(G);
                set<PII> eV(ALL(edgesV)), eG(ALL(edgesG));
                VPII res;
                for( PII p : edgesV ) if( eG.count(p) == 0 ) res.push_back(p);
                for( PII p : edgesG ) if( eV.count(p) == 0 ) res.push_back(p);
                return res;
            };


            bool write_heur_rule_results = false;

            bool write_auxiliary = (kernelization_phase == kernelize);

            VI mapper = indResV.nodes;

            if(write_auxiliary){
                cout << "c mapper " << indResV.nodes.size(); // writing mapper to remap ids.
                for(int i=0; i<mapper.size(); i++) cout << " " << mapper[i];
                cout << endl;

                cout << "c additional_clusters " << additional_clusters.size() << endl; // writing additional_clusters, not mapped
                for( auto & cl : additional_clusters ){
                    cout << "c " << cl.size();
                    for( int d : cl ) cout << " " << d;
                    cout << endl;
                }
            }

            resV = indResV.V;

            auto transformToG = [&]( VVI & newV ){
                VVI res(V.size());
                for( int i=0; i<newV.size(); i++ ){
                    for( int d : newV[i] ){
                        res[ mapper[i] ].push_back( mapper[d] );
                    }
                }
                return res;
            };

            auto liftSolution = [&]( VVI newV ){
                newV = transformToG(newV);
                for( int i=0; i<additional_clusters.size(); i++ ){
                    for( int a : additional_clusters[i] ){
                        for( int b : additional_clusters[i] ){
                            if(a==b) continue;
                            newV[a].push_back(b);
                        }
                    }
                }
                return newV;
            };

            int modifications_done;
            {
                VVI tempV = transformToG(resV);
                for (int i = 0; i < additional_clusters.size(); i++) {
                    for (int a : additional_clusters[i]) {
                        for (int b : additional_clusters[i]) {
                            if (a == b) continue;
                            tempV[a].push_back(b);
                        }
                    }
                }
                VPII temp_mods = getEdgeMods(V, tempV);
                modifications_done = temp_mods.size();
            }

            cout << modifications_done << endl;
            int resN = resV.size(), resE = GraphUtils::countEdges(resV);
            cout << "p cep " << resN << " " << resE << endl;
            for( int i=0; i<resN; i++ ){
                for( int p : resV[i] ){
                    if( p > i ) cout << i+1 << " " << p+1 << endl;
                }
            }

            int N = resV.size();

            auto heur1 = [&](){
                VVI newV(resV.size());

//                newV = transformToG(newV);
//                for( int i=0; i<additional_clusters.size(); i++ ){
//                    for( int a : additional_clusters[i] ){
//                        for( int b : additional_clusters[i] ){
//                            if(a==b) continue;
//                            newV[a].push_back(b);
//                        }
//                    }
//                }
                newV = liftSolution(newV);

                auto mods = getEdgeMods(V, newV);
                cout << mods.size() << endl;
                for(PII p : mods) cout << p.first+1 << " " << p.second+1 << endl;
            };

            auto heur2 = [&](){
                VVI newV(resV.size());
                for( int i=0; i<newV.size(); i++ ){
                    for( int k=0; k<newV.size(); k++ ){
                        if(i != k) newV[i].push_back(k);
                    }
                }

//                newV = transformToG(newV);
//                for( int i=0; i<additional_clusters.size(); i++ ){
//                    for( int a : additional_clusters[i] ){
//                        for( int b : additional_clusters[i] ){
//                            if(a==b) continue;
//                            newV[a].push_back(b);
//                        }
//                    }
//                }
                newV = liftSolution(newV);

                auto mods = getEdgeMods(V, newV);
                cout << mods.size() << endl;
                for(PII p : mods) cout << p.first+1 << " " << p.second+1 << endl;
            };

            auto heur3 = [&](){
                VB was(N,false);
                VVI clusters;
                for( int i=0; i<N; i++ ){
                    if( was[i] ) continue;
                    was[i] = true;
                    clusters.push_back({});
                    clusters.back().push_back(i);

                    for( int d : resV[i] ){
                        if( was[d] ) continue;
                        was[d] = true;
                        clusters.back().push_back(d);
                    }
                }

                VVI newV(N);
                for( int i=0; i<clusters.size(); i++ ){
                    for( int a : clusters[i] ){
                        for( int b : clusters[i] ){
                            if(a==b) continue;
                            newV[a].push_back(b);
                        }
                    }
                }

//                newV = transformToG(newV);
//                for( int i=0; i<additional_clusters.size(); i++ ){
//                    for( int a : additional_clusters[i] ){
//                        for( int b : additional_clusters[i] ){
//                            if(a==b) continue;
//                            newV[a].push_back(b);
//                        }
//                    }
//                }
                newV = liftSolution(newV);

                auto mods = getEdgeMods(V, newV);
                cout << mods.size() << endl;
                for(PII p : mods) cout << p.first+1 << " " << p.second+1 << endl;
            };

            auto heur4 = [&](){
                VB was(N,false);
                VVI clusters;
                for( int i=N-1; i>=0; i-- ){
                    if( was[i] ) continue;
                    was[i] = true;
                    clusters.push_back({});
                    clusters.back().push_back(i);

                    for( int d : resV[i] ){
                        if( was[d] ) continue;
                        was[d] = true;
                        clusters.back().push_back(d);
                    }
                }

                VVI newV(N);
                for( int i=0; i<clusters.size(); i++ ){
                    for( int a : clusters[i] ){
                        for( int b : clusters[i] ){
                            if(a==b) continue;
                            newV[a].push_back(b);
                        }
                    }
                }

//                newV = transformToG(newV);
//                for( int i=0; i<additional_clusters.size(); i++ ){
//                    for( int a : additional_clusters[i] ){
//                        for( int b : additional_clusters[i] ){
//                            if(a==b) continue;
//                            newV[a].push_back(b);
//                        }
//                    }
//                }
                newV = liftSolution(newV);

                auto mods = getEdgeMods(V, newV);
                cout << mods.size() << endl;
                for(PII p : mods) cout << p.first+1 << " " << p.second+1 << endl;
            };

            if(write_heur_rule_results){
                heur1();
                heur2();
                heur3();
                heur4();
            }

            cerr << "value of d: " << modifications_done << endl;
        }

        cerr << "Total real time: " << Global::secondsFromStart() << endl;
    }

    if(!Global::disable_all_logs) {
        TimeMeasurer::stop("Total time");
        TimeMeasurer::write();
    }

    return;
}
