//
// Created by sylwester on 3/18/21.
//

#include <CONTESTS/PACE21/test_graphs.h>
#include "CONTESTS/PACE21/kernelization/CEKernelizer.h"
#include "gtest/gtest.h"


TEST(CEKernelizer, critical_cliques) {
    clog << "Testing creating critical cliques" << endl;

    {
        clog << "\ttest: big_test_50_104" << endl;
        VVI V = CE_test_graphs::big_test_50_104;
        CEKernelizer kern(V);

        auto cc = kern.createCriticalCliques();
        sort(ALL(cc), [](auto &c1, auto &c2) { return c1.nodes.size() > c2.nodes.size(); });
        ASSERT_EQ(cc.size(), 50 - 4);
    }

    {
        clog << "\ttest: path P9" << endl;
        VVI V = CE_test_graphs::pathP9;
        CEKernelizer kern(V);

        auto cc = kern.createCriticalCliques();
        sort(ALL(cc), [](auto &c1, auto &c2) { return c1.nodes.size() > c2.nodes.size(); });
        ASSERT_EQ(cc.size(), 9);
    }
}

TEST(CEKernelizer, editing_degrees) {

    clog << "Testing editing degrees, CE_test_graphs::kern_ed" << endl;
    VVI V = CE_test_graphs::kern_ed;

    CEKernelizer kern(V);

    VI K, N1;
    VB inK(V.size(), false), inN1(V.size(), false);

    {
        K = {0, 1, 2};
        VI N1 = kern.getNeighborhood(K, 1);

        for (int d : K) inK[d] = true;
        for (int d : N1) inN1[d] = true;
        ASSERT_EQ(kern.getEditingDegree(3, inK, inN1, N1.size()) , 5);
        ASSERT_EQ(kern.getEditingDegree(4, inK, inN1, N1.size()) , 2);
        for (int d : K) inK[d] = false;
        for (int d : N1) inN1[d] = false;
    }

    {
        K = {7, 8};
        VI N1 = kern.getNeighborhood(K, 1);

        for (int d : K) inK[d] = true;
        for (int d : N1) inN1[d] = true;
        ASSERT_EQ(kern.getEditingDegree(3, inK, inN1, N1.size()) , 5);
        ASSERT_EQ(kern.getEditingDegree(5, inK, inN1, N1.size()) , 2);
        ASSERT_EQ(kern.getEditingDegree(6, inK, inN1, N1.size()) , 1);
        for (int d : K) inK[d] = false;
        for (int d : N1) inN1[d] = false;
    }

    {
        K = {9};
        VI N1 = kern.getNeighborhood(K, 1);

        for (int d : K) inK[d] = true;
        for (int d : N1) inN1[d] = true;
        ASSERT_EQ(kern.getEditingDegree(3, inK, inN1, N1.size()) , 7);
        ASSERT_EQ(kern.getEditingDegree(4, inK, inN1, N1.size()) , 6);
        ASSERT_EQ(kern.getEditingDegree(5, inK, inN1, N1.size()) , 4);
        ASSERT_EQ(kern.getEditingDegree(6, inK, inN1, N1.size()) , 3);
        for (int d : K) inK[d] = false;
        for (int d : N1) inN1[d] = false;
    }

    clog << "Editing degrees passed" << endl;
}

TEST(CEKernelizer, rule1) {
    clog << "\tRule1: kern1" << endl;
    VVI V = CE_test_graphs::kern1;

    CEKernelizer kern(V);
    kern.rule1();

    ASSERT_EQ( kern.inCl[0] , kern.inCl[1] );

    ASSERT_EQ( kern.inCl[2] , kern.inCl[3] );

    ASSERT_EQ( kern.inCl[4] , kern.inCl[5] );

    ASSERT_EQ( kern.inCl[6] , kern.inCl[7] );
    ASSERT_EQ( kern.inCl[7] , kern.inCl[8] );

    ASSERT_EQ( kern.inCl[9] , 9 );
    ASSERT_EQ( kern.inCl[10] , 10 );
    ASSERT_EQ( kern.inCl[11] , 11 );
}

TEST(CEKernelizer, rule2) {
    clog << "Rule2: CE_test_graphs::kern2" << endl;

    {
        VVI V = CE_test_graphs::kern2;
        CEKernelizer kern(V);
        kern.rule2();

        ASSERT_TRUE(kern.inCl[0] == kern.inCl[1] &&
               kern.inCl[0] == kern.inCl[2]
        );
    }

    {
        VVI V = CE_test_graphs::kern2;

        V[3].push_back(5);
        V.push_back(VI(1, 3)); // added node should prevent kernelization rule 2

        CEKernelizer kern(V);
        kern.rule2();

        ASSERT_EQ(kern.inCl[0] , kern.inCl[1]);
        ASSERT_EQ(kern.inCl[0] , kern.inCl[2]);

        ASSERT_EQ(kern.inCl[3] , 3);
        ASSERT_EQ(kern.inCl[4] , 4);
    }


}

TEST(CEKernelizer, rule3) {
    {
        clog << "Rule3: CE_test_graphs::pathP3" << endl;
        VVI V = CE_test_graphs::pathP3;
        CEKernelizer kern(V);
        kern.rule3();

//        for (int i = 0; i < V.size(); i++) ASSERT_EQ(kern.inCl[i] , i);
        ASSERT_EQ( kern.inCl[1], kern.inCl[2] );
        ASSERT_NE( kern.inCl[0], kern.inCl[2] );
    }

    {
        clog << endl << endl << "Rule3: CE_test_graphs::kern3" << endl;
        VVI V = CE_test_graphs::kern3;
        CEKernelizer kern(V);
        kern.rule3();

        for (int i = 0; i < V.size(); i++) {
            if (i == 6 || i == 11) ASSERT_EQ(V[i].size() , kern.V[i].size() + 1);
            else
                ASSERT_EQ(V[i].size() , kern.V[i].size());
        }
    }
}

TEST(CEKernelizer, rule4) {
    {
        clog << endl << endl << "Rule4: CE_test_graphs::kern4positive" << endl;
        VVI V = CE_test_graphs::kern4positive;
        CEKernelizer kern(V);
        kern.rule4();

        for (int i = 0; i < V.size(); i++) assert(kern.inCl[i] == kern.inCl[0]); // just one clique
    }

    {
        clog << endl << endl << "Rule4: CE_test_graphs::kern4negative" << endl;
        VVI V = CE_test_graphs::kern4negative;
        CEKernelizer kern(V);
        kern.rule4();

        for (int i = 0; i < V.size(); i++) {
            ASSERT_TRUE(equal(ALL(V[i]), ALL(kern.V[i]))); // nothing done by the rule
        }
    }

}

TEST(CEKernelizer, rule6_lemma3) {
    {
        clog << "Rule6_lemma3: CE_test_graphs::kern6positive" << endl;
        VVI V = CE_test_graphs::kern6positive;

        CEKernelizer kern(V);
        kern.rule6_lemma3();

        for (int i = 0; i < 9; i++) assert(kern.inCl[i] == kern.inCl[0]);
        ASSERT_EQ(kern.inCl[9], kern.inCl[10]);
        ASSERT_EQ(kern.inCl[11], kern.inCl[12]);

    }

    {
        clog << endl << endl << "Rule6_lemma3: CE_test_graphs::kern6negative" << endl;
        VVI V = CE_test_graphs::kern6negative;

        CEKernelizer kern(V);
        kern.rule6_lemma3();

        for (int i = 0; i < V.size(); i++) {
            ASSERT_TRUE(equal(ALL(V[i]), ALL(kern.V[i]))); // nothing done by the rule
        }
    }
}

TEST(CEKernelizer, rule7) {
    {
        clog << "Rule7, CE_test_graphs::kern7positive" << endl;
        VVI V = CE_test_graphs::kern7positive;
        CEKernelizer kern(V);
        kern.inV[7] = false; // node 7 is missing
        kern.rule7();

        for (int i = 1; i < 7; i++) ASSERT_EQ(kern.inCl[i], kern.inCl[i - 1]); // clique {0,1,2,3,4,5,6}
    }

    {
        clog << endl << endl << "Rule7, CE_test_graphs::kern7negative" << endl;
        VVI V = CE_test_graphs::kern7negative;
        CEKernelizer kern(V);
        kern.inV[7] = false; // node 7 is missing
        kern.rule7();

        for (int i = 0; i < V.size(); i++) {
            ASSERT_TRUE(equal(ALL(V[i]), ALL(kern.V[i]))); // nothing done by the rule
        }
    }

    {
        clog << endl << endl << "Rule7, CE_test_graphs::kern8positive" << endl;
        VVI V = CE_test_graphs::kern8positive;
        CEKernelizer kern(V);
        kern.rule7();

        for (int i = 0; i < V.size(); i++) {
            ASSERT_TRUE(equal(ALL(V[i]), ALL(kern.V[i]))); // nothing done by the rule
        }
    }

}

TEST(CEKernelizer, rule8) {

    {
        clog << "Rule8, CE_test_graphs::kern8positive" << endl;
        VVI V = CE_test_graphs::kern8positive;
        CEKernelizer kern(V);
        kern.rule8();

        for (int i = 0; i < V.size(); i++) {
            if (i == 3 || i == 5) ASSERT_EQ(kern.V[i].size(), V[i].size() + 1);
            else
                ASSERT_TRUE(equal(ALL(V[i]), ALL(kern.V[i])));
        }
    }

    {
        clog << endl << endl << "Rule8, CE_test_graphs::kern8positive2" << endl;
        VVI V = CE_test_graphs::kern8positive2;
        CEKernelizer kern(V);
        kern.inV[9] = false;
        kern.rule8();

        ASSERT_EQ(find(ALL(kern.V[2]), 8), kern.V[2].end());
        ASSERT_NE(find(ALL(kern.V[4]), 7), kern.V[4].end());
        ASSERT_NE(find(ALL(kern.V[7]), 4), kern.V[7].end());
        for (int d : {0, 1, 11, 12, 3, 5, 6, 10}) ASSERT_TRUE(equal(ALL(V[d]), ALL(kern.V[d])));
    }

    {
        clog << endl << endl << "Rule8, CE_test_graphs::kern8negative" << endl;
        VVI V = CE_test_graphs::kern8negative;
        CEKernelizer kern(V);
        kern.inV[9] = false;
        kern.rule8();

        for (int i = 0; i < V.size(); i++) {
            ASSERT_TRUE(equal(ALL(V[i]), ALL(kern.V[i]))); // nothing done by the rule
        }
    }

}

TEST(CEKernelizer, rule9) {
    {
        clog << "Rule9, CE_test_graphs::kern9positive" << endl;
        VVI V = CE_test_graphs::kern9positive;
        CEKernelizer kern(V);
        kern.rule9();

        for (int i = 0; i < 5; i++) ASSERT_EQ(kern.inCl[i], kern.inCl[0]);

        kern.rule9(); // applying again - this will mark {3,5} as the same cluster

        for (int i = 0; i < V.size(); i++) ASSERT_EQ(kern.inCl[i], kern.inCl[0]);

        ASSERT_TRUE(kern.V[0].empty());
        ASSERT_TRUE(kern.V[1].empty());

        int cnt = 0;
        for (int i = 2; i <= 4; i++) if (kern.V[i].empty()) cnt++;
        ASSERT_EQ(cnt, 2);
    }

    {
        clog << endl << endl << "Rule9, CE_test_graphs::kern9negative" << endl;
        VVI V = CE_test_graphs::kern9negative;
        CEKernelizer kern(V);
        kern.rule9();

        for (int i = 0; i < V.size(); i++) {
            ASSERT_TRUE(equal(ALL(V[i]), ALL(kern.V[i]))); // nothing done by the rule
        }
    }
}

TEST(CEKernelizer, rule15) {
    {
        clog << "Rule15: CE_test_graphs::kern15" << endl;
        VVI V = CE_test_graphs::kern15;
        CEKernelizer kern(V);
        bool changes_done = kern.rule15();
        ASSERT_TRUE(changes_done);

        ASSERT_EQ( kern.inCl[2], kern.inCl[1] );
        ASSERT_EQ( kern.inCl[3], kern.inCl[6] );
        ASSERT_EQ( kern.inCl[4], kern.inCl[7] );
    }
}


TEST(CEKernelizer, rule_heur_1) {
    clog << "Rule_heur_1: CE_test_graphs::kern_heur_1" << endl;
    VVI V = CE_test_graphs::kern_heur_1;
    CEKernelizer kern(V);

    kern.rule_heur_1();

    ASSERT_TRUE(kern.V[11].empty());
    ASSERT_EQ(kern.inCl[9] , kern.inCl[10]);
    ASSERT_TRUE(kern.V[5].empty());
    ASSERT_EQ(kern.inCl[3] , kern.inCl[2]);
    ASSERT_EQ(kern.inCl[6] , kern.inCl[7]);
    ASSERT_EQ(kern.inCl[6] , kern.inCl[8]);
    ASSERT_EQ(kern.inCl[0] , kern.inCl[1]);
    ASSERT_NE(kern.inCl[4] , kern.inCl[3]);
    ASSERT_NE(kern.inCl[4] , kern.inCl[2]);

    kern.rule_heur_1(); // running again

    ASSERT_TRUE(kern.inCl[4] == kern.inCl[3]);
    ASSERT_TRUE(kern.inCl[4] == kern.inCl[2]);
}

TEST(CEKernelizer, rule_heur_2) {
    {
        clog << "Rule_heur_2: CE_test_graphs::kern_heur_2" << endl;
        VVI V = CE_test_graphs::kern_heur_2;
        CEKernelizer kern(V);
        kern.rule_heur_2();

        ASSERT_TRUE(kern.V[4].empty());
        ASSERT_TRUE(kern.V[5].empty());
    }

    {
        clog << "Rule_heur_2: CE_test_graphs::kern_heur_2" << endl;
        VVI V = CE_test_graphs::kern_heur_2_neg;
        CEKernelizer kern(V);
        kern.rule_heur_2();

        for (int i = 0; i < V.size(); i++) {
            ASSERT_TRUE(equal(ALL(V[i]), ALL(kern.V[i]))); // nothing done by the rule
        }
    }
}

TEST(CEKernelizer, rule_heur_3) {
    {
        clog << "Rule_heur_3: CE_test_graphs::kern_heur_3, test1" << endl;
        VVI V = CE_test_graphs::kern_heur_3;
        CEKernelizer kern(V);
        kern.rule_heur_3(0,2);

        for (int i = 0; i < V.size(); i++) {
            // nothing done by the rule, because {0,1,2,3} is not a cc since {7,8,9} was not removed
            ASSERT_TRUE(equal(ALL(V[i]), ALL(kern.V[i])));
        }
    }

    {
        clog << endl << "Rule_heur_3: CE_test_graphs::kern_heur_3, test2 with strict_condition" << endl;
        VVI V = CE_test_graphs::kern_heur_3;
        CEKernelizer kern(V);
        kern.rule_heur_3(0,3);

        ASSERT_EQ( kern.inCl[5], kern.inCl[6] );

        ASSERT_EQ( kern.inCl[7], kern.inCl[8] );
        ASSERT_EQ( kern.inCl[7], kern.inCl[9] );

        function<bool(int,int)> isEqual = []( int a, int b){ return a==b; };
        ASSERT_PRED2(isEqual, kern.inCl[0], kern.inCl[1]);
        ASSERT_PRED2(isEqual,  kern.inCl[2], kern.inCl[3]);

        for( int i : {4,5,6} ) ASSERT_TRUE( kern.V[i].empty() );
        for( int i : {7,8,9} ) ASSERT_EQ( kern.V[i].size(),4 );

    }

    {
        clog << endl << "Rule_heur_3: CE_test_graphs::kern_heur_3, test2 WITHOUT strict_condition" << endl;
        VVI V = CE_test_graphs::kern_heur_3;
        CEKernelizer kern(V);
        kern.rule_heur_3(0,3,false);

        ASSERT_EQ( kern.inCl[5], kern.inCl[6] );

        ASSERT_EQ( kern.inCl[7], kern.inCl[8] );
        ASSERT_EQ( kern.inCl[7], kern.inCl[9] );

        function<bool(int,int)> isEqual = []( int a, int b){ return a==b; };
        ASSERT_PRED2(isEqual, kern.inCl[0], kern.inCl[1]);
        ASSERT_PRED2(isEqual,  kern.inCl[2], kern.inCl[3]);

        for( int i : {4,5,6,7,8,9} ) ASSERT_TRUE( kern.V[i].empty() );
    }

    {
        clog << endl << "Rule_heur_3: CE_test_graphs::kern_heur_3, test3" << endl;
        VVI V = CE_test_graphs::kern_heur_3;
        CEKernelizer kern(V);
        kern.rule_heur_3(2,3);

        for (int i = 0; i < V.size(); i++) {
            // nothing done by the rule, because {0,1,2,3} is not a cc since {4} was not removed
            ASSERT_TRUE(equal(ALL(V[i]), ALL(kern.V[i])));
        }
    }

}


TEST(CEKernelizer, rule_heur_4) {
    {
        VVI V = CE_test_graphs::kern_heur_4;
        CEKernelizer kern(V);
        kern.rule_heur_4(1,false);

        ASSERT_EQ(kern.rule_heur_4_paired_nodes.size(), 2);

        ASSERT_EQ(kern.inCl[4], kern.inCl[5]);
        ASSERT_NE(kern.inCl[4], kern.inCl[9]);

        ASSERT_EQ(kern.inCl[6], kern.inCl[7]);
        ASSERT_NE(kern.inCl[6], kern.inCl[8]);
    }

    {
        VVI V = CE_test_graphs::kern_heur_4;
        CEKernelizer kern(V);
        kern.rule_heur_4(0.67,false);

        ASSERT_EQ(kern.rule_heur_4_paired_nodes.size(), 1);

        ASSERT_NE(kern.inCl[4], kern.inCl[5]);

        ASSERT_EQ(kern.inCl[6], kern.inCl[7]);
        ASSERT_NE(kern.inCl[6], kern.inCl[8]);
    }

    {
        VVI V = CE_test_graphs::kern_heur_4;
        CEKernelizer kern(V);
        kern.rule_heur_4(0.66,false);

        ASSERT_EQ(kern.rule_heur_4_paired_nodes.size(), 0);

        ASSERT_NE(kern.inCl[4], kern.inCl[5]);

        ASSERT_NE(kern.inCl[6], kern.inCl[7]);
        ASSERT_NE(kern.inCl[6], kern.inCl[8]);
    }
}

TEST(CEKernelizer, rule_heur_5) {
    {
//        clog.rdbuf(cout.rdbuf());
        VVI V = CE_test_graphs::kern_heur_5;
        CEKernelizer kern(V);
        kern.rule_heur_5(0.5,2,2);

        int cl_id = kern.inCl[6];
        set<int> zb = {0,1,5,6,7};
        for( int i=0; i<V.size(); i++ ){
            if( zb.count(i) ) ASSERT_EQ( kern.inCl[i], cl_id );
            else ASSERT_EQ(kern.inCl[i],i);
        }
    }

    {
//        clog.rdbuf(cout.rdbuf());
        VVI V = CE_test_graphs::kern_heur_5;
        CEKernelizer kern(V);
        kern.rule_heur_5(0.75,2,2);

        int cl_id = kern.inCl[6];
        set<int> zb = {6,7};
        for( int i=0; i<V.size(); i++ ){
            if( zb.count(i) ) ASSERT_EQ( kern.inCl[i], cl_id );
            else ASSERT_EQ(kern.inCl[i],i);
        }
    }
}

TEST(CEKernelizer, rule_heur_6) {
    {
        clog.rdbuf(cout.rdbuf());
        VVI V = CE_test_graphs::kern_heur_6;
        CEKernelizer kern(V);

        int max_out_deg_sum = 3, max_out_deg = 2, max_out_deg_to_merge = 2;
        kern.rule_heur_6(max_out_deg_sum, max_out_deg, max_out_deg_to_merge);

        for( int x : VI({1,2,3}) ) ASSERT_EQ(kern.inCl[x], kern.inCl[1]);
        ASSERT_NE(kern.inCl[1], kern.inCl[0]);
        ASSERT_NE(kern.inCl[1], kern.inCl[4]);
        ASSERT_NE(kern.inCl[1], kern.inCl[5]);
        ASSERT_NE(kern.inCl[1], kern.inCl[6]);
        ASSERT_NE(kern.inCl[1], kern.inCl[7]);
    }

    {
        ENDL(1);
        VVI V = CE_test_graphs::kern_heur_6;
        CEKernelizer kern(V);

        int max_out_deg_sum = 3, max_out_deg = 2, max_out_deg_to_merge = 1;
        kern.rule_heur_6(max_out_deg_sum, max_out_deg, max_out_deg_to_merge);

        for( int x : VI({1,2}) ) ASSERT_EQ(kern.inCl[x], kern.inCl[1]);
        ASSERT_NE(kern.inCl[1], kern.inCl[0]);
        ASSERT_NE(kern.inCl[1], kern.inCl[3]);
        ASSERT_NE(kern.inCl[1], kern.inCl[4]);
        ASSERT_NE(kern.inCl[1], kern.inCl[5]);
        ASSERT_NE(kern.inCl[1], kern.inCl[6]);
        ASSERT_NE(kern.inCl[1], kern.inCl[7]);
    }

    {
        ENDL(1);
        VVI V = CE_test_graphs::kern_heur_6;
        CEKernelizer kern(V);

        int max_out_deg_sum = 4, max_out_deg = 2, max_out_deg_to_merge = 2;
        kern.rule_heur_6(max_out_deg_sum, max_out_deg, max_out_deg_to_merge);

        for( int x : VI({1,2,3,4,5}) ) ASSERT_EQ(kern.inCl[x], kern.inCl[1]);
        ASSERT_NE(kern.inCl[1], kern.inCl[0]);
        ASSERT_NE(kern.inCl[1], kern.inCl[6]);
        ASSERT_NE(kern.inCl[1], kern.inCl[7]);
    }

    {
        ENDL(1);
        VVI V = CE_test_graphs::createVVIfromVPII(CE_test_graphs::cluster_graph_test2_edges);
        CEKernelizer kern(V);

        int max_out_deg_sum = 4, max_out_deg = 2, max_out_deg_to_merge = 1;
        kern.rule_heur_6(max_out_deg_sum, max_out_deg, max_out_deg_to_merge);
    }

}

TEST(CEKernelizer, rule_heur_7_test1) {
    {
        clog.rdbuf(cout.rdbuf());
        VVI V = CE_test_graphs::createVVIfromVPII(CE_test_graphs::kern_heur_7_edges);
        CEKernelizer kern(V);

        int max_out_deg_sum = 5, max_out_deg = 3, max_out_deg_to_merge = 1;
        bool admit_diamonds = true;
        kern.rule_heur_7(max_out_deg_sum, max_out_deg, max_out_deg_to_merge, admit_diamonds);

        ASSERT_EQ( kern.inCl[9], kern.inCl[10] );
        ASSERT_EQ( kern.inCl[9], kern.inCl[14] );
        ASSERT_NE( kern.inCl[9], kern.inCl[11] );

        ASSERT_EQ( kern.inCl[4], kern.inCl[5] );
        ASSERT_EQ( kern.inCl[4], kern.inCl[8] );
        ASSERT_NE( kern.inCl[4], kern.inCl[6] );
        ASSERT_NE( kern.inCl[4], kern.inCl[7] );

        ASSERT_EQ( kern.inCl[12], kern.inCl[13] );
        ASSERT_NE( kern.inCl[12], kern.inCl[15] );
        ASSERT_NE( kern.inCl[12], kern.inCl[11] );

        ASSERT_EQ( kern.inCl[0], kern.inCl[1] );
        ASSERT_EQ( kern.inCl[0], kern.inCl[3] );
        ASSERT_NE( kern.inCl[0], kern.inCl[2] );

        ASSERT_EQ( kern.inCl[16], kern.inCl[19] );
        ASSERT_EQ( kern.inCl[16], kern.inCl[19] );
        ASSERT_NE( kern.inCl[16], kern.inCl[15] );
    }

    { // the same as above, but do not admit diamonds
        clog.rdbuf(cout.rdbuf());
        VVI V = CE_test_graphs::createVVIfromVPII(CE_test_graphs::kern_heur_7_edges);
        CEKernelizer kern(V);

        int max_out_deg_sum = 5, max_out_deg = 3, max_out_deg_to_merge = 1;
        bool admit_diamonds = false;
        kern.rule_heur_7(max_out_deg_sum, max_out_deg, max_out_deg_to_merge, admit_diamonds);

        ASSERT_EQ( kern.inCl[9], kern.inCl[10] );
        ASSERT_EQ( kern.inCl[9], kern.inCl[14] );
        ASSERT_NE( kern.inCl[9], kern.inCl[11] );

        ASSERT_EQ( kern.inCl[4], kern.inCl[5] );
        ASSERT_EQ( kern.inCl[4], kern.inCl[8] );
        ASSERT_NE( kern.inCl[4], kern.inCl[6] );
        ASSERT_NE( kern.inCl[4], kern.inCl[7] );

        ASSERT_EQ( kern.inCl[12], kern.inCl[13] );
        ASSERT_NE( kern.inCl[12], kern.inCl[15] );
        ASSERT_NE( kern.inCl[12], kern.inCl[11] );

        ASSERT_NE( kern.inCl[0], kern.inCl[1] );
        ASSERT_NE( kern.inCl[0], kern.inCl[3] );
        ASSERT_NE( kern.inCl[0], kern.inCl[2] );

        ASSERT_NE( kern.inCl[16], kern.inCl[19] );
        ASSERT_NE( kern.inCl[16], kern.inCl[19] );
        ASSERT_NE( kern.inCl[16], kern.inCl[15] );
    }
}

TEST(CEKernelizer, rule_heur_7_test2) {

    clog.rdbuf(cout.rdbuf());
    VVI V = CE_test_graphs::createVVIfromVPII(CE_test_graphs::kern_heur_7_edges);
    CEKernelizer kern(V);

    int max_out_deg_sum = 6, max_out_deg = 3, max_out_deg_to_merge = 2; // these are 'HARD' parameters
    bool admit_diamonds = true;
    kern.rule_heur_7(max_out_deg_sum, max_out_deg, max_out_deg_to_merge, admit_diamonds);

    ASSERT_EQ(kern.inCl[9], kern.inCl[10]);
    ASSERT_EQ(kern.inCl[9], kern.inCl[14]);
    ASSERT_NE(kern.inCl[9], kern.inCl[11]);

    ASSERT_EQ( kern.inCl[4], kern.inCl[5] );
    ASSERT_EQ( kern.inCl[4], kern.inCl[6] );
    ASSERT_EQ( kern.inCl[4], kern.inCl[8] );
    ASSERT_EQ( kern.inCl[4], kern.inCl[7] ); // with max_out_deg_sum = 8 we have {4,5,6,7,8} as a cluster

    ASSERT_EQ(kern.inCl[12], kern.inCl[13]);
    ASSERT_EQ(kern.inCl[12], kern.inCl[15]);
    ASSERT_NE(kern.inCl[12], kern.inCl[11]);

    ASSERT_EQ(kern.inCl[0], kern.inCl[1]);
    ASSERT_EQ(kern.inCl[0], kern.inCl[3]);
    ASSERT_EQ(kern.inCl[0], kern.inCl[2]);

    ASSERT_EQ(kern.inCl[16], kern.inCl[19]);
    ASSERT_EQ(kern.inCl[16], kern.inCl[19]);
    ASSERT_NE(kern.inCl[16], kern.inCl[15]);

}