//
// Created by sylwester on 3/13/21.
//

#include <CONTESTS/PACE21/kernelization/CEKernelizer.h>
#include <CONTESTS/PACE21/test_graphs.h>
#include <graphs/GraphUtils.h>
#include <Constants.h>
#include <utils/TimeMeasurer.h>
#include <graphs/components/ConnectedComponents.h>
#include <combinatorics/CombinatoricUtils.h>
#include <CONTESTS/PACE21/heur/ClusterGraph.h>
#include <graphs/cliques/CliqueUtils.h>
#include <utils/StandardUtils.h>
#include <CONTESTS/PACE21/heur/Global.h>
#include "utils/RandomNumberGenerators.h"
#include "CollectionOperators.h"

CEKernelizer::CEKernelizer(VVI &V, unsigned seed) :
    fau(V.size()), rnd(UniformIntGenerator(0, 1ll * 1'000'000'000 * 1'000'000'000, seed))
    {
    this->V = V;
    N = V.size();
    createHashes();

    inCl = VI(N);
    iota(ALL(inCl),0);

    inV = VB(N,true);

    helper1 = VB(N,false);
    ruleAppliedCnt = VI(MAX_RULES+MAX_HEUR_RULES+2,0);
    disabled_rules = VB(MAX_RULES+1,false);
    disabled_heur_rules = VB(MAX_HEUR_RULES+1,false);

    drng.seed( seed );
}

void CEKernelizer::createHashes() {
    hashes = VLL(N);
    LL seed = 1'718'237'121;
    for(int i=0; i<N; i++) hashes[i] = rnd.rand();
}

void CEKernelizer::markInClVector() {
    for(int i=0; i<N; i++) inCl[i] = fau.Find(i);
}

vector<CriticalClique> CEKernelizer::createCriticalCliques(VI nodes) {
    TimeMeasurer::startMeasurement("CriticalCliques");
    if( nodes.empty() ){
        nodes.reserve(N);
        for(int i=0; i<N; i++) if(inV[i]) nodes.push_back(i);
    }

    unordered_map<LL, VI, fib_hash> cc; // should be reproducible after changing fib_hash to deterministic
//    map<LL, VI> cc;

    for( int i : nodes ){
        if( V[i].empty() ){ inV[i] = false; continue; }
        LL h = hashes[i];
        for( int w : V[i] ) h ^= hashes[w];
        cc[h].push_back(i);
    }

    vector<CriticalClique> res; //(cc.size());
    res.reserve( cc.size() );
    int cnt = 0;
    for( auto & p : cc ) res.emplace_back( V, p.second, cnt++ );

    for( auto& cc : res ) for( int i=0; i<cc.nodes.size(); i++ ) fau.Union( cc.nodes[0], cc.nodes[i] );
    markInClVector();

    StandardUtils::shuffle(res,drng);

   /* for( int i=0; i<res.size(); i++ ) res[i].id = i; // random_shuffling res and changing ids
    sort( ALL(res), []( auto& c1, auto& c2 ){
        if( c1.size() != c2.size() ) return c1.size() > c2.size(); else return c1.id < c2.id;
    });*/
    { // count sort
        int M = 0;
        for( auto& c : res ) M = max(M,c.size());
        vector<vector<CriticalClique*>> v( 1 + M );
        for( int i=0; i<res.size(); i++ ) v[res[i].size()].push_back( &res[i] );

        vector<CriticalClique> res2(res.size());
        int cnt = 0;
        for( int i=M; i>=1; i-- ){
            for( auto* c : v[i] ){
                res2[cnt++] = move(*c);
            }
        }
        swap(res,res2);
    }

    for( int i=0; i<res.size(); i++ ) res[i].id = i; // assigning ids after change

    TimeMeasurer::stopMeasurement("CriticalCliques");
    return res;
}

vector<CriticalClique> CEKernelizer::createCriticalCliques(VVI &V, VI nodes) {
    TimeMeasurer::startMeasurement("CriticalCliques_2");
    int N = V.size();
    if( nodes.empty() ){
        nodes.resize(N);
        iota(ALL(nodes),0);
    }

    VLL hashes(N);
    for(int i=0; i<N; i++) hashes[i] = rnd.rand();

    unordered_map<LL, VI, fib_hash> cc; // should be reproducible after changing fib_hash to deterministic
//    map<LL, VI> cc; // make the results reproducible

    for( int i : nodes ){
        if( V[i].empty() ){ continue; }
        LL h = hashes[i];
        for( int w : V[i] ) h ^= hashes[w];
        cc[h].push_back(i);
    }

    vector<CriticalClique> res;
    res.reserve( cc.size() );
    int cnt = 0;
    for( auto & p : cc ) res.emplace_back( V, p.second, cnt++ );

    StandardUtils::shuffle(res,drng);

//    for( int i=0; i<res.size(); i++ ) res[i].id = i; // assigning ids - no need in case of count-sort
//    sort( ALL(res), []( auto& c1, auto& c2 ){
//        if( c1.size() != c2.size() ) return c1.size() > c2.size();  else return c1.id < c2.id;
//    });
    { // count sort
        int M = 0;
        for( auto& c : res ) M = max(M,c.size());
        vector<vector<CriticalClique*>> v( 1 + M );
        for( int i=0; i<res.size(); i++ ) v[res[i].size()].push_back( &res[i] );

        vector<CriticalClique> res2(res.size());
        int cnt = 0;
        for( int i=M; i>=1; i-- ){
            for( auto* c : v[i] ){
                res2[cnt++] = move(*c);
            }
        }
        swap(res,res2);
    }
    for( int i=0; i<res.size(); i++ ) res[i].id = i; // random_shuffling res

    TimeMeasurer::stopMeasurement("CriticalCliques_2");
    return res;
}



VI CEKernelizer::getNeighborhood(VI &nodes, int dst) {
    VI res = nodes;
    for( int p : res ) helper1[p] = true;

    VI next_l, prev_l = res;
    for( int i=0; i<dst; i++ ){
        next_l.clear();
        for( int p : prev_l ){
            for(int d : V[p]){
                if(!helper1[d]){
                    helper1[d] = true;
                    next_l.push_back(d);
                }
            }
        }
        prev_l = next_l;
        res += next_l;
    }

    for( int p : res ) helper1[p] = false; // clearing everything

    return next_l;
}

void CEKernelizer::makeCliqueAndRemoveFromGraph(VI nodes) {
    for( int i=0; i<nodes.size(); i++ ){
        inCl[ nodes[i] ] = inCl[nodes[0]];
        inV[nodes[i]] = false;
        fau.Union( nodes[i], nodes[0] );
    }

    GraphUtils::removeNodesFromGraph( V, nodes );
}


int CEKernelizer::getEditingDegree( int v, VB & inK, VB & inN1, int N1size ){
    assert( inN1[v] );
    int ed = N1size - 1; // the editting degree

    for( int d : V[v] ){
        if( inN1[d] ) ed--;
        else if( !inK[d] ) ed++;
    }

    return ed;
}


VVI CEKernelizer::createEmptyCriticalCliques() {
//    unordered_map<LL,VI, fib_hash> accs;
    map<LL,VI> accs;

    for( int i=0; i<N; i++ ){
        if( V[i].empty() || !inV[i] ) continue;
        LL h = 0;
        for( int w : V[i] ) h ^= hashes[w];
        accs[h].push_back(i);
    }

    VVI res;
    res.reserve(accs.size());
    for( auto& p : accs ) res.push_back(p.second);
    return res;
}


vector<PII> CEKernelizer::getEmptyCriticalCliqueSizes() {
   VVI accs = createEmptyCriticalCliques();

    unordered_map<int,int> cnt;
    for( auto& v : accs ) cnt[v.size()]++;

    VPII res( ALL(cnt) );
    sort( res.rbegin(), res.rend() );
    return res;
}

bool CEKernelizer::rule1() {
    if(disabled_rules[1]) return false;
    TimeMeasurer::startMeasurement("CEKernelizer_rule1");

    const bool debug = false;
    if(ccs.empty())
        ccs = createCriticalCliques();
    bool changes = false;

    for( auto& cc : ccs ){
//        auto neigh = GraphUtils::getNeighborhoodExclusive(V, cc.nodes, helper1);
        /*bool emptyN1 = true;
        for( int d : cc.nodes ) helper1[d] = true;
        for( int d : cc.nodes ){
            if(!emptyN1) break;
            for( int w : V[d] ){
                if( !helper1[w] ){
                    emptyN1 = false;
                    break;
                }
            }
        }
        for( int d : cc.nodes ) helper1[d] = false;*/

        // since all nodes in cc have the same closed neighborhoods, we infer that N1 is empty only if
        // V[v].size() >= K.size() for every v \in K
        bool emptyN1 = ( (V[cc.nodes[0]].size() + 1) == cc.size());

//        if( neigh.empty() ){
        if( emptyN1 ){
            if(debug) clog << "Rule1 applies to cc: " << cc << endl;
            makeCliqueAndRemoveFromGraph( cc.nodes );
            ruleAppliedCnt[1]++;
            changes = true;
        }
    }

    if(debug) DEBUG(inCl);

    if(changes) ccs.clear(); // no need to clear and create CC's again, can just remove those that were made a cluster

    markInClVector();
    TimeMeasurer::stopMeasurement("CEKernelizer_rule1");
    return changes;
}


bool CEKernelizer::rule2() {
    if(disabled_rules[2]) return false;
    TimeMeasurer::startMeasurement("CEKernelizer_rule2");

    const bool debug = false;

    if(ccs.empty())
        ccs = createCriticalCliques();
    bool changes = false;

    VB affected(N,false);

    VI temp_neigh;
    temp_neigh.reserve(N);

    for( auto& cc : ccs ){
        bool checkCC = true;
        for( int d : cc.nodes ) if(affected[d]){ checkCC = false;break; }
        if(!checkCC) continue;

        auto K = cc.nodes;

        // below is the same as condition K.size() < N1.size(), but no need to find N1 beforehand. This is because K is a
        // critical clique.  This is just    K.size() < V[K[0]].size() - ( K.size()-1 )
        if( (K.size()<<1) <= 1 + V[K[0]].size() ) continue;
//        if( K.size() <= N1.size() ) continue; // condition  K.size() > N1.size() + N2size will not be met

        auto N1 = getNeighborhood( K,1 );

        int N2size = 0;
        { // this here is just to   NOT use getNeighborhood(K,2), but just count N2size.
            temp_neigh.clear();
            auto KN1 = K+N1;
            for (int d : KN1) helper1[d] = true;
            for (int d : N1){
                for(int w : V[d]) if (!helper1[w]) {
                    helper1[w] = true; N2size++;
                    temp_neigh.push_back(w);
                }

                // in debug mode we want to calculate exact values - in non debug we can make the speedup,
                if( !debug && (K.size() <= N1.size() + N2size) ) break; // condition  K.size() > N1.size() + N2size will not be met
            }
//            for (int d : N1) for(int w : V[d]) helper1[w] = false; // clearing
            for (int w : temp_neigh) helper1[w] = false; // clearing
            for (int d : KN1) helper1[d] = false;
        }

        if( K.size() > N1.size() + N2size ){
            if(debug)  clog << "Rule2 applies to cc: " << cc << endl;
            makeCliqueAndRemoveFromGraph( K + N1 );
            changes = true;
            ruleAppliedCnt[2]++;

            // all nodes from this set are affected, cc's in N1 may become incorporated into K+N1
            for( int d : K+N1 ) affected[d] = true;
        }
    }

    if(changes) ccs.clear();
    markInClVector();

    TimeMeasurer::stopMeasurement("CEKernelizer_rule2");
    return changes;
}

bool CEKernelizer::rule3() {
    if(disabled_rules[3]) return false;
    TimeMeasurer::startMeasurement("CEKernelizer_rule3");

    const bool debug = false;

    if(ccs.empty())
        ccs = createCriticalCliques();
    bool changes = false;
    VB affected(N,false);

    VI inCC(N,-1); // id of the CC that contains node d
    for( auto& cc : ccs ){
        for( int d : cc.nodes ){
            inCC[d] = cc.id;
        }
    }

    if(debug) DEBUG(ccs);

    /**
     * cc_inters[i] is the size of the intersection of ccs[i].nodes and N1 (see below)
     */
    VI cc_inters(ccs.size(),0);

    VB inN1(N,false);
    VB inN2(N,false);
    VB inKp(N,false);

    for( auto& cc : ccs ){
        bool checkCC = true;
        for( int d : cc.nodes ) if(affected[d]){ checkCC = false;break; }
        if(!checkCC) continue;


        auto K = cc.nodes;

        // below is the same as condition K.size() < N1.size(), but no need to find N1 beforehand. This is because K is a
        // critical clique.  This is just    K.size() < V[K[0]].size() - ( K.size()-1 )
        if( (K.size()<<1) < 1 + V[K[0]].size() ){ if(debug) clog << "K.size() < N1.size()" << endl;continue; }
//        if( K.size() < N1.size() ){ if(debug) clog << "K.size() < N1.size()" << endl;continue; }


        auto N1 = getNeighborhood(K,1);
        if(debug) clog << endl << endl << "Considering cc K: " << K << endl;

        /**
         * internal critical cliques - those that have node set contained in N(K)
         */
        vector<CriticalClique*> ccs2;

        for( int d : N1 ){
            int cc_id = inCC[d];
            cc_inters[cc_id]++;

            if(cc_inters[cc_id] == ccs[cc_id].nodes.size()){
                ccs2.push_back( &ccs[cc_id] );
            }
        }

        // clearing cc_inters
        for( int d : N1 ){
            int cc_id = inCC[d];
            cc_inters[cc_id]--;
        }

        if( ccs2.empty() ) continue;

        sort( ALL(ccs2), []( auto cc1, auto cc2 ){
            if( cc1->size() != cc2->size() ) return cc1->size() > cc2->size();
            else return cc1->id < cc2->id;
        } );

        for(int d : N1) inN1[d] = true;

        auto N2 = getNeighborhood(K,2);
        for( int d : N2 ) inN2[d] = true;

        if(debug){
            clog << "Internal cliques found" << endl;
            clog << "K: " << K << endl << "N1: " << N1 << endl << "N2: " << N2 << endl;
            clog << "Critical cliques in N1: " << endl;
            for( auto* c : ccs2) clog << *c << endl;
        }

        for( auto& cc2 : ccs2 ){
            auto Kp = cc2->nodes;
            auto NKp = getNeighborhood(Kp,1);

            for( int d : Kp ) inKp[d] = true;

            if(debug) clog << endl << "Considering Kp: " << Kp << "   NKp: " << NKp << endl;

            assert( N1.size() >= Kp.size() );
            int KpN1absent = Kp.size() * (N1.size() - Kp.size() ); // number of edges needed to make Kp and N1 a clique
            int KpN2present = 0; // size of set  N(K') \cap N2(K)

            for( int d : Kp ){
                for( int w : V[d] ){
                    if( inN1[w] && !inKp[w] ) KpN1absent--;
                    if( inN2[w] ) KpN2present++;
                }
            }

            for( int d : Kp ) inKp[d] = false;

            if(KpN2present == 0 ) continue;

            if( K.size() * Kp.size() >= KpN2present + KpN1absent ){
                VPII edges_to_remove;
                for( int d : Kp ){
                    for( int w : V[d] ){
                        if( inN2[w] ) edges_to_remove.push_back({d, w} );
                    }
                }

                if(debug){
                    clog << "Rule3 applies to cliques K: " << K << " and Kp: " << Kp << endl;
                    DEBUG2(KpN1absent,KpN2present);
                    clog << "L = " << K.size() * Kp.size() << "  >=   " << KpN1absent + KpN2present << " = R" << endl;
                    clog << "Removing edges: " << edges_to_remove << endl;
                }

                GraphUtils::removeEdges(V, edges_to_remove );

                changes = true;
                ruleAppliedCnt[3]++;
                for( PII e : edges_to_remove ) affected[e.first] = affected[e.second] = true;
                { // making union of K and Kp
                    fau.Union( K[0], Kp[0] ); // this will make it, since K and Kp are critical cliques
                }
            }else{
                if(debug){
                    clog << "Rule3 does not apply to Kp: " << Kp << endl;
                    DEBUG2(KpN1absent,KpN2present);
                    clog << "L = " << K.size() * Kp.size() << "  <   " << KpN1absent + KpN2present << " = R" << endl;
                }
            }
        }

        for(int d : N1) inN1[d] = false;
        for( int d : N2 ) inN2[d] = false;

        if(debug) clog << endl;
    }

    if(changes) ccs.clear();

    markInClVector();

    TimeMeasurer::stopMeasurement("CEKernelizer_rule3");
    return changes;
}

bool CEKernelizer::rule4() {
    if(disabled_rules[4]) return false;
    TimeMeasurer::startMeasurement("CEKernelizer_rule4");

    const bool debug = false;

    if(ccs.empty())
        ccs = createCriticalCliques();
    bool changes = false;

    VB affected(N,false);

    for( auto& cc : ccs ){
        bool checkCC = true;
        for( int d : cc.nodes ) if(affected[d]){ checkCC = false;break; }
        if(!checkCC) continue;


        auto K = cc.nodes;


        // below is the same as condition K.size() < N1.size(), but no need to find N1 beforehand. This is because K is a
        // critical clique.  This is just    K.size() < V[K[0]].size() - ( K.size()-1 )
        if( (K.size()<<1) < 1 + V[K[0]].size() ) continue;
//        if( K.size() < N1.size() ) continue;

        auto N1 = getNeighborhood(K,1);
        auto N2 = getNeighborhood(K,2);

        if( N2.empty() ){
            makeCliqueAndRemoveFromGraph(K+N1);
            changes = true;
            ruleAppliedCnt[4]++;
            for( int d : K+N1 ) affected[d] = true;
        }
    }

    if(changes) ccs.clear();

    markInClVector();

    TimeMeasurer::stopMeasurement("CEKernelizer_rule4");
    return changes;
}

bool CEKernelizer::rule6_lemma3() {
    if(disabled_rules[6]) return false;
    TimeMeasurer::startMeasurement("CEKernelizer_rule6");

    const bool debug = false;

    if(ccs.empty())
        ccs = createCriticalCliques();
    bool changes = false;

    VB affected(N,false);
    VB inK(N,false), inN1(N,false);

    if(debug)
        DEBUG(ccs);

    for( auto& cc : ccs ){
        bool checkCC = true;
        for( int d : cc.nodes ) if(affected[d]){ checkCC = false;break; }
        if(!checkCC) continue;


        auto K = cc.nodes;

        if(debug){
            auto N1 = getNeighborhood(K,1);
            clog << endl; DEBUG2(K,N1);
        }

        // below is the same as condition K.size() < N1.size(), but no need to find N1 beforehand. This is because K is a
        // critical clique.  This is just    K.size() < V[K[0]].size() - ( K.size()-1 )
        if( (K.size()<<1) < 1 + V[K[0]].size() ) continue;
//        if( K.size() < N1.size() ) continue;

        auto N1 = getNeighborhood(K,1);

        for(int d : K) inK[d] = true; for(int d : N1) inN1[d] = true;

        bool can = true;
        for( int v : N1 ){
            int ed = getEditingDegree( v, inK, inN1, N1.size() );
            if(debug){
                clog << "ed(" << v << ") = " << ed << endl;
            }
            if( ed > K.size() ){
                if(debug) clog << "editing degree of node " << v << " = " << ed << " > " << K.size() << endl;

                can = false;
                break;
            }
        }

        if( can ){
            // all nodes v from N1 have degree <= K.size()
            if(debug){
                clog << "Reul6_lemma3 applies to K: " << K << "  N1: " << N1 << endl;
            }

            makeCliqueAndRemoveFromGraph(K+N1);
            for(int d : K+N1) affected[d] = true;

            ruleAppliedCnt[6]++;
            changes = true;
        }else{
            if(debug)
                clog << "Rule6_lemma3 does not apply to K: " << K << "   N1: " << N1 << endl;
        }

        for(int d : K) inK[d] = false; for(int d : N1) inN1[d] = false;
    }

    if(changes) ccs.clear();

    markInClVector();

    TimeMeasurer::stopMeasurement("CEKernelizer_rule6");
    return changes;
}

bool CEKernelizer::rule7() {
    if(disabled_rules[7]) return false;
    TimeMeasurer::startMeasurement("CEKernelizer_rule7");

    const bool debug = false;

    if(ccs.empty())
        ccs = createCriticalCliques();
    bool changes = false;

    VB affected(N,false);
    VB inK(N,false), inN1(N,false);
    VI uDegs(N,0);

    if(debug) DEBUG(ccs);

    for( auto& cc : ccs ){
        bool checkCC = true;
        for( int d : cc.nodes ) if(affected[d]){ checkCC = false;break; }
        if(!checkCC) continue;

        auto K = cc.nodes;

        // below is the same as condition K.size() >= N1.size(), but no need to find N1 beforehand. This is because K is a
        // critical clique.  This is just    K.size() >= V[K[0]].size() - ( K.size()-1 )
        if( (K.size()<<1) >= 1 + V[K[0]].size() ) continue;

        auto N1 = getNeighborhood(K,1);
//        if( K.size() >= N1.size() ) continue;

        int KNKsize = K.size() + N1.size();

        if(debug) DEBUG2(K,N1);

        for(int d : K) inK[d] = true; for(int d : N1) inN1[d] = true;

        long long sum_ed = 0;

        for( int i=0; i<N1.size(); i++ ){
            sum_ed += getEditingDegree( N1[i], inK, inN1, N1.size() );

            if(!debug && KNKsize <= sum_ed){
                // for speedup - no need to calculate editing degrees of the rest of the nodes if the threshold
                // is already achieved
                break;
            }
        }

        if(debug){
            clog << "sum_ed = " << sum_ed << ( sum_ed >= KNKsize ? " >= " : " < " )
                 << KNKsize << " = KNKsize" << endl;
        }

        if( KNKsize <= sum_ed ){
            for(int d : K) inK[d] = false; for(int d : N1) inN1[d] = false; // clearing
            continue;
        }

        VI largeDegs;
        for( int d : N1 ){
            for( int w : V[d] ){
                if( !inK[w] && !inN1[w] ){
                    uDegs[w]++;
                    if( (uDegs[w] << 1) > KNKsize ){
                        largeDegs.push_back(w);
                        uDegs[w] = -Constants::INF; // assigning large negative value so that w is not added twice
                    }
                }
                if( !debug && !largeDegs.empty() ) break; // no need to check further if largeDegs has already size = 1
            }
            if( !debug && !largeDegs.empty() ) break; // no need to check further if largeDegs has already size = 1
        }

        for( int d : N1 ) for( int w : V[d] ) if( !inK[w] && !inN1[w] ) uDegs[w] = 0; // clearing arrays

        // there mey be at most one vertex with given property. This assertion is useless with speedup optimization
        // in the form of    if( !debug && !largeDegs.empty() ) break;
        assert( largeDegs.size() <= 1 );

        if( largeDegs.empty() ){
            if(debug) clog << "Rule7 applies to K: " << K << "  N1: " << N1 << endl;
            makeCliqueAndRemoveFromGraph(K+N1);
            ruleAppliedCnt[7]++;
            changes = true;
            for(int d : K+N1) affected[d] = true;
        }
        else{
            if(debug) clog << "Rule7 does not apply to K: " << K << "   N1: " << N1 <<
                            "   largeDegs: " << largeDegs << endl;
        }

        for(int d : K) inK[d] = false; for(int d : N1) inN1[d] = false; // clearing arrays
    }

    if(changes) ccs.clear();

    markInClVector();

    TimeMeasurer::stopMeasurement("CEKernelizer_rule7");
    return changes;
}

bool CEKernelizer::rule8() {
    if(disabled_rules[8]) return false;
    TimeMeasurer::startMeasurement("CEKernelizer_rule8");

    const bool debug = false;

    if(ccs.empty())
        ccs = createCriticalCliques();
    bool changes = false;

    VB affected(N,false);
    VB inK(N,false), inN1(N,false);
    VI uDegs(N,0);

    if(debug) DEBUG(ccs);

    for( auto& cc : ccs ){
        bool checkCC = true;
        for( int d : cc.nodes ) if(affected[d]){ checkCC = false;break; }
        if(!checkCC) continue;


        auto K = cc.nodes;


        // below is the same as condition K.size() >= N1.size(), but no need to find N1 beforehand. This is because K is a
        // critical clique.  This is just    K.size() >= V[K[0]].size() - ( K.size()-1 )
        if( (K.size()<<1) >= 1 + V[K[0]].size() ) continue;
//        if( K.size() >= N1.size() ) continue;

        auto N1 = getNeighborhood(K,1);
        int KNKsize = K.size() + N1.size();

        for(int d : K) inK[d] = true; for(int d : N1) inN1[d] = true;

        if(debug) DEBUG2(K,N1);

        long long sum_ed = 0;
        for( int d : N1 ){
            sum_ed += getEditingDegree( d, inK, inN1, N1.size() );
            if(!debug && KNKsize <= sum_ed){
                // for speedup - no need to calculate editing degrees of the rest of the nodes if the threshold
                // is already achieved
                break;
            }
        }


        if(debug) clog << "KNKsize = " << KNKsize <<
                ( KNKsize > sum_ed ? " > " : " <= " ) << sum_ed << " = sum_ed" << endl;

        if( KNKsize <= sum_ed ){
            for(int d : K) inK[d] = false; for(int d : N1) inN1[d] = false; // clearing arrays
            continue;
        }

        VI largeDegs;
        for( int d : N1 ){
            for( int w : V[d] ){
                if( !inK[w] && !inN1[w] ){
                    uDegs[w]++;
                    if( 2 * uDegs[w] > KNKsize ){
                        largeDegs.push_back(w);
                        uDegs[w] = -Constants::INF; // assigning large negative value so that w is not added twice
                    }
                }
                if( !debug && !largeDegs.empty() ) break; // no need to check further if largeDegs has already size = 1
            }
            if( !debug && !largeDegs.empty() ) break; // no need to check further if largeDegs has already size = 1
        }

        for( int d : N1 ) for( int w : V[d] ) if( !inK[w] && !inN1[w] ) uDegs[w] = 0; // clearing arrays

        // there mey be at most one vertex with given property. This assertion is useless with speedup optimization
        // in the form of    if( !debug && !largeDegs.empty() ) break;
        assert( largeDegs.size() <= 1 ); // there mey be at most one vertex with given property

        if( !largeDegs.empty() ){
            int u = largeDegs[0];
            if(debug) clog << "Rule8 applies to K: " << K << "   N1: " << N1 << "   u: " << u << endl;

            VPII edges_to_add;
            VPII edges_to_remove;
            for( int d : N1 ){
                for( int w : V[d] ){
                    helper1[w] = true; // marking neighbors of d
                    if( !inK[w] && !inN1[w] && w != u ) edges_to_remove.push_back( {d,w} );
                }
                for( int d2 : N1 ) if( d < d2 && !helper1[d2] ) edges_to_add.push_back({d,d2});
                for( int w : V[d] ) helper1[w] = false; // clearing array
            }

            if(debug){
                DEBUG2(edges_to_remove, edges_to_add);
            }

            GraphUtils::removeEdges( V, edges_to_remove );
            for( auto [a,b] : edges_to_add ) GraphUtils::addEdge(V,a,b);

            for(int d : K + N1) affected[d] = true;
            for( auto [a,b] : edges_to_remove ) affected[a] = affected[b] = true;
            if( edges_to_add.size() > 0 || edges_to_remove.size() > 0 ){
                changes = true;
                ruleAppliedCnt[8]++;
            }
        }
        else{
            if(debug) clog << "Rule8 does not apply to K: " << K << "   N1: " << N1 << endl;
        }

        for(int d : K) inK[d] = false; for(int d : N1) inN1[d] = false; // clearing arrays
    }

    if(changes) ccs.clear();

    markInClVector();

    TimeMeasurer::stopMeasurement("CEKernelizer_rule8");
    return changes;
}


bool CEKernelizer::rule9() {
    if(disabled_rules[9]) return false;
    TimeMeasurer::startMeasurement("CEKernelizer_rule9");

    const bool debug = false;

    if(ccs.empty())
        ccs = createCriticalCliques();
    bool changes = false;

    VB affected(N,false);
    VB inK(N,false), inN1(N,false);

    unordered_set<LL,fib_hash> cc_hashes;
    for( auto& cc : ccs ){
        LL h = 0;
        for( int d : cc.nodes ) h ^= hashes[d];
        cc_hashes.insert(h);
    }

    if(debug) DEBUG(ccs);

    for( auto& cc : ccs ){
        bool checkCC = true;
        for( int d : cc.nodes ) if(affected[d]){ checkCC = false;break; }
        if(!checkCC) continue;


        auto K = cc.nodes;
        auto N1 = getNeighborhood(K,1);

        if(debug) DEBUG(K);

        LL N1hash = 0;
        for(int d : N1) N1hash ^= hashes[d];


        bool isCC = (cc_hashes.count(N1hash) > 0);

        if(debug && !isCC) clog << "N1 is not a critical clique" << endl;

        if(!isCC) continue; // N1 must be a single critical clique

        bool isN2Size1 = ( V[N1[0]].size() == K.size() + N1.size() ); // this is because N1 is a critical clique
//        if( N2.size() != 1 ) continue; // N2 must be just single vertex
        if( !isN2Size1 ) continue; // N2 must be just single vertex

        if(debug && !isN2Size1){
            auto N2 = getNeighborhood(K,2);
            clog << "N2 is not a single vertex, N2: " << N2 << endl;
        }

        { // checking condition about editing degrees : K.size() + N1.size() > sum_ed
            int KNKsize = K.size() + N1.size();
            for (int d : K) inK[d] = true;
            for (int d : N1) inN1[d] = true;

            if (debug) DEBUG2(K, N1);
            long long sum_ed = 0;
            for (int d : N1){
                sum_ed += getEditingDegree(d, inK, inN1, N1.size());
                if(!debug && KNKsize <= sum_ed){
                    // for speedup - no need to calculate editing degrees of the rest of the nodes if the threshold
                    // is already achieved
                    break;
                }
            }

            if (debug) clog << "KNKsize = " << KNKsize <<
                    (KNKsize > sum_ed ? " > " : " <= ") << sum_ed << " = sum_ed" << endl;

            for (int d : K) inK[d] = false; for (int d : N1) inN1[d] = false; // clearing arrays
            if (KNKsize <= sum_ed) continue;
        }

        if( K.size() >= N1.size() ){
            if(debug)
                clog << "Rule9 applies to K: " << K << "  N1: " << N1 << endl;

           // assertion below may not be true if rule9 has previously modified some neighborhood of K+N1
//            assert( (1 == 0) && "This condition should not be met, because rule6_lemma3 should apply first" );

            makeCliqueAndRemoveFromGraph(K+N1);
        }else{
            VI U = CombinatoricUtils::getRandomSubset( N1.size()-1, K.size(), rnd.nextInt(1e9) );
            for( int& d : U ) d = N1[d];
            makeCliqueAndRemoveFromGraph(K+U);
            makeCliqueAndRemoveFromGraph(K+U);
        }

        ruleAppliedCnt[9]++;
        changes = true;
        for(int d : K+N1) affected[d] = true;
    }

    if(changes) ccs.clear();
    markInClVector();

    TimeMeasurer::stopMeasurement("CEKernelizer_rule9");
    return changes;
}

bool CEKernelizer::rule15() {
    if(disabled_rules[15]) return false;
    TimeMeasurer::startMeasurement("CEKernelizer_rule15");

    const bool debug = false;
    bool changes = false;

    VI order;
    for( int i=0; i<N; i++ ){
        if( inV[i] && V[i].size() >= 2 ){
            sort( ALL( V[i] ), [&]( int a, int b ){ return V[a].size() > V[b].size(); } );
        }

        if(inV[i] && !V[i].empty()) order.push_back(i);
    }

    sort( ALL(order), [&]( int a, int b ){ return V[a].size() > V[b].size(); } );

    if(debug) DEBUG(order);

    VB was(N,false);

    for( int i=0; i<order.size(); i++ ){
        int a = order[i];

        if(debug) clog << endl << "Considering node " << a << endl;

        for( int d : V[a] ) was[d] = true;

        for( int b : V[a] ){
            if( V[b].size() < V[a].size()-1 ) break; // condition will no longer be met
            if( V[b].size() >= V[a].size() ) continue; // condition not yet met

            if( fau.Find(a) == fau.Find(b) ) continue; // if a and b are already in the same cluster.

            int cnt = 0;
            for( int w : V[b] ) if( was[w] ) cnt++;

            if(debug) clog << "\tb: " << b << ", cnt: " << cnt << endl;

            assert( cnt+2 <= V[a].size() ); // +2 because was[a] is also false

            if( cnt+2 == V[a].size() ){
                changes = true;
                fau.Union(a,b); // nodes belong to the same component
                if(debug){
                    clog << "\tRule15 applies to edge (" << a << "," << b << "), V[a]: "
                         << V[a] << ", V[b]: " << V[b] << endl;
                }
                ruleAppliedCnt[15]++;
            }else{
                if(debug) clog << "\tNode " << b << " does not have included neighborhood" << endl;
            }
        }

        for( int d : V[a] ) was[d] = false;
    }

    markInClVector();

    {
        // Filling all missing edges between nodes marked to the same cluster
        vector< unordered_set<int> > edges(N);
        for( int i=0; i<N; i++ ) edges[i].insert(ALL(V[i]));

        VVI clusters(N);
        for( int d : order ) clusters[ fau.Find(d) ].push_back(d);
        for( auto& cl : clusters ){
            if( cl.size() <= 1 ) continue;
            for( int i=0; i<cl.size(); i++ ){
                int a = cl[i];
                for( int j=i+1; j<cl.size(); j++ ){
                    int b = cl[j];
                    if( edges[a].count(b) == 0 ){ // edge (a,b) is not in the graph
                        GraphUtils::addEdge( V,a,b );
                        if(debug)
                            clog << "In rule15, adding edge (" << a << "," << b << ")" << endl;
//                        edges[a].insert(b); edges[b].insert(a); // no need in fact to add nodes a and b to [edges]
                    }
                }
            }
        }
    }

    TimeMeasurer::stopMeasurement("CEKernelizer_rule15");
    return changes;
}

bool CEKernelizer::rule16(bool debugRuleApplicationOnTheFly) {
    if( disabled_rules[16] ) return false;

    bool old_dis_h_6 = disabled_heur_rules[6];
    bool old_dis_h_7 = disabled_heur_rules[7];

    setDisableHeurRule(6,false);
    setDisableHeurRule(7,false);

    bool changes = false;
    if( rule_heur_6(2,2,2) ){
        if(debugRuleApplicationOnTheFly) clog << "Almost critical clique triangle" << endl;
        changes = true;
    }
    if( rule_heur_7(3,3,3,false) ){
        if(debugRuleApplicationOnTheFly) clog << "Almost critical clique K4" << endl;
        changes = true;
    }
    if( rule_heur_7(2,2,2,true) ){
        if(debugRuleApplicationOnTheFly) clog << "Almost critical clique Diamond" << endl;
        changes = true;
    }

    setDisableHeurRule(6,old_dis_h_6);
    setDisableHeurRule(7,old_dis_h_7);

    if(changes) ruleAppliedCnt[16]++;
    return changes;
}


bool CEKernelizer::rule_heur_1() {
    if(disabled_heur_rules[1]) return false;
    TimeMeasurer::startMeasurement("CEKernelizer_rule_heur1");

    const bool debug = false;
    bool changes = false;

    VLL closed_neigh_hash(N,0); // hash of closed neighborhood N[i]
    unordered_map<LL,int> hash_count; // number of same hashes from closed_neigh_hashes

    unordered_map<LL,VI> hash_count_debug; // number of same hashes from closed_neigh_hashes

    for( int i=0; i<N; i++ ){
        closed_neigh_hash[i] = hashes[i];
        for(int d : V[i]) closed_neigh_hash[i] ^= hashes[d];
        hash_count[ closed_neigh_hash[i] ]++;

        hash_count_debug[ closed_neigh_hash[i] ].push_back(i);
    }

    if( ccs.empty() )
        ccs = createCriticalCliques();
    if(debug) DEBUG(ccs);

    // cc_clusters[hash] is a set of CC's that have the same N1
    unordered_map<LL,vector<CriticalClique*>> cc_clusters; // should be reproducible after changing fib_hash to deterministic
//    map<LL,vector<CriticalClique*>> cc_clusters; // make results reproducible

    for( auto& cc : ccs ){
        auto K = cc.nodes;
        auto N1 = getNeighborhood(K,1);
        LL h = 0;
        for( int d : N1 ) h ^= hashes[d];

        cc_clusters[h].push_back(&cc);
    }

    if(debug){
        clog << "cc_clusters:" << endl;
        for( auto & [hash, vec] : cc_clusters ){
            for( auto* c : vec ) clog << *c << "  |  ";  clog << endl;
        }
        ENDL(1);
    }

    VB inN1(N,false);

    // making the results reproducible (on map instead on unordered_map) but considering CCs in random order
    vector< vector<CriticalClique*> > cc_clusters_vec;
    cc_clusters_vec.reserve(cc_clusters.size());
    for (auto &[hash, vec] : cc_clusters) cc_clusters_vec.push_back(vec);
    StandardUtils::shuffle(cc_clusters_vec,drng);

    for( auto & vec : cc_clusters_vec ){
        VI N1 = getNeighborhood(vec[0]->nodes,1);

        if(debug){
            clog << "cc_cluster:  "; for( auto* c : vec ) clog << *c << "  |  ";  clog << endl;
            DEBUG(N1);
        }

        for( int d : N1 ) inN1[d] = true;

        for( CriticalClique * cc : vec ){
            for( int d : cc->nodes ){
                for( int w : V[d] ){
                    if( inN1[w] ) closed_neigh_hash[w] ^= hashes[d];
                }
            }
        }

        unordered_set<LL,fib_hash> hashes_after_change;
        for( int d : N1 ){
            hashes_after_change.insert( closed_neigh_hash[d] );
            if( hashes_after_change.size() >= 2 ) break;
        }

        if( hashes_after_change.size() == 1 ){ // all nodes from N1 have the same closed neighborhoods
            LL h = *hashes_after_change.begin();

            int KpSize = N1.size() + hash_count[h]; // size of the critical clique after removing current cc_cluster

            if(debug){
                clog << "N1 is contained in a critical clique of size " << KpSize << endl;
                clog << "That critical clique is " << (N1 + hash_count_debug[h]) << endl;
            }

            for( auto* c : vec ){
                if( 2 * c->size() <= KpSize && 2 * N1.size() <= KpSize ){
                    if(debug){
                        clog << "!! RULE APPLIES!  Making a clique from " << *c << endl;
                    }

                    makeCliqueAndRemoveFromGraph( c->nodes );

                    // need to update hashes for nodes in N1 and hash_count map, as c->nodes is permanently removed
                    // and the hashes will be once again changed at the end of the main loop
                    for( int d : N1 ) hash_count[ closed_neigh_hash[d] ]--;
                    for( int d : c->nodes ){
                        hash_count[ closed_neigh_hash[d] ]--;
                        for( int w : V[d] ){
                            if( inN1[w] ) closed_neigh_hash[w] ^= hashes[d];
                        }
                    }
                    for( int d : N1 ) hash_count[ closed_neigh_hash[d] ]++;

                    changes = true;
                    ruleAppliedCnt[MAX_RULES+1]++;
                }
            }
        }else{
            if(debug){
                clog << "cc_cluster:  "; for( auto* c : vec ) clog << *c << "  |  ";  clog << endl;
                DEBUG(N1);
                clog << "N1 is not a critical clique after removing cc_cluster" << endl << endl;
            }
        }

        for( CriticalClique * cc : vec ) for( int d : cc->nodes ) for( int w : V[d] )
            if( inN1[w] ) closed_neigh_hash[w] ^= hashes[d]; // resetting hash to its original value

        for( int d : N1 ) inN1[d] = false;
    }

    if(changes) ccs.clear();
    markInClVector();

    TimeMeasurer::stopMeasurement("CEKernelizer_rule_heur1");
    return changes;
}

bool CEKernelizer::rule_heur_2() {
    if(disabled_heur_rules[2]) return false;
    TimeMeasurer::startMeasurement("CEKernelizer_rule_heur2");

    const bool debug = false;
    bool changes = false;

    VI A;

    {
        // finding all nodes u for which N(u) is a clique
        for( int i=0; i<N; i++ ){
            if( V[i].empty() || !inV[i] ) continue;
            bool isClq = CliqueUtils::isClique( V, V[i], helper1 );
            if( isClq ) A.push_back(i);
        }
    }

    if(debug) DEBUG(A);

    auto G1 = V;
    GraphUtils::removeNodesFromGraph(G1,A);
    if(debug) DEBUG(G1);

    auto ccs_g1 = createCriticalCliques( G1 );


    if(debug) DEBUG(ccs_g1);

    VB inCCg1(N,-1);
    for( int i=0; i<ccs_g1.size(); i++ ){
        for( int d : ccs_g1[i].nodes ) inCCg1[d] = i;
    }

    for( int a : A ){
        int cc_id = inCCg1[ V[a][0] ];
        bool has_neighborhood_in_cc = true;
        for( int w : V[a] ){
            if( inCCg1[w] != cc_id ){ has_neighborhood_in_cc = false; break; }
        }

        if(!has_neighborhood_in_cc || cc_id == -1) continue;

        int KpSize = ccs_g1[cc_id].size();

        if( 2 * V[a].size() <= KpSize ){
            if(debug) clog << "Rule_heur_2 applies! Making node " << a << " a single cluster" << endl;

            makeCliqueAndRemoveFromGraph( {a} );
            changes = true;
            ruleAppliedCnt[MAX_RULES+2]++;
        }else{
            if(debug)clog << "Node " << a << " has V[" << a << "].size() = "
                        << V[a].size() << " > " << KpSize << " / 2  - it cannot be made a single cluster" << endl;
        }
    }

    if(changes) ccs.clear();

    TimeMeasurer::stopMeasurement("CEKernelizer_rule_heur2");
    return changes;
}

bool CEKernelizer::rule_heur_3(int min_cc_size, int max_cc_size, bool apply_strict_size_condition) {
    if(disabled_heur_rules[3]) return false;
    TimeMeasurer::startMeasurement("CEKernelizer_rule_heur3");

    const bool debug = false;
    bool changes = false;


    if( ccs.empty() ) ccs = createCriticalCliques();

    if(debug) DEBUG(ccs);

    vector<CriticalClique*> size_filtered_ccs;
   /* {
        for (auto &cc : ccs) size_filtered_ccs.push_back(&cc);
        size_filtered_ccs.resize(
                remove_if(ALL(size_filtered_ccs),
                          [&max_cc_size, &min_cc_size](auto &cc) {
                              return !(cc->size() >= min_cc_size && cc->size() <= max_cc_size);
                          })
                - size_filtered_ccs.begin());
    }*/
    {
        // this is just another option - all cc's are already sorted by their size, can use binary search to quickly find if exist
        auto it = lower_bound( ALL(ccs), max_cc_size, []( auto c1, int s ){ return c1.size() > s; } );
        for( ; it != ccs.end(); ++it ){
            if( (*it).size() < min_cc_size ) break;
            else size_filtered_ccs.push_back(&(*it));
        }
    }

    if( size_filtered_ccs.empty() ){
        TimeMeasurer::stopMeasurement("CEKernelizer_rule_heur3");
        return false;
    }


    vector<CriticalClique*> A;
    {
        // finding all nodes u for which N(u) is a clique
        for( auto& cc : size_filtered_ccs ){
            VI N1 = getNeighborhood( cc->nodes,1 );
            bool isClq = CliqueUtils::isClique( V, N1, helper1 );
            if(isClq) A.push_back( cc );
        }
    }

    if(debug){ clog << "A: ";for( auto a : A ) clog << *a << ", "; clog << endl; }

    auto G1 = V;

    VI to_remove;
    for_each( ALL(A), [&to_remove]( auto cc ){ to_remove += cc->nodes; } );
    GraphUtils::removeNodesFromGraph(G1,to_remove); // removing nodes from all cc's in A
    if(debug) DEBUG(G1);

    auto ccs_g1 = createCriticalCliques( G1 );


    if(debug) DEBUG(ccs_g1);

    VB inCCg1(N,-1);
    for( int i=0; i<ccs_g1.size(); i++ ){
        for( int d : ccs_g1[i].nodes ) inCCg1[d] = i;
    }

    for( CriticalClique* c : A ){
        VI N1 = getNeighborhood( c->nodes,1 );
        if( N1.empty() ) continue; // this should not happen if rule1 was applied before

        if(debug) clog << "Considering cc: " << *c << endl;

        int cc_id = inCCg1[ N1[0] ];
        bool has_neighborhood_in_cc = true;
        for( int w : N1 ){
            if( inCCg1[w] != cc_id ){ has_neighborhood_in_cc = false; break; }
        }

        if(!has_neighborhood_in_cc || cc_id == -1){
            if(debug) clog << "cc " << *c << " does not have all nodes in a critical clique in G1" << endl;
            continue;
        }

        int KpSize = ccs_g1[cc_id].size();

        if( apply_strict_size_condition ){
            if(2*c->size() > KpSize ){
                if (debug)
                    clog << "N1 is contained in critical clique " << ccs_g1[cc_id] << ", but condition "
                         << 2 * c->size() << " = 2*c->size() <= KpSize = " << KpSize << " is not met" << endl;
                continue;
            }
        }else{
            if(c->size() > KpSize ){
                if (debug)
                    clog << "N1 is contained in critical clique " << ccs_g1[cc_id] << ", but condition "
                         << c->size() << " = c->size() <= KpSize = " << KpSize << " is not met" << endl;
                continue;
            }
        }

        if( 2 * N1.size() <= KpSize ){
            if(debug) clog << "Rule_heur_3 applies! Making critical clique " << *c << " a single cluster" << endl;

            makeCliqueAndRemoveFromGraph( c->nodes );
            changes = true;
            ruleAppliedCnt[MAX_RULES+3]++;
        }else{
            if(debug)clog << "Critical clique " << *c << " has N1.size() = "
                          << N1.size() << " > " << KpSize << " / 2  - it cannot be made a single cluster" << endl;
        }
    }

    if(changes) ccs.clear();

    TimeMeasurer::stopMeasurement("CEKernelizer_rule_heur3");
    return changes;
}

bool CEKernelizer::iterative_rule_heur_3() {
    const bool debugRuleApplicationOnTheFly = false;

    const int MAXIMAL_MIN_S = 3; // this value can be larger (e.g. 50) for smaller graphs (e.g. N < 1000)
    const int MAX_S_RELATIVE_RANGE = 3; // this value can be larger (e.g. 50) for smaller graphs (e.g. N < 1000)

    auto createMaxsSet = [v = &this->V](int min_s){
        VI res(MAX_S_RELATIVE_RANGE);
        iota(ALL(res),min_s);
        res.push_back((*v).size());
        return res;
    };

    { // heur_rule3 for a few params, with |K_i| <= |K| / 2 condition
        for (int min_s = MAXIMAL_MIN_S; min_s > 0; min_s--) {
            auto MSET = createMaxsSet(min_s);
            for (int max_s : MSET) {
                if (rule_heur_3(min_s, max_s, true)) {
                    if (debugRuleApplicationOnTheFly)
                        clog << "Rule heur_3 applies, size_condition: true, min_s: " << min_s << ", max_s: " << max_s << endl;
                    return true;
                }
            }
        }
    }


//    bool USE_WITHOUT_STRICT_SIZE_CONDITION = false;
    if(USE_HEUR_RULE_3_WITHOUT_STRICT_SIZE_CONDITION){ // heur_rule3 for a few params, WITHOUT |K_i| <= |K| / 2 condition
        for (int min_s = MAXIMAL_MIN_S; min_s >= 0; min_s--) {
            auto MSET = createMaxsSet(min_s);
            for (int max_s : MSET) {
                if (rule_heur_3(min_s, max_s, false)) {
                    if (debugRuleApplicationOnTheFly)
                        clog << "Rule heur_3 applies, size_condition: false, min_s: " << min_s << ", max_s: " << max_s << endl;
                    return true;
                }
            }
        }
    }

    return false;
}

bool CEKernelizer::rule_heur_4(double perc, const bool add_edge) {
    if(disabled_heur_rules[4]) return false;
    TimeMeasurer::startMeasurement("CEKernelizer_rule_heur4");
    assert( (0 <= perc) && (perc <= 1) );

    const bool debug = false;
    bool changes = false;

    unordered_map<LL,VI,fib_hash> open_neigh_hash;
    for( int i=0; i<N; i++ ){
        if( V[i].empty() || !inV[i] ) continue;
        LL h = 0;
        for( int d : V[i] ) h ^= hashes[d];
        open_neigh_hash[h].push_back(i);
    }


    for( auto & [hash,A] : open_neigh_hash ){
        int x = perc * A.size();

        if(debug) clog << "Considering set A: " << A << ", size: " << A.size() << ", x: " << x << endl;

        if( x == 0 ) continue;

        int p=1;
        while( p < x ){
            int u = A[p-1], v = A[p];
            if( fau.Find(u) != fau.Find(v) ) {
                if (debug) clog << "Pairing nodes {" << u << "," << v << "}" << endl;
                fau.Union(u, v);
                rule_heur_4_paired_nodes.emplace_back(u,v);
                if (add_edge) GraphUtils::addEdge(V, u, v);
                changes = true;
                ruleAppliedCnt[MAX_RULES+4]++;
            }
            p += 2;
        }
    }

    if(changes) ccs.clear();
    markInClVector();

    TimeMeasurer::stopMeasurement("CEKernelizer_rule_heur4");
    return changes;
}

bool CEKernelizer::rule_heur_5(double perc, int max_absolute_diff, int min_common_neighbors) {
    if(disabled_heur_rules[5]) return false;
    TimeMeasurer::startMeasurement("CEKernelizer_rule_heur5");

    const bool debug = false;
    bool changes = false;

    VI order;
    for( int i=0; i<N; i++ ){
        if( inV[i] && V[i].size() >= 2 ){
            sort( ALL( V[i] ), [&]( int a, int b ){ return V[a].size() > V[b].size(); } );
        }
        if(inV[i] && !V[i].empty()) order.push_back(i);
    }

    sort( ALL(order), [&]( int a, int b ){ return V[a].size() > V[b].size(); } );

    if(debug) DEBUG(order);

    VB was(N,false);
    VB was2(N,false);

    for( int i=0; i<order.size(); i++ ){
        int a = order[i];

        if(debug) clog << endl << "Considering node " << a << endl;

//        was[a] = true;
        for( int d : V[a] ) was[d] = true;
        was2[a] = true;

        for( int b : V[a] ){
            if(was2[b]) continue;
            if( fau.Find(a) == fau.Find(b) ) continue; // if a and b are already in the same cluster.

            if(debug) DEBUG(b);

            int cnt = 0;
            for( int w : V[b] ) if( was[w] ) cnt++;
            if( cnt+1 < V[b].size() ){ continue; /* requiring N(b) \subset N[a]*/ }
            if( cnt < min_common_neighbors ){ continue; /* there must be at least min_common_neighbors*/ }

            assert(V[b].size() <= V[a].size());
            int diff = ( V[a].size() - V[b].size() );
            double thr = perc * (V[a].size()-1);

            if(debug){ clog << "\tb: " << b << ", cnt: " << cnt << endl;DEBUG(diff);DEBUG(thr); }

            { // perhaps some restriction on the structure of N(a) or N(b) or the interaction between those sets ?
                // TO BE IMPLEMENTED - perhaps the number of edges with one end in N(a) and the other in N(b) should be
                // at least diff?
            }

            if( (cnt > thr ) && (diff <= max_absolute_diff) ){
                changes = true;
                fau.Union(a,b); // nodes belong to the same component
                if(debug){
                    clog << "\tRule_heur_5 applies to edge (" << a << "," << b << "), V[a]: "
                         << V[a] << ", V[b]: " << V[b] << endl;
                    DEBUG(diff);DEBUG(thr);
                }
                ruleAppliedCnt[MAX_RULES+5]++;
            }

            if(debug) ENDL(1);
        }

        if(debug) ENDL(5);
//        was[a] = false;
        for( int d : V[a] ) was[d] = false;
    }

    markInClVector();

    TimeMeasurer::stopMeasurement("CEKernelizer_rule_heur5");
    return changes;
}

bool CEKernelizer::rule_heur_6(int max_out_deg_sum, int max_out_deg, int max_out_deg_to_merge) {
    if(disabled_heur_rules[6]) return false;
    bool changes = false;
    const bool debug = false;
    TimeMeasurer::startMeasurement("CEKernelizer_rule_heur6");

    max_out_deg_sum += 6;
    max_out_deg += 2;
    max_out_deg_to_merge += 2;

    VVI V = this->V; // making a copy
    int N = (int)V.size();

    VI order(N);
    iota(ALL(order),0);
    sort(ALL(order), [&]( int a, int b ){
        if(V[a].size() != V[b].size()) return V[a].size() < V[b].size();
        else return a < b;
    } );
    VI inOrder(N);
    for(int i=0; i<N; i++) inOrder[order[i]] = i;

    int sq = (int)ceil( sqrt( 2 * GraphUtils::countEdges(V) ) ); // just to check, can be removed later
    for( int i=0; i<N; i++ ){ // removing all edges (a,b) such that b is before a in order
        for( int k = (int)V[i].size() - 1; k>=0; k-- ){
            if( inOrder[V[i][k]] < inOrder[i] ){
                swap( V[i][k], V[i].back() );
                V[i].pop_back();
            }
        }
        sort( ALL(V[i]), [&](int a, int b){
            if(inOrder[a] != inOrder[b]) return inOrder[a] < inOrder[b];
            else return a < b;
        } );

        assert( V[i].size() <= sq );
    }

    VB was(N,false);
    for( int a : order ){
        if( this->V[a].size() > max_out_deg ) continue;

        for( int d : V[a] ) was[d] = true;

        for( int b : V[a] ){
            if( this->V[b].size() > max_out_deg ) continue;
            if( this->V[a].size() + this->V[b].size() > max_out_deg_sum ) continue;

            for( int c : V[b] ){
                if(!was[c]) continue;
                if( this->V[c].size() > max_out_deg ) continue;
                if( this->V[a].size() + this->V[b].size() + this->V[c].size() > max_out_deg_sum ) continue;

                int last = -1;
//                bool remove_structure_from_graph = false;
                for( int x : VI({a,b,c}) ){
                    if( this->V[x].size() <= max_out_deg_to_merge ){
                        if( last == -1 ) last = x;
                        else if( fau.Find(x) != fau.Find(last) ){
                            fau.Union(last,x);
                            ruleAppliedCnt[MAX_RULES+6]++;
                            changes = true;
//                            remove_structure_from_graph = true;
                            if(debug){
                                clog << "Rule heur6 applies to (a,b,c) : (" << a << "," << b << "," << c
                                     << ") --> Merging nodes " << x << " and " << last << endl;
                            }
                        }
                    }
                }

//                if( remove_structure_from_graph ){ // #TEST #CAUTION #FIXME:remove temporary for check
//                    GraphUtils::removeNodesFromGraph(V,VI({a,b,c}));
//                }
            }
        }

        for( int d : V[a] ) was[d] = false;
    }

    markInClVector();

    TimeMeasurer::stopMeasurement("CEKernelizer_rule_heur6");
    return changes;
}

bool CEKernelizer::rule_heur_7(int max_out_deg_sum, int max_out_deg, int max_out_deg_to_merge, bool enable_diamonds) {
    if(disabled_heur_rules[7]) return false;
    bool changes = false;
    const bool debug = false;
    TimeMeasurer::startMeasurement("CEKernelizer_rule_heur7");

    VI order(N,0);
    iota(ALL(order),0);

    VB was(N,false), was2(N,false);
    for( int a : order ){
        if( V[a].size() > 3 + max_out_deg ) continue;

        if(debug) clog << "a: " << a << endl;

        for( int d : V[a] ) was[d] = true; // was marks neighbors of a

        for( int b : V[a] ){ // a and b form the 'diagonal' of the diamond or K4
            if( a > b ) continue; // consider 'diagonal' edge just once
            if( V[b].size() > 3 + max_out_deg ) continue;
            if( V[a].size() + V[b].size() > 6 + max_out_deg_sum ) continue;

            if(debug) clog << "\tb: " << b << endl;

            VI common_neigh;
            for( int c : V[b] ) {
                if (!was[c]) continue;
                else common_neigh.push_back(c);
            }

            if(debug) clog << "\tcommon_neigh: " << common_neigh << endl;
            if(common_neigh.size() <= 1) continue;

            for( int i=0; i<common_neigh.size(); i++ ){
                int c = common_neigh[i];
                if( V[c].size() > 3 + max_out_deg ) continue;

                if(debug) clog << "\t\tc: " << c << endl;

                for( int d : V[c] ) was2[d] = true;

                for( int j=i+1; j<common_neigh.size(); j++ ){
                    int c2 = common_neigh[j];
                    if( V[c2].size() > 3 + max_out_deg ) continue;

                    if(debug) clog << "\t\tc2: " << c2 << endl;

                    int outdegsum = V[a].size() + V[b].size() + V[c].size() + V[c2].size() - 12;
                    bool is_diamond = false;

                    if(!was2[c2]){ // a,c,b,c2 is a diamond with diagonal a,b
                        if( V[c].size() > 2 + max_out_deg ) continue;
                        if( V[c2].size() > 2 + max_out_deg ) continue;
                        outdegsum+=2;
                        is_diamond = true;
                    }

                    if( !enable_diamonds && is_diamond ) continue;

                    if( outdegsum > max_out_deg_sum ){
                        if(debug) clog << "\t\toutdegsum = " << outdegsum << " > " << max_out_deg_sum
                                    << " = max_out_deg_sum" << endl;
                        continue;
                    }else if(debug)  clog << "\t\toutdegsum = " << outdegsum << endl;

                    if( debug ){
                        clog << "\t\t (" << a << "," << c << "," << b << "," << c2 << ") "
                            << (is_diamond ? "is DIAMOND with diagonal a-b" : "is K4") << endl;
                    }

                    int last = -1;
//                    bool remove_structure_from_graph = false;
                    for( int x : VI({a,b,c,c2}) ){
                        if( this->V[x].size() <= max_out_deg_to_merge + ( is_diamond ? ( ( x == c || x == c2 ) ? 2 : 3 ) : 3 ) ){
                            if( last == -1 ) last = x;
                            else if( fau.Find(x) != fau.Find(last) ){
                                fau.Union(last,x);
                                ruleAppliedCnt[MAX_RULES+7]++;
                                changes = true;
//                                remove_structure_from_graph = true;
                                if(debug){
                                    clog << "\t\t\tRule heur7 applies to (a,c,b,c2) : (" << a << "," << c << "," << b
                                         << ", " << c2 << ") --> Merging nodes " << x << " and " << last << endl;
                                }
                            }
                        }
                    }

//                    if( remove_structure_from_graph ){ // #TEST #CAUTION #FIXME:remove temporary for check
//                        GraphUtils::removeNodesFromGraph(V,VI({a,b,c}));
//                    }
                }

                for( int d : V[c] ) was2[d] = false;
            }
        }

        for( int d : V[a] ) was[d] = false;
    }

    markInClVector();

    TimeMeasurer::stopMeasurement("CEKernelizer_rule_heur7");
    return changes;
}

void CEKernelizer::fullKernelization(bool use_heuristic_rules, int additional_randomized_iterations) {
    bool debug = true;
    bool debugRuleApplicationOnTheFly = false;

    if(Global::disable_all_logs) debug = debugRuleApplicationOnTheFly = false; // #TEST - using Global in Kernelizer

    if(debug){
        clog << "Kernelization -->  heur: " << use_heuristic_rules
             << ",   iters: " << additional_randomized_iterations << endl;
        clog << "Before: Nodes: " << N << "\t\tEdges: " << GraphUtils::countEdges(V) << endl;
    }

    TimeMeasurer::startMeasurement("FullKernelization");

    bool changes = true;
    while( changes  ){
        if(Global::checkTle()) break;

        changes = false;
//        for(int i=0; i<N; i++) assert( helper1[i] == false );

        if( rule1() ){
            if(debugRuleApplicationOnTheFly) clog << "Rule 1 applies" << endl;
            changes = true; continue;
        }
        if( rule2() ){
            if(debugRuleApplicationOnTheFly) clog << "Rule 2 applies" << endl;
            changes = true; continue;
        }
        if( rule16(debugRuleApplicationOnTheFly) ){ // almost critical clique for C3, D3 and K4.
            // This rule is very fast, hence it can be done before other rules are checked.
            if(debugRuleApplicationOnTheFly) clog << "Rule 16 applies" << endl;
            changes = true; continue;
        }
        if( rule3() ){
            if(debugRuleApplicationOnTheFly) clog << "Rule 3 applies" << endl;
            changes = true; continue;
        }
        if( rule4() ){
            if(debugRuleApplicationOnTheFly) clog << "Rule 4 applies" << endl;
            changes = true; continue;
        }
        if(!disabled_rules[3] && rule6_lemma3() ){
            // if rule3 did not apply, then rule 6 will not either, so there is no need to check it
            if(debugRuleApplicationOnTheFly) clog << "Rule 6_lemma3 applies" << endl;
            changes = true; continue;
        }
        if( rule7() ){
            if(debugRuleApplicationOnTheFly) clog << "Rule 7 applies" << endl;
            changes = true; continue;
        }
        if( rule8() ){
            if(debugRuleApplicationOnTheFly) clog << "Rule 8 applies" << endl;
            changes = true; continue;
        }
        if( rule9() ){
            if(debugRuleApplicationOnTheFly) clog << "Rule 9 applies" << endl;
            changes = true; continue;
        }
        if( rule15() ){
            if(debugRuleApplicationOnTheFly) clog << "Rule 15 applies" << endl;
            changes = true; continue;
        }
        if( rule16(debugRuleApplicationOnTheFly) ){ // almost critical clique for C3, D3 and K4.
            if(debugRuleApplicationOnTheFly) clog << "Rule 16 applies" << endl;
            changes = true; continue;
        }


        // applying heuristic rules only if normal rules made no more progress
        if(use_heuristic_rules) {
           /* if (rule_heur_1()) {
                if (debugRuleApplicationOnTheFly) clog << "Rule heur_1 applies" << endl;
                changes = true; continue;
            }
            if (rule_heur_2()) {
                if (debugRuleApplicationOnTheFly) clog << "Rule heur_2 applies" << endl;
                changes = true; continue;
            }
            if (iterative_rule_heur_3()) {
                if (debugRuleApplicationOnTheFly) clog << "Rule heur_3 applies" << endl;
                changes = true; continue;
            }
            if (rule_heur_4()) {
                if (debugRuleApplicationOnTheFly) clog << "Rule heur_4 applies" << endl;
                changes = true; continue;
            }
            if (rule_heur_5()) {
                if (debugRuleApplicationOnTheFly) clog << "Rule heur_5 applies" << endl;
                changes = true; continue;
            }
            if (rule_heur_6()) {
                if (debugRuleApplicationOnTheFly) clog << "Rule heur_6 applies" << endl;
                changes = true; continue;
            }
            if (rule_heur_7()) {
                // default: max_out_deg_sum = 4, max_out_deg = 2, max_out_deg_to_merge = 1, enable_diamonds = true
                // for harder kernelization use (6,3,2,true)
                if (debugRuleApplicationOnTheFly) clog << "Rule heur_7 applies" << endl;
                changes = true; continue;
            }*/

           changes = applyHeuristicKernelization(debugRuleApplicationOnTheFly);
        }
    }

    auto getName = [](int i){
        string name = "CEKernelizer_rule";
        int a = i;
        if( a > MAX_RULES ){
            name += "_heur";
            a -= MAX_RULES;
        }
        name += to_string(a);
        return name;
    };

    if(debug){
        VPII apps;
        for( int i=1; i<ruleAppliedCnt.size(); i++ ){
            if( ruleAppliedCnt[i] > 0 || TimeMeasurer::getMeasurementTimeInSeconds( getName(i) ) > 0.1 )
                apps.push_back( {i, ruleAppliedCnt[i]} );
        }

        int AN = 0; for(int i=0; i<N; i++) if(inV[i]) AN++;
        clog << "After:  Nodes: " << AN << "\t\tEdges: " << GraphUtils::countEdges(V) << endl;

        { // ClusterGraph data
            VI part = inCl;
            for( int i=0; i<N; i++ ) if( V[i].empty() ) part[i] = -1;

            ClusterGraph g(&V,part);
            clog << "Cl. g.: Nodes: " << g.V.size() << "\t\tEdges: " << GraphUtils::countEdges(g.V) << endl;
            auto cluster_sizes_count = g.getClusterSizesCount();
            clog << "cluster sizes count (size,cnt): " << cluster_sizes_count << endl;
        }

        if( !apps.empty() ){
            auto old_prec = clog.precision();
            clog << "Rules used (id,cnt,time): ";
            for( auto [a,b] : apps ){
                string name = "CEKernelizer_rule";
                if( a > MAX_RULES ){
                    name += "_heur";
                    name += to_string(a-MAX_RULES);
                }else{
                    name += to_string(a);
                }
                double sec = TimeMeasurer::getMeasurementTimeInSeconds(name);
                clog << fixed; clog.precision(1);
                clog << "(" << a << " , " << b << " , " << sec << ")  |  ";
            }
            clog << endl;
            clog.precision(old_prec);
        }
    }

    for( int i=0; i<MAX_RULES+MAX_HEUR_RULES+2; i++ ){
        string name = getName(i);
        TimeMeasurer::resetOption(name);
    }

    TimeMeasurer::stopMeasurement("FullKernelization");

    if(!Global::disable_all_logs) {
        clog << "Time CriticalCliques (s): " << TimeMeasurer::getMeasurementTimeInSeconds("CriticalCliques") << endl;
        clog << "Time (s): " << TimeMeasurer::getMeasurementTimeInSeconds("FullKernelization") << endl;
        clog << "From which heuristic part: (s): "
             << TimeMeasurer::getMeasurementTimeInSeconds("HeuristicKernelization") << endl << endl;
    }

    { // stopping and clearing
        TimeMeasurer::resetOption("CriticalCliques");
        TimeMeasurer::resetOption("FullKernelization");
        TimeMeasurer::resetOption("HeuristicKernelization");
    }
}

void CEKernelizer::improveKernelizationUsingHeuristicRules() {
    const bool debug = true;
    const bool debugRuleApplicationOnTheFly = true;
    bool changes = false;
    bool iter_changes = true;

    int AN = 0; for(int i=0; i<N; i++) if(inV[i]) AN++;
    int EE = GraphUtils::countEdges(V);


//    TimeMeasurer::startMeasurement("HeuristicKernelization");

   /* while( iter_changes ) {
        iter_changes = false;

        if (rule_heur_1()) {
            if (debugRuleApplicationOnTheFly) clog << "Rule heur_1 applies" << endl;
            iter_changes = true;
        }
        if (!iter_changes && rule_heur_2()) {
            if (debugRuleApplicationOnTheFly) clog << "Rule heur_2 applies" << endl;
            iter_changes = true;
        }
        if (!iter_changes && iterative_rule_heur_3()) {
            if (debugRuleApplicationOnTheFly) clog << "Rule heur_3 applies" << endl;
            iter_changes = true;
        }
        if (!iter_changes && rule_heur_4()) {
            if (debugRuleApplicationOnTheFly) clog << "Rule heur_4 applies" << endl;
            iter_changes = true;
        }
        if (!iter_changes && rule_heur_5()) {
            if (debugRuleApplicationOnTheFly) clog << "Rule heur_5 applies" << endl;
            iter_changes = true;
        }
        if (!iter_changes && rule_heur_6()) {
            if (debugRuleApplicationOnTheFly) clog << "Rule heur_6 applies" << endl;
            iter_changes = true;
        }
        if (!iter_changes && rule_heur_7()) {
            if (debugRuleApplicationOnTheFly) clog << "Rule heur_7 applies" << endl;
            iter_changes = true;
        }

//        changes = iter_changes;
        if(iter_changes) changes = true;
        if( !use_heuristic_rules_separately ) break;
    }*/

   changes = applyHeuristicKernelization(debugRuleApplicationOnTheFly);

//    TimeMeasurer::stopMeasurement("HeuristicKernelization");

    if(changes){
        if(debug) clog << "Re-running full kernelization (check in "
                    << TimeMeasurer::getMeasurementTimeInSeconds("HeuristicKernelization") << " sec.)" << endl;

        TimeMeasurer::resetOption("HeuristicKernelization");

        fullKernelization(true, 0);
    }
    else if(debug){
        clog << "Kernelization -->  heur: " << 1 << ",   iters: " << 0 << endl;
        clog << "Before: Nodes: " << AN << "\t\tEdges: " << EE << endl;

        int AN = 0; for(int i=0; i<N; i++) if(inV[i]) AN++;
        clog << "After:  Nodes: " << AN << "\t\tEdges: " << GraphUtils::countEdges(V) << endl;
        { // ClusterGraph data
            VI part = inCl;
            for( int i=0; i<N; i++ ) if( V[i].empty() ) part[i] = -1;

            ClusterGraph g(&V,part);
            clog << "Cl. g.: Nodes: " << g.V.size() << "\t\tEdges: " << GraphUtils::countEdges(g.V) << endl;
            auto cluster_sizes_count = g.getClusterSizesCount();
            clog << "cluster sizes count (size,cnt): " << cluster_sizes_count << endl;
        }
        clog << "No heuristic improvement (check in "
            << TimeMeasurer::getMeasurementTimeInSeconds("HeuristicKernelization") << " sec.)" << endl;
        TimeMeasurer::resetOption("HeuristicKernelization");
    }



}


bool CEKernelizer::applyHeuristicKernelization(bool debugRuleApplicationOnTheFly) {
    TimeMeasurer::startMeasurement("HeuristicKernelization");

    bool iter_changes = true, changes = false;
    while( iter_changes ) {
        iter_changes = false;

        if (rule_heur_1()) {
            if (debugRuleApplicationOnTheFly) clog << "Rule heur_1 applies" << endl;
            iter_changes = true;
        }
        if (!iter_changes && rule_heur_2()) {
            if (debugRuleApplicationOnTheFly) clog << "Rule heur_2 applies" << endl;
            iter_changes = true;
        }
        if (!iter_changes && iterative_rule_heur_3()) {
            if (debugRuleApplicationOnTheFly) clog << "Rule heur_3 applies" << endl;
            iter_changes = true;
        }
        if (!iter_changes && rule_heur_4()) {
            if (debugRuleApplicationOnTheFly) clog << "Rule heur_4 applies" << endl;
            iter_changes = true;
        }
        if (!iter_changes && rule_heur_5()) {
            if (debugRuleApplicationOnTheFly) clog << "Rule heur_5 applies" << endl;
            iter_changes = true;
        }
        if (!iter_changes && rule_heur_6()) {
            if (debugRuleApplicationOnTheFly) clog << "Rule heur_6 applies" << endl;
            iter_changes = true;
        }
        if (!iter_changes && rule_heur_7()) {
            if (debugRuleApplicationOnTheFly) clog << "Rule heur_7 applies" << endl;
            iter_changes = true;
        }

        if(iter_changes) changes = true;
        if( !use_heuristic_rules_separately ) break;
    }

    TimeMeasurer::stopMeasurement("HeuristicKernelization");
    return changes;
}



void CEKernelizer::test() {

    {
        ENDL(2);
        {
            clog << "Testing creating critical cliques" << endl;

            {
                clog << "\ttest: big_test_50_104" << endl;
                VVI V = CE_test_graphs::big_test_50_104;
                CEKernelizer kern(V);

                auto cc = kern.createCriticalCliques();
                sort(ALL(cc), [](auto &c1, auto &c2) { return c1.nodes.size() > c2.nodes.size(); });
                assert(cc.size() == 50 - 4);
            }

            ENDL(3);

            {
                clog << "\ttest: path P9" << endl;
                VVI V = CE_test_graphs::pathP9;
                CEKernelizer kern(V);

                auto cc = kern.createCriticalCliques();
                sort(ALL(cc), [](auto &c1, auto &c2) { return c1.nodes.size() > c2.nodes.size(); });
                assert(cc.size() == 9);
                clog << "P9 test passed" << endl;
            }

            clog << "Creating critical cliques passed" << endl;
        }

        ENDL(3);

        {
            clog << "Testing editing degrees, CE_test_graphs::kern_ed" << endl;
            VVI V = CE_test_graphs::kern_ed;

            CEKernelizer kern(V);

            VI K, N1;
            VB inK(V.size(), false), inN1(V.size(),false);

            {
                K = {0,1,2};
                VI N1 = kern.getNeighborhood(K,1);

                for( int d : K ) inK[d] = true; for( int d : N1 ) inN1[d] = true;
                assert(kern.getEditingDegree( 3, inK, inN1,N1.size() ) == 5);
                assert(kern.getEditingDegree( 4, inK, inN1, N1.size()) == 2);
                for( int d : K ) inK[d] = false; for( int d : N1 ) inN1[d] = false;
            }

            {
                K = {7,8};
                VI N1 = kern.getNeighborhood(K,1);

                for( int d : K ) inK[d] = true; for( int d : N1 ) inN1[d] = true;
                assert(kern.getEditingDegree( 3, inK, inN1, N1.size()) == 5);
                assert(kern.getEditingDegree( 5, inK, inN1, N1.size()) == 2);
                assert(kern.getEditingDegree( 6, inK, inN1, N1.size()) == 1);
                for( int d : K ) inK[d] = false; for( int d : N1 ) inN1[d] = false;
            }

            {
                K = {9};
                VI N1 = kern.getNeighborhood(K,1);

                for( int d : K ) inK[d] = true; for( int d : N1 ) inN1[d] = true;
                assert(kern.getEditingDegree( 3, inK, inN1, N1.size()) == 7);
                assert(kern.getEditingDegree( 4, inK, inN1, N1.size()) == 6);
                assert(kern.getEditingDegree( 5, inK, inN1, N1.size()) == 4);
                assert(kern.getEditingDegree( 6, inK, inN1, N1.size()) == 3);
                for( int d : K ) inK[d] = false; for( int d : N1 ) inN1[d] = false;
            }

            clog << "Editing degrees passed" << endl;
        }

        ENDL(3);

        {
            clog << "\tRule1: kern1" << endl;
            VVI V = CE_test_graphs::kern1;

            CEKernelizer kern(V);
            kern.rule1();

            assert( kern.inCl[0] == kern.inCl[1] );

            assert( kern.inCl[2] == kern.inCl[3] );

            assert( kern.inCl[4] == kern.inCl[5] );

            assert( kern.inCl[6] == kern.inCl[7] );
            assert( kern.inCl[7] == kern.inCl[8] );

            assert( kern.inCl[9] == 9 );
            assert( kern.inCl[10] == 10 );
            assert( kern.inCl[11] == 11 );
        }
    }

    ENDL(3);

    {
        clog << "Rule2: CE_test_graphs::kern2" << endl;

        {
            VVI V = CE_test_graphs::kern2;
            CEKernelizer kern(V);
            kern.rule2();

            assert(kern.inCl[0] == kern.inCl[1] &&
                   kern.inCl[0] == kern.inCl[2]
            );
        }

        {
            VVI V = CE_test_graphs::kern2;

            V[3].push_back(5);
            V.push_back( VI(1,3) ); // added node should prevent kernelization rule 2

            CEKernelizer kern(V);
            kern.rule2();

            assert(kern.inCl[0] == kern.inCl[1] );
            assert(kern.inCl[0] == kern.inCl[2] );

            assert(kern.inCl[3] == 3 );
            assert(kern.inCl[4] == 4 );
        }

        clog << "Rule2 passed" << endl;
    }

    ENDL(3);

    {

        {
            clog << "Rule3: CE_test_graphs::pathP3" << endl;
            VVI V = CE_test_graphs::pathP3;
            CEKernelizer kern(V);
            kern.rule3();

            for( int i=0; i<V.size(); i++ ) assert( kern.inCl[i] == i );
        }

        {
            clog << endl << endl << "Rule3: CE_test_graphs::kern3" << endl;
            VVI V = CE_test_graphs::kern3;
            CEKernelizer kern(V);
            kern.rule3();

            for( int i=0; i<V.size(); i++ ){
                if( i == 6 || i == 11 ) assert( V[i].size() == (kern.V[i].size() + 1) );
                else assert( V[i].size() == kern.V[i].size() );
            }
        }


        clog << "Rule3 passed" << endl;
    }

    ENDL(3);

    {
        {
            clog << endl << endl << "Rule4: CE_test_graphs::kern4positive" << endl;
            VVI V = CE_test_graphs::kern4positive;
            CEKernelizer kern(V);
            kern.rule4();

            for( int i=0; i<V.size(); i++ ) assert( kern.inCl[i] == kern.inCl[0] ); // just one clique
        }

        {
            clog << endl << endl << "Rule4: CE_test_graphs::kern4negative" << endl;
            VVI V = CE_test_graphs::kern4negative;
            CEKernelizer kern(V);
            kern.rule4();

            for( int i=0; i<V.size(); i++ ){
                assert( equal( ALL(V[i]), ALL(kern.V[i]) ) ); // nothing done by the rule
            }
        }
    }

    ENDL(3);

    {
        {
            clog << "Rule6_lemma3: CE_test_graphs::kern6positive" << endl;
            VVI V = CE_test_graphs::kern6positive;

            CEKernelizer kern(V);
            kern.rule6_lemma3();

            for( int i=0; i<9; i++ ) assert( kern.inCl[i] == kern.inCl[0] );
            assert( kern.inCl[9] == kern.inCl[10] );
            assert( kern.inCl[11] == kern.inCl[12] );

        }

        {
            clog << endl << endl << "Rule6_lemma3: CE_test_graphs::kern6negative" << endl;
            VVI V = CE_test_graphs::kern6negative;

            CEKernelizer kern(V);
            kern.rule6_lemma3();

            for( int i=0; i<V.size(); i++ ){
                assert( equal( ALL(V[i]), ALL(kern.V[i]) ) ); // nothing done by the rule
            }
        }
    }

    ENDL(3);

    {
        {
            clog << "Rule7, CE_test_graphs::kern7positive" << endl;
            VVI V = CE_test_graphs::kern7positive;
            CEKernelizer kern(V);
            kern.inV[7] = false; // node 7 is missing
            kern.rule7();

            for( int i=1; i<7; i++ ) assert( kern.inCl[i] == kern.inCl[i-1] ); // clique {0,1,2,3,4,5,6}
        }

        {
            clog << endl << endl << "Rule7, CE_test_graphs::kern7negative" << endl;
            VVI V = CE_test_graphs::kern7negative;
            CEKernelizer kern(V);
            kern.inV[7] = false; // node 7 is missing
            kern.rule7();

            for( int i=0; i<V.size(); i++ ){
                assert( equal( ALL(V[i]), ALL(kern.V[i]) ) ); // nothing done by the rule
            }
        }

        {
            clog << endl << endl << "Rule7, CE_test_graphs::kern8positive" << endl;
            VVI V = CE_test_graphs::kern8positive;
            CEKernelizer kern(V);
            kern.rule7();

            for( int i=0; i<V.size(); i++ ){
                assert( equal( ALL(V[i]), ALL(kern.V[i]) ) ); // nothing done by the rule
            }
        }
    }

    ENDL(3);

    {
        {
            clog << "Rule8, CE_test_graphs::kern8positive" << endl;
            VVI V = CE_test_graphs::kern8positive;
            CEKernelizer kern(V);
            kern.rule8();

            for( int i=0; i<V.size(); i++ ){
                if( i == 3 || i == 5 ) assert( kern.V[i].size() == (V[i].size()+1) );
                else assert( equal(ALL(V[i]), ALL(kern.V[i])) );
            }
        }

        {
            clog << endl << endl << "Rule8, CE_test_graphs::kern8positive2" << endl;
            VVI V = CE_test_graphs::kern8positive2;
            CEKernelizer kern(V);
            kern.inV[9] = false;
            kern.rule8();

            assert( find( ALL(kern.V[2]),8 ) == kern.V[2].end() );
            assert( find( ALL(kern.V[4]),7 ) != kern.V[4].end() );
            assert( find( ALL(kern.V[7]),4 ) != kern.V[7].end() );
            for( int d : {0,1,11,12,3,5,6,10} ) assert( equal( ALL(V[d]), ALL(kern.V[d]) ) );
        }

        {
            clog << endl << endl << "Rule8, CE_test_graphs::kern8negative" << endl;
            VVI V = CE_test_graphs::kern8negative;
            CEKernelizer kern(V);
            kern.inV[9] = false;
            kern.rule8();

            for( int i=0; i<V.size(); i++ ){
                assert( equal( ALL(V[i]), ALL(kern.V[i]) ) ); // nothing done by the rule
            }
        }
    }

    ENDL(3);

    {
        {
            clog << "Rule9, CE_test_graphs::kern9positive" << endl;
            VVI V = CE_test_graphs::kern9positive;
            CEKernelizer kern(V);
            kern.rule9();

            for( int i=0; i<5; i++ ) assert( kern.inCl[i] == kern.inCl[0] );

            kern.rule9(); // applying again - this will mark {3,5} as the same cluster

            for( int i=0; i<V.size(); i++ ) assert( kern.inCl[i] == kern.inCl[0] );

            assert( kern.V[0].empty() );
            assert( kern.V[1].empty() );

            int cnt = 0;
            for(int i=2; i<=4; i++) if( kern.V[i].empty() ) cnt++;
            assert(cnt == 2);
        }

        {
            clog << endl << endl << "Rule9, CE_test_graphs::kern9negative" << endl;
            VVI V = CE_test_graphs::kern9negative;
            CEKernelizer kern(V);
            kern.rule9();

            for( int i=0; i<V.size(); i++ ){
                assert( equal( ALL(V[i]), ALL(kern.V[i]) ) ); // nothing done by the rule
            }
        }
    }

    ENDL(3);

    {
        {
            clog << "Rule_heur_1: CE_test_graphs::kern_heur_1" << endl;
            VVI V = CE_test_graphs::kern_heur_1;
            CEKernelizer kern(V);

            kern.rule_heur_1();

            assert( kern.V[11].empty() );
            assert( kern.inCl[9] == kern.inCl[10] );
            assert( kern.V[5].empty() );
            assert( kern.inCl[3] == kern.inCl[2] );
            assert( kern.inCl[6] == kern.inCl[7] );
            assert( kern.inCl[6] == kern.inCl[8] );
            assert( kern.inCl[0] == kern.inCl[1] );
            assert( kern.inCl[4] != kern.inCl[3] );
            assert( kern.inCl[4] != kern.inCl[2] );

            kern.rule_heur_1(); // running again

            assert( kern.inCl[4] == kern.inCl[3] );
            assert( kern.inCl[4] == kern.inCl[2] );
        }
    }

    ENDL(3);

    {
        {
            clog << "Rule_heur_2: CE_test_graphs::kern_heur_2" << endl;
            VVI V = CE_test_graphs::kern_heur_2;
            CEKernelizer kern(V);
            kern.rule_heur_2();

            assert(kern.V[4].empty());
            assert(kern.V[5].empty());
        }
    }

    clog << "CEKernelizer tested" << endl;
}

bool CEKernelizer::fillMissingEdgesBetweenClusterNodes() {
    const bool debug = false;
    // Filling all missing edges between nodes marked to the same cluster
    vector< unordered_set<int> > edges(N);
    for( int i=0; i<N; i++ ) edges[i].insert(ALL(V[i]));

    int edges_added = 0;
    VVI clusters(N);
    for( int d =0; d<N; d++ ) if( !inV[d] && !V[d].empty() ) clusters[ fau.Find(d) ].push_back(d);
    for( auto& cl : clusters ){
        if( cl.size() <= 1 ) continue;
        for( int i=0; i<cl.size(); i++ ){
            int a = cl[i];
            for( int j=i+1; j<cl.size(); j++ ){
                int b = cl[j];
                if( edges[a].count(b) == 0 ){ // edge (a,b) is not in the graph
                    GraphUtils::addEdge( V,a,b );
                    edges_added++;
                    if(debug)
                        clog << "In rule15, adding edge (" << a << "," << b << ")" << endl;
//                        edges[a].insert(b); edges[b].insert(a); // no need in fact to add nodes a and b to [edges]
                }
            }
        }
    }

    return edges_added;
}















