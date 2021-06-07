//
// Created by sylwester on 3/9/21.
//

#include <Constants.h>
#include "CONTESTS/PACE21/heur/EOCreators/ComponentExpansion.h"



ComponentExpansion::ComponentExpansion( Cluster & cl ){
    this->cl = &cl;
    this->clg = &cl.g;
    V = &clg->V;
    N = V->size();
    eInS = VI(N, 0);

    sumEW = VI(N, 0);
    for( int i=0; i<N; i++ ) for( auto & [d,w] : (*V)[i] ) sumEW[i] += w;

    inS = VB(N,false);
    sumNWinS = 0;

    cmp_rules = { 3,2,1 }; // by default rules are executed in order 3,2,1

    {
        heap = Heap<int>( N,0, [&](int a, int b ){ return cmpFun(a, b); } );
        for( int i=0; i<N; i++ ){
            heap.removeFromHeap(i); // removing all elements from heap - only those that will be necessary will be
            // pushed again
        }
    }
}

ExpansionOrder ComponentExpansion::getExpansionOrder(VI A, const bool use_heap, const int max_nodes_in_eo) {
    const bool debug = false;

    S = A;
//    fill(ALL(inS),false);
//    fill(ALL(eInS), 0);
//
//    fill(ALL(sumEW), 0);
//    for( int i=0; i<N; i++ ) for( auto & [d,w] : (*V)[i] ) sumEW[i] += w;

    cut_value = 0;
    sumNWinS = 0;

    for( int d : S ){
        inS[d] = true;
        sumNWinS += clg->node_weights[d];
    }
    for( int d : S ){
        for( auto & [a,w] : (*V)[d] ){
            eInS[a] += w;
            if( !inS[a] ){
                cut_value += w;
            }
        }
    }

    if(debug){
        DEBUG(sumEW);
        DEBUG(S); DEBUG(sumNWinS); DEBUG(eInS); DEBUG(cut_value);
    }

    if(use_heap){
//        heap = Heap<int>( N,0, [&](int a, int b ){ return cmpFun(a, b); } );
//        for( int i=0; i<N; i++ ){
//            if(inS[i] || eInS[i] == 0) heap.removeFromHeap(i);
//            else heap.set( i, i );
//        }

        for( int p : S ){
//            heap.set(p,p);
            for( auto & [d,w] : (*V)[p] ){
                if( inS[d] ) continue;
                heap.set(d,d);
            }
        }
    }


    while( (!use_heap || !heap.empty() ) && S.size() < min(N,max_nodes_in_eo) ){
        int d = 0;
        if(use_heap) d = heap.extract_min();
        else{
            d = -1;
            for( int i=0; i<N; i++ ){
                if( !inS[i] ){
                    if( d == -1 ) d = i;
                    else if( eInS[i] > 0 && cmpFun( i,d ) ) d = i; // we consider only nodes that are neighbors to S
                }
            }
        }

        if( terminate_on_cluster_violator && eInS[d] < ( 3 * sumNWinS / 4)  ) break;

        moveNodeToS(d, use_heap);

        if(debug){
            DEBUG(S); DEBUG(sumNWinS); DEBUG(eInS); DEBUG(cut_value);
            clog << endl;
        }
    }

    for( int d : S ){
        inS[d] = false; // clearing
        eInS[d] = 0;
        heap.removeFromHeap(d);

        for( auto & [p,w] :(*V)[d] ){
            eInS[p] = 0;
            heap.removeFromHeap(p);
        }
    }

    return ExpansionOrder(S,cl);
}

void ComponentExpansion::moveNodeToS(int d, bool use_heap) {
    const bool debug = false;
    if(debug) clog << "Adding node " << d << endl;

    S.push_back(d);
    inS[d] = true;
    sumNWinS += clg->node_weights[d];

    for( auto & [a,w] : (*V)[d] ){
        if( !inS[a] ){
            eInS[a] += w;
            cut_value += w;
           if( use_heap ) heap.set( a,a ); // value of node is just its id, because we use comparator base on its id
        }
        else cut_value -= w;
    }
}

bool ComponentExpansion::cmpFun(int a, int b) {
    for( int i=0; i<cmp_rules.size(); i++ ){
        int rule = cmp_rules[i];
        switch(rule){
            case 1:{
                // 1. First rule: select node with tightest (sum of edge weights) connection to S.
                int val1 = eInS[a];
                int val2 = eInS[b];
                if( val1 != val2 ) return val1 > val2;
                break;
            }
            case 2:{
                // 2. Second rule: select node with least weight of edges to  neighbors outside S.
                int val1 = (sumEW[a] - eInS[a]);
                int val2 = (sumEW[b] - eInS[b]);
                if(val1 != val2 ) return val1 < val2;
                break;
            }
            case 3:{
                //3. Third rule: select node that minimizes X - Y, where X is weight of edges to neighbors outside S
                // and Y is sum of weights to neighbors in S. Because X = deg - Y, we have X-Y = deg - 2*Y
                int val1 = sumEW[a] - ( eInS[a] << 1 );
                int val2 = sumEW[b] - ( eInS[b] << 1 );
                if(val1 != val2 ) return val1 < val2;
                break;
            }
            case 4:{
                // 4. Fourth rule: select node with AVERAGE tightest connection to S.
                int val1 = eInS[a] * clg->node_weights[b];
                int val2 = eInS[b] * clg->node_weights[a];
                if( val1 != val2 ) return val1 > val2;
                break;
            }
            case 5:{
                // 5. Fifth rule: select node with least AVERAGE number of neighbors outside S.
                int val1 = (sumEW[a] - eInS[a]) * clg->node_weights[b];
                int val2 = (sumEW[b] - eInS[b]) * clg->node_weights[a];
                if(val1 != val2 ) return val1 < val2;
                break;
            }
            case 6:{
                //6. Sixth rule: select node that minimizes AVERAGE  X - Y, where X is number of neighbors outside S and Y is number
                //  of neighbors in S. Because X = deg - Y, we have X-Y = deg - 2*Y
                int val1 = (sumEW[a] - ( eInS[a] << 1 )) * clg->node_weights[b];
                int val2 = (sumEW[b] - ( eInS[b] << 1 )) * clg->node_weights[a];
                if(val1 != val2 ) return val1 < val2;
                break;
            }
            case 10:{
                // avoids division and doubles
                LL cut1 = cut_value + (sumEW[a] - ( eInS[a] << 1 ));
                LL cut2 = cut_value + (sumEW[b] - ( eInS[b] << 1 ));
                int w1 = clg->node_weights[a];
                int w2 = clg->node_weights[b];
                cut1 *= ((sumNWinS + w2) * (cl->cluster_weight - sumNWinS - w2));
                cut2 *= ((sumNWinS + w1) * (cl->cluster_weight - sumNWinS - w1));
                if (cut1 != cut2) return cut1 < cut2;
                break;
            }
            case 11:{
                int cut1 = cut_value + (sumEW[a] - ( eInS[a] << 1 ));
                int cut2 = cut_value + (sumEW[b] - ( eInS[b] << 1 ));
                int w1 = clg->node_weights[a];
                int w2 = clg->node_weights[b];
                int all_edges_1 = ((sumNWinS + w1) * (cl->cluster_weight - sumNWinS - w1));
                int all_edges_2 = ((sumNWinS + w2) * (cl->cluster_weight - sumNWinS - w2));
                int diff1 = (cut1 << 1) - all_edges_1;
                int diff2 = (cut2 << 1) - all_edges_2;
                if( diff1 != diff2 ) return diff1 < diff2;
                break;
            }
            default:{ return a < b; }
        }
    }

    return a < b;
}


