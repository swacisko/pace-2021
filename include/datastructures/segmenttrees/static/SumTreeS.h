//
// Created by sylwester on 1/9/21.
//

#ifndef ALGORITHMSPROJECT_SUMTREES_H
#define ALGORITHMSPROJECT_SUMTREES_H

#include "Makros.h"
#include "datastructures/segmenttrees/SegmentTreeUtils.h"
#include "utils/RandomNumberGenerators.h"

/**
 * This is a segment tree. It enalbes in O(logN) time to:
 * - find sum of elements on given interval
 * - add any value to all elements in a given interval
 * @tparam _T type of elements, usually integer
 */
template<class _T>
class SumTreeS{

    using VT = vector<_T>;
    static const int L(int i){ return i<<1; } // id of left son of i
    const int R(int i){ return (i<<1)+1; } // id of right son of i
    const int PAR(int i){ return i>>1; } // if od parent of i

public:
    /**
     * Creates segment tree with maximal size S.
     * @param S
     */
    SumTreeS(int S);

    /**
     * Adds given value [val] to all elements in interval [x,y]
     */
    void add( int x, int y, _T val ){ addHelper( x,y,1,0,N-1,val ); }

    /**
     * Finds and returns sum of elements in given interval [x,y]
     */
    _T sum( int x, int y ){ return sumHelper( x,y, 1,0, N-1 ); }

    /**
     * Runs tests.
     */
    static void test();

protected:
    /**
     * Elements that can be accessed are from range [0,N-1]. N is the maximal size given in constructor rounded up
     * to the nearest power of two.
     */
    int N;

    VT covered; // covered[i] is the sum of segments covering i-th-node-interval totally.
    VT sums; // sums[i] is the sum on interval represented by i-th node

    /**
     * Adds given value [val] to all elements in the overlapping of [x,y] and [a,b].
     */
    void addHelper(int x, int y, int id, int a, int b, _T& val);

    /**
     * Finds and returns sum of all elements in the overlapping of [x,y] and [a,b].
     */
    _T sumHelper(int x, int y, int id, int a, int b);
};

template<class _T>
void SumTreeS<_T>::addHelper(int x, int y, int id, int a, int b, _T& val){
    PII ovl = SegmentTreeUtils::getOverlapSegment(PII(x,y), {a,b});
    if( ovl.second < ovl.first ) return;    // if [x,y] and [a,b] have no common elements
    if( ovl.second - ovl.first == b-a ){ covered[id] += val; sums[id] += val * (b-a+1); return; } // if [a,b] is entirely in [x,y]
    addHelper( x,y, L(id), a, PAR(a+b), val ); // PAR(a+b) is just (a+b)/
    addHelper( x,y, R(id), PAR(a+b)+1, b, val );
    sums[id] = sums[L(id)] + sums[R(id)] + covered[id] * (b-a+1);
}

template<class _T>
_T SumTreeS<_T>::sumHelper(int x, int y, int id, int a, int b){
    PII ovl = SegmentTreeUtils::getOverlapSegment(PII(x,y), {a,b});
    if( ovl.second < ovl.first ) return 0;  // if [x,y] and [a,b] have no common elements
    if( ovl.second - ovl.first == b-a ){ return sums[id]; } // if [a,b] is entirely in [x,y]
    return sumHelper( x,y, L(id), a, PAR(a+b) ) + sumHelper( x,y, R(id), PAR(a+b)+1, b ) +
        covered[id] * (ovl.second - ovl.first+1); // PAR(a+b) is just (a+b)/
}

template<class _T>
SumTreeS<_T>::SumTreeS(int S) : N(1){
    while( N < S ) N <<= 1;
    sums = VT(N<<1,0);
    covered = VT(N<<1,0);
}

template<class _T>
void SumTreeS<_T>::test(){

    int N = 1'000;

    VLL V(N,0);
    SumTreeS<long long> stree(N);

    int T = 10'000;
    auto gen = UniformIntGenerator(0,1e9);

    while(T--){
        { // modification
            int a = gen.rand() % N;
            int b = gen.rand() % N;
            if (a > b) swap(a, b);
            if (T % 100 == 0) a = b;

            int v = gen.rand();
            for (int i = a; i <= b; i++) V[i] += v;

            stree.add( a,b,v );
        }

        { // query
            int a = gen.rand() % N;
            int b = gen.rand() % N;
            if (a > b) swap(a, b);
            if (T % 100 == 0) a = b;

            long long S = accumulate( V.begin() + a, V.begin() + b+1,0ll );
            long long S2 = stree.sum(a,b);

            if(S != S2){
                DEBUG(T);
                DEBUG2(S,S2);
                DEBUG3(V,a,b);
                assert(S == S2);
            }

        }

    }

    clog << "TESTS PASSED" << endl;
}


#endif //ALGORITHMSPROJECT_SUMTREES_H
