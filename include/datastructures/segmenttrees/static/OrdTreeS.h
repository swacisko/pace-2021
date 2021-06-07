//
// Created by sylwester on 1/9/21.
//

#ifndef ALGORITHMSPROJECT_ORDTREES_H
#define ALGORITHMSPROJECT_ORDTREES_H

#include <utils/TimeMeasurer.h>
#include "SumTreeS.h"

/**
 * OrdTreeS is a static segment tree that enables in O(logN) to do the following operation:
 * - add [val] occurences of each element from interval [x,y] to the set
 * - count the k-th smallest element in the set.
 *
 * To add elements use function add from SumTreeS.
 *
 * Use _T = int if the query k and number of elements in the set will be always less than 2^32.
 * Otherwise use _T = long long
 * @tparam _T
 */
template<class _T>
class OrdTreeS : public SumTreeS<_T>{
    static const int L(int i){ return i<<1; } // id of left son of i
    const int R(int i){ return (i<<1)+1; } // id of right son of i
    const int PAR(int i){ return i>>1; } // id of parent of i, or just half of the given value
public:
    /**
     * Creates a segment tree with elements from range [0,S-1].
     * @param S
     */
    OrdTreeS(int S) : SumTreeS<_T>(S){}

    /**
     * Finds and returns the K-th statistic, that is K-th smallest element in the set.
     */
    _T findKthStatistic(_T K){ return findKthStatisticHelper(K, 0, this->N-1, 1,0); }

    /**
     * Runs tests
     */
    static void test();

protected:

    /**
     * Finds and returns K-th statistic. Each element from the segment [a,b] occurs factor times in the 'upper nodes'.
     * @param a beginning of the interval of the current node
     * @param b end of the interval of the current node
     * @param factor since values could be added on segments, not only points, this is necessary to keep track of
     * their numbers
     */
    int findKthStatisticHelper( _T K, int a, int b, int id, _T factor );

};

template<class _T>
int OrdTreeS<_T>::findKthStatisticHelper( _T K, int a, int b, int id, _T factor ){
    if(id == 1 && K > this->sums[id]) return -1;
    if(K == 0 ) return -1;
    if( a == b ) return a;
    int l = L(id);
    int ovla = PAR(b-a+1);
    factor += this->covered[id];

    if( K > this->sums[l] + factor * ovla ){
        return findKthStatisticHelper( K - this->sums[l] - factor * ovla, PAR(a+b)+1,b,R(id), factor );
    }else{
        return findKthStatisticHelper( K, a,PAR(a+b),l, factor );
    }

}

template<class _T>
void OrdTreeS<_T>::test(){

    int N = 1'000'000;
    int T = 50'000;
    VLL V(N,0);
    auto gen = UniformIntGenerator(0,1ll * 1e6 * 1e6);
    OrdTreeS<long long> otree(N);

    long long sumV = 0;
    long long sumCheck = 0;
    long long sumFind = 0;
    for(int testid=0; testid < T; testid++){
        { // modification
            int a = gen.rand() % N;
            int b = a + gen.rand() % (N-a);
            if(a>b) swap(a,b);
            if (testid % 100 == 0) a = b;


            int v = gen.rand() % 1000;

            TimeMeasurer::startMeasurement("brute_add");
            for(int i=a; i<=b; i++) V[i] += v;
            sumV += v * (b-a+1);
            sumCheck += b-a+1;
            TimeMeasurer::stopMeasurement("brute_add");


            TimeMeasurer::startMeasurement("tree_add");
            otree.add( a,b,v );
            TimeMeasurer::stopMeasurement("tree_add");

        }

        { // query
            long long K = gen.rand() % sumV;
            if (testid % 100 == 0) K <<= 1;

            TimeMeasurer::startMeasurement("brute_find");
            long long kk = K;
            int res1 = -1;
            for( int i=0; i<N; i++ ){
                if( kk <= V[i] ){
                    res1 = i;
                    sumFind += i;
                    break;
                }
                kk -= V[i];
            }
            if(K == 0) res1 = -1;
            TimeMeasurer::stopMeasurement("brute_find");

            TimeMeasurer::startMeasurement("tree_find");
            int res2 = otree.findKthStatistic(K);
            TimeMeasurer::stopMeasurement("tree_find");

            if( res1 != res2 ){
                DEBUG2(testid,K);
                cerr << "V:" << endl;
                for(int i=0; i<N; i++) if( V[i]>0 ) cerr << "(" << i << "," << V[i] << "), ";
                cerr << endl;
                DEBUG2(res1,res2);
                assert(res1 == res2);
            }

        }

    }

    TimeMeasurer::writeAllMeasurements();
    DEBUG(sumCheck);
    DEBUG(sumFind);
    DEBUG(sumV);
    clog << "TESTS PASSED" << endl;
}

#endif //ALGORITHMSPROJECT_ORDTREES_H
