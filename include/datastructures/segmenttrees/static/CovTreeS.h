//
// Created by sylwester on 1/9/21.
//

#ifndef ALGORITHMSPROJECT_COVTREE_H
#define ALGORITHMSPROJECT_COVTREE_H

#include "Makros.h"
#include "datastructures/segmenttrees/SegmentTreeUtils.h"
#include "utils/RandomNumberGenerators.h"

/**
 * CovTreeS is a static segment tree that enables in O(logN) to do the following operation:
 * - cover all elements on an interval [x,y] with a segment with given id  (in O(logN) )
 * - find number of segments that cover given point x
 * - find all ids of segments that cover given point x (in O(logN + res.size()) )
 * @tparam _T
 */
class CovTreeS{
    static const int L(int i){ return i<<1; } // id of left son of i
    const int R(int i){ return (i<<1)+1; } // id of right son of i
    const int PAR(int i){ return i>>1; } // if od parent of i

public:
    /**
     * Creates segment tree with maximal size S.
     * @param S
     */
    CovTreeS(int S);

    /**
     * Adds a segment [x,y] with given [id] to the set of segments.
     */
    void add( int x, int y, int id ){ addHelper(x,y,id,1,0,N-1); }

    /**
     * Finds and returns number of segments covering totally point x
     */
    int count( int x );

    /**
     * Finds and returns ids of all segments covering totally point x
     * Works in logN + K, where K is the number of found segment id's.
     */
    VI find( int x );

    /**
     * Runs tests.
     */
    static void test();

private:
    /**
     * Elements that can be accessed are from range [0,N-1]. N is the maximal size given in constructor rounded up
     * to the nearest power of two.
     */
    int N;

    VVI segments; // covered[i] is the sum of segments covering i-th-node-interval totally.

    /**
     * 'Covers' all elements from interval [x,y] \inters [a,b] with segment with given [ids].
     * @param ids id of the segment added
     * @param idn id of the currently processed node
     */
    void addHelper( int x, int y, int ids, int idn, int a, int b );
};

#endif //ALGORITHMSPROJECT_COVTREE_H
