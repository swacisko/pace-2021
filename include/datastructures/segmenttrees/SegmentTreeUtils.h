//
// Created by sylwester on 1/9/21.
//

#ifndef ALGORITHMSPROJECT_SEGMENTTREEUTILS_H
#define ALGORITHMSPROJECT_SEGMENTTREEUTILS_H

#include "Makros.h"

namespace SegmentTreeUtils{

    /**
     * Calculates overlapping segment between two segments [a.first,a.second] and [x.first, x.second].
     * The resulting segment [x,y] may have x > y. If so, then there is no overlapping segment.
     * @param a
     * @param x
     * @return
     */
    template<class _T>
    extern const pair<_T,_T> getOverlapSegment( pair<_T,_T> a, pair<_T,_T> x ) {
        return { max(a.first, x.first), min(a.second, x.second) };
    }

    /**
     * Calculates and return number of common elements (integers) in given two segments
     * [a.first,a.second] and [x.first, x.second].
     * @tparam _T
     * @param a
     * @param x
     * @return
     */
    template<class _T>
    extern const _T getOverlapLength( pair<_T,_T> a, pair<_T,_T> x ){
        auto s = getOverlapSegment(a,x);
        if(s.first <= s.second) return s.second - s.first + 1;
        else return 0;
    }
}

#endif //ALGORITHMSPROJECT_SEGMENTTREEUTILS_H
