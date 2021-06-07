//
// Created by sylwester on 9/24/19.
//

#ifndef ALGORITHMSPROJECT_SORTING_H
#define ALGORITHMSPROJECT_SORTING_H

#include "Makros.h"

namespace Sorting{

    /**
     * Count sort. Works in time O(W + N), where W is difference between largest and smallest element and N is the number of elements.
     * @tparam _T
     * @param V
     * @return
     */
    template<class _T>
    VI countSort( vector<_T> & V ){

        int n = *min_element(ALL(V));
        int N = *max_element(ALL(V));

        VI tab(n+N+2,0);
        for(int d : V) tab[d-n]++;
        VI res;
        res.reserve(V.size());
        for( int i=0; i<= n+N; i++ ){
            while( tab[i] > 0 ){
                res.push_back( i + n );
                tab[i]--;
            }
        }
        return res;
    }



}

#endif //ALGORITHMSPROJECT_SORTING_H
