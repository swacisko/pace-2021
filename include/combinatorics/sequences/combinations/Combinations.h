//
// Created by sylwester on 5/11/20.
//

#ifndef ALGORITHMSPROJECT_COMBINATIONS_H
#define ALGORITHMSPROJECT_COMBINATIONS_H

#include "Makros.h"

namespace Combinations{

    /**
     * Rearranges vector a to the next k-element combination of universe [1,...,n], where
     * k = a.size()
     * @param a
     * @param n
     * @return true, if next combination was found, false otherwise
     */
    extern bool next_combination(vector<int>& a, int n);
}

#endif //ALGORITHMSPROJECT_COMBINATIONS_H
