//
// Created by sylwester on 9/6/19.
//

#ifndef ALGORITHMSPROJECT_STANDARDUTILS_H
#define ALGORITHMSPROJECT_STANDARDUTILS_H

#include "Makros.h"
#include "RandomNumberGenerators.h"

namespace StandardUtils{

    /**
     * this function returns smallest integer d int [p,q] such that fun returns true for d as argument or -1 if such d does not exist.
     * in other words if the table has form {0,0,0,0,0,0,0,0,1,1,1,1,1,1} bijary search will return index of the first true element
     */
    extern int binarySearch( int p, int q, function<bool (int)> objComp );


    /**
     *
     * @param tab
     * @return vector containing indices for which tab[i] is true
     */
    extern VI extractIndices(VB & tab);


    extern VI toVI(VB& v);

    /**
     *
     * @param N
     * @param v
     * @return VB of soze N in which v[i] is set for all i in {0,v.size()-1}
     */
    extern VB toVB(int N, VI& v);

    /**
     * Given a list partition. creates layers such that res[partition[i]] contains i.
     * @param partition
     * @return
     */
    extern VVI partitionToLayers( VI & partition );

    /**
     * Given layers, creates a prtition such that i is in layer partition[i].
     */
    extern VI layersToPartition( VVI & layers );

    /**
     * Binary power x^a.
     * @tparam _T
     * @param x
     * @param a this need to be at least 1
     * @param multiply this function is used to multiply two objects
     * @param square this function is used to square two objects. By default it is multiplication, but in some cases squaring may be faster than usual product.
     * @return
     */
    template<class _T>
    _T binaryPower( _T x, int a, function< _T(_T,_T) > multiply = [](_T x, _T y ){ return x*y; }, function< _T(_T) > square = [](_T x){ return x*x; } ){
        if( a == 1 ) return x;
        else if( a & 1 ){

            _T t = binaryPower(x, (a-1)>>1, multiply, square );
            return multiply( x, square(t) );
        }else{
            _T t = binaryPower(x, a>>1, multiply, square );
            return square(t);
        }
    }

    /**
     * Writes given integer number in binary format (right bit is least significant).
     * @tparam _T
     * @param x
     */
    template<class _T>
    void writeInBinary( _T x, int bits = -1 ){
        int i = 0;
        while( (1<<(i+1)) <= x ) i++;

        for( int j=0; j<min(bits,bits - i-1); j++ ) cerr << 0;
        while( i >= 0 ){
            cerr << __builtin_popcount( (x & (1<<i)) );
            i--;
        }
        cerr << endl;
    }

    /**
     * Creates and returns a subarray of array [s] with beginning in index beg, and end in index end
     */
    template<class _T>
    vector<_T> getSubarray(vector<_T> s, int beg, int end){
        vector<_T> res;
        while(beg <= end){
            res.push_back(s[beg]);
            beg++;
        }
        return res;
    }


    /**
     * Implementation of function similar to that of python split
     */
    vector<string> split( string & s, string pat = " " );

    /**
     * The same as v1.insert( v1.end(), ALL(v2) );
     * @tparam _T
     * @param v1
     * @param v2
     */
    template<class _T>
    void append( vector<_T> & v1, vector<_T> & v2 ){ v1.insert( v1.end(), ALL(v2) ); }

    /**
     * Equivalent to python slicing. Extracts and returns in a vector all elements from [a,b) that are at positions
     * a, a+step, a+2*step, ...
     * It is also possible to pass a > b and step < 0
     */
    template<class _T>
    vector<_T> slice( vector<_T> & v, int a, int b, int step ){
        cerr << "SLICE not tested yet" << endl;
        assert(step != 0);
        vector<_T> res;
        res.reserve( abs( (b-a) / step ) );
        if(a<=b){
            assert(step > 0);
            while(a<b){
                res.push_back( a );
                a += step;
            }
        }
        else{
            assert(step < 0);
            while(a>b){
                res.push_back(a);
                a += step;
            }
        }
        return res;
    }

    template<class _T>
    void shuffle( vector<_T> & V ){
        UniformIntGenerator rnd(0,1e9);
        shuffle(V, rnd.getRNG());
    }

    template<class _T, class _rnd>
    void shuffle( vector<_T> & V, _rnd rnd ){
        std::uniform_int_distribution<long long> unif( 0, 10ll * V.size() );

        for( int i=(int)V.size()-1; i>=0; i-- ){
            int ind = unif(rnd) % (i+1);
            if( ind != i ) swap( V[i], V[ind] );
        }
    }

    /**
     * Creates a vector of pairs from two vectors. E.g. from [0,4,7], ['a', 's'] it will make [ {0,'a'}, {4,'s'} ]
     */
    template<class _t, class _s>
    vector< pair<_t,_s> > zip( vector<_t> a, vector<_s> b ){
        vector<pair<_t,_s>> res(min(a.size(), b.size()));
        for( int i=0; i<min(a.size(), b.size()); i++ ) res[i] = {a[i], b[i]};
        return res;
    }

    /**
     * Reverse to the [zip] function - creates two vectors from a vector of pairs.
     */
    template<class _t, class _s>
    pair<vector<_t>, vector<_s> > unzip( vector< pair<_t,_s> > &&a ){
        pair<vector<_t>, vector<_s> > res;
        res.first.reserve(a.size());
        res.second.reserve(a.size());
        for( int i=0; i<a.size(); i++ ){
            res.first.push_back(a[i].first);
            res.second.push_back(a[i].second);
        }
        return res;
    }

}

#endif //ALGORITHMSPROJECT_STANDARDUTILS_H
