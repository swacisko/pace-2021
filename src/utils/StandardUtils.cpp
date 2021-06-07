//
// Created by sylwester on 9/6/19.
//

#include <utils/StandardUtils.h>

#include "utils/StandardUtils.h"


namespace StandardUtils{

    int binarySearch( int p, int q, function<bool (int)> objComp ){
        int P = p;
        int Q = q+1;
        int d;

        while( P < Q ){
            d = ( (P+Q) >> 1 );
            if( objComp(d) ) Q = d;
            else P = d+1;
        }

        if(Q > q) return -1;
        else return Q;
    }

    VI extractIndices(VB &tab) {
        return toVI(tab);
      /*  VI res;
        for(int i=0; i<tab.size(); i++){
            if(tab[i]) res.push_back(i);
        }
        return res;*/
    }

    VI toVI(VB &v) {
        VI res; res.reserve( max( 5, (int)v.size() / 2 ) );
        for(int i=0; i <v.size(); i++) if( v[i] ) res.push_back(i);
        return res;
    }

    VB toVB(int N, VI &v) {
        VB res(N,false);
        for(int d : v) res[d] = true;
        return res;
    }

    vector<string> split(string &s, string pat) {
        cerr << "not implemented yet" << endl;
        exit(18);
    }

    VVI partitionToLayers(VI &partition) {
        int M = *max_element(ALL(partition));
        VVI v(M+1);
        for( int i=0; i<partition.size(); i++ ) if( partition[i] != -1 ) v[partition[i]].push_back(i);
        return v;
    }

    VI layersToPartition(VVI &layers) {
        int M = -1;
        for( int i=0; i<layers.size(); i++ ) M = max( M, *max_element( ALL(layers[i]) ) );

        VI part(M+1,-1);
        for( int i=0; i<layers.size(); i++ ) for(int t : layers[i]) part[t] = i;
        return part;
    }


}