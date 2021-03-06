//
// Created by sylwester on 3/25/20.
//

#include <graphs/cliques/CliqueUtils.h>
#include <graphs/GraphUtils.h>

#include "graphs/cliques/CliqueUtils.h"


#include "Makros.h"

namespace CliqueUtils{


    bool isClique( VVI& V, VI& clq){
        unordered_set<int> zb(ALL(clq));
        for( int p : clq ){
            int t = 0;
            for(int d : V[p]){
                if( zb.count(d) ) t++;
            }
            if( t != clq.size()-1 ) return false;
        }
        return true;
    }

    bool isClique(VVI &V, VI &clq, VB &helper) {
        for(int d : clq) helper[d] = true;
        for( int p : clq ){
            int t = 0;
            for(int d : V[p]){
                if( helper[d] ) t++;
            }
            if( t != clq.size()-1 ){
                for(int d : clq) helper[d] = false;
                return false;
            }
        }
        for(int d : clq) helper[d] = false;
        return true;
    }

    void fillToClique(VVI &V, VI nodes) {
        vector< unordered_set<int> > presEdges(nodes.size());
        for( int i=0; i<nodes.size(); i++ ){
            int p = nodes[i];
            for( int d : V[p] ) if( d > p ) presEdges[i].insert(d);
        }

        for( int i=0; i<nodes.size(); i++ ){
            int p = nodes[i];
            for( int d : nodes ){
                if( d > p && presEdges[i].count(d) == 0 ) GraphUtils::addEdge( V,p,d );
            }
        }

        assert( isClique( V, nodes ) );
    }




}
