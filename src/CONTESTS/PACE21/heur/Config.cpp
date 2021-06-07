//
// Created by sylwester on 3/8/21.
//

#include <CONTESTS/PACE21/heur/Config.h>
#include <combinatorics/CombinatoricUtils.h>
#include <graphs/GraphUtils.h>

Config::Config() {
    createConfiguration();
}

void Config::createConfiguration() {
    { // default configuration uses all swap candidate creators
        swpCndCreatorsToUse = {node, edge_all, exp_ord, triangle, exp_ord_attr, exp_ord_rep};
    }
}

VVI Config::expansionOrderInitialSetProvider(Cluster *cl) {

    { // there are sqrt(cl->size()) initial sets, each set is a random edge from cluster cl
        int T = 4 * ceil(sqrt(cl->size())); // #TEST - many expansion orders considered

        assert(T > 0);
        if( T == 1 ) return {{0}};

        VVI res;
        auto edges = GraphUtils::getGraphEdges( cl->g.V );
        int E =  (int)edges.size();
        VI subs;
        if(E>0) subs = CombinatoricUtils::getRandomSubset( E-1, min( E, T )  ); // DEFAULT SEED'ING USED

        for( int i=0; i<min(T,(int)subs.size()); i++ ){
            int e = subs[i];
            res.push_back(VI());
            res.back().push_back( get<0>(edges[e]) );
            res.back().push_back( get<1>(edges[e]) );
            if( i&1 ){
                // edge (a,b) or (b,a) given in a sort of 'random' order
                swap( res.back()[0], res.back()[1] );
            }
        }

        return res;
    }
}

void Config::setSpeedMode(SPEED_MODE mode) {
    speed_mode = mode;

    if( mode <= very_fast ){
        use_only_common_neighbors_in_swp_cnd_edge = true;
        for(int i=0; i<swpCndCreatorsToUse.size(); i++){ // removing SwpCndTriangle
            if( swpCndCreatorsToUse[i] == triangle ){
                swpCndCreatorsToUse.erase(swpCndCreatorsToUse.begin()+i);
                break;
            }
        }
    }

    if( mode <= fast ){
        use_only_empty_cluster_in_swp_cnd_triangle = true;
    }

    if( mode <= medium_fast ){
        for(int i=0; i<swpCndCreatorsToUse.size(); i++){ // removing ComponentExpansionAttraction
            if( swpCndCreatorsToUse[i] == exp_ord_attr ){
                swpCndCreatorsToUse.erase(swpCndCreatorsToUse.begin()+i);
                break;
            }
        }
    }
}

