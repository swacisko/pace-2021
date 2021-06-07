//
// Created by sylwester on 3/20/21.
//

#include <Constants.h>
#include <utils/RandomNumberGenerators.h>
#include "CONTESTS/PACE21/heur/SwapCandidates/SwpCndNode.h"



vector<SwpCndNode> SwpCndNodeCreator::createSwapCandidates() {
    vector<SwpCndNode> res;

    for( auto & cl : state->clusters ){
        auto cnds = create_MoveTo_SwapCandidatesForCluster(cl);
        res += cnds;
    }

    return res;
}

vector<SwapCandidate *> SwpCndNodeCreator::createSwapCandidatesRaw() {
    return vector<SwapCandidate *>();
}

vector<SwpCndNode> SwpCndNodeCreator::create_MoveTo_SwapCandidatesForCluster(Cluster &cl) {
    const bool debug = false;

    vector<SwpCndNode> res;

    for( int d : cl.g.nodes ){
        int in_cl_d = state->inCl[d];
        int deg_in_cl = state->degInCl[d];
        int nw = clg->node_weights[d];
        int clw = cl.cluster_weight;
        int swpval1 = nw * (clw - nw) - deg_in_cl; // number of edges that had to be added to keep d in its cluster

        if(debug){
            clog << "Considering node d: " << d << endl;
            DEBUG(in_cl_d);DEBUG(deg_in_cl);DEBUG(nw);DEBUG(clw);DEBUG(swpval1);
        }

        VI cluster_neighbors( 1, state->getIdOfEmptyCluster() );

        updateEdgesInNeighboringClustersForNode( d, in_cl_d, cluster_neighbors );

        if(debug) DEBUG(cluster_neighbors);

       if( find_possible_movements_to_all_clusters ){
           clog << "SwpCndNode::create_MoveTo_SwapCandidatesForCluster, NOT IMPLEMENTED YET" << endl;
           exit(1);
       }else{

           createAllCandidatesForSwpCnd( swpval1, nw, deg_in_cl, cluster_neighbors, res,
                                         [dd = d]( int swpval, int trg_cl ){
                                             return SwpCndNode( swpval, dd, trg_cl );
                                         } );

           for( int c : cluster_neighbors) edges_to_cluster[c] = 0; // clearing used array

           if(debug){
               DEBUG(res);
               ENDL(2);
           }
       }
    }

    // if there is a swap candidate that contains whole cluster and is moved with swap value 0 to an empty cluster, then
    // we remove it from the set of swap candidates
    for( int i=(int)res.size()-1; i>=0; i-- ) {
        int v = res[i].node;
        if( state->clusters[ state->inCl[v] ].size() == 1 && res[i].move_node_to == state->getIdOfEmptyCluster() ){
            swap( res[i], res.back() );
            res.pop_back();
        }
    }

    return res;
}

ostream& operator<<(ostream& str, SwpCndNode& cnd){
    str << "(nodes --> move_to: " << cnd.node << "  -->  " << cnd.move_node_to
        << ", swp_val: " << cnd.swpVal() << ")";
    return str;
}

SwpCndNodeCreator::SwpCndNodeCreator(State &s) : SwapCandidateCreatorAdapter(s) {}

SwpCndNode::SwpCndNode(int swpval, int v, int trg_cl) {
    swap_value = swpval;
    node = v;
    move_node_to = trg_cl;
}

