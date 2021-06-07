//
// Created by sylwester on 3/13/21.
//

#include <CONTESTS/PACE21/kernelization/CriticalClique.h>

CriticalClique::CriticalClique(VVI &V, VI nodes, int id) {
    this->V = &V;
//    this->nodes = nodes;
    swap(this->nodes, nodes); // faster? probably compiler does that anyway
    this->id = id;
}

ostream& operator<<( ostream& str, CriticalClique& cc ){
    str << "{id: " << cc.id << ", nodes: " << cc.nodes << "}";
    return str;
}

void CriticalClique::test() {

}

