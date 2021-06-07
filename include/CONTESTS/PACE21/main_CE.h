//
// Created by sylwester on 3/8/21.
//

#ifndef ALGORITHMSPROJECT_MAIN_CE_H
#define ALGORITHMSPROJECT_MAIN_CE_H

#include <graphs/GraphReader.h>
#include <CONTESTS/PACE21/heur/PaceUtils.h>
#include <CONTESTS/PACE21/heur/ExpansionOrder.h>
#include <CONTESTS/PACE21/heur/Cluster.h>
#include <CONTESTS/PACE21/heur/Config.h>
#include <CONTESTS/PACE21/heur/State.h>
#include <CONTESTS/PACE21/heur/SwapCandidates/SwapCandidate.h>
#include <CONTESTS/PACE21/heur/EOCreators/ComponentExpansion.h>
//#include <CONTESTS/PACE21/heur/EOCreators/FlowCutter.h>
#include <CONTESTS/PACE21/heur/Global.h>
#include <CONTESTS/PACE21/kernelization/CEKernelizer.h>
#include <CONTESTS/PACE21/kernelization/CriticalClique.h>
#include <graphs/GraphUtils.h>
#include <utils/RandomNumberGenerators.h>
#include <graphs/GraphTrimmer.h>
#include <CONTESTS/PACE21/heur/Solver.h>

#include "Makros.h"

void kernelizationCompare();

void main_CE(int argc, char **argv);

#endif //ALGORITHMSPROJECT_MAIN_CE_H
