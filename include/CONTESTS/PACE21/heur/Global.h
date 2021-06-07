//
// Created by sylwester on 3/9/21.
//

#ifndef ALGORITHMSPROJECT_GLOBAL_H
#define ALGORITHMSPROJECT_GLOBAL_H

#include <csignal>
#include <mutex>
#include <sys/resource.h>
#include "Makros.h"

namespace Global{

    extern volatile sig_atomic_t tle;

    extern int max_runtime_in_seconds;

    /**
     * If true, then no logs will be written to [clog]
     */
    extern const bool disable_all_logs;

    /**
     *
     * @return true if TLE, false otherwise
     */
    extern bool checkTle();

    extern void startAlg();

    /**
     * @return number of seconds from the start of the algorithm
     */
    extern int secondsFromStart();

    /**
     *
     * @param signum
     */
    extern void terminate(int signum);


    /**
     * Adds handling SIGTERM signals
     */
    void addSigtermCheck();

    /**
     * Sets stack size
     */
    void increaseStack();

    /**
     * Maximum execution time in seconds
     */
//    extern int max_execution_time = 600;

    /**
     * Percentage of maximal time that can be spent on creating known solutions
     */
//    extern double max_perc_time_for_known_solutions = 0.5;

    /**
     * If true, then
     */
    extern const bool CONTEST_MODE;
}

#endif //ALGORITHMSPROJECT_GLOBAL_H
