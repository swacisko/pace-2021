//
// Created by sylwester on 3/9/21.
//

#include "CONTESTS/PACE21/heur/Global.h"
#include <chrono>

namespace Global{
    using namespace chrono;

    volatile sig_atomic_t tle = 0;

    int max_runtime_in_seconds = 580;

    const bool CONTEST_MODE = false;

//    const bool disable_all_logs = CONTEST_MODE; // by default should be equal to CONTEST_MODE
    const bool disable_all_logs = true; // by default should be equal to CONTEST_MODE

    bool checkTle() {
        return tle;
//        return secondsFromStart() > max_runtime_in_seconds;
    }

    void terminate(int signum) {
//        cout << "#TEST: Time from start: " << secondsFromStart() << endl;
        tle = 1;
    }

    using namespace chrono;
    high_resolution_clock::time_point start_time;

    void startAlg(){
        start_time = chrono::high_resolution_clock::now();
    }

    int secondsFromStart(){
        chrono::high_resolution_clock::time_point now = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(now - start_time);
        return time_span.count();
    }

    void addSigtermCheck(){
        struct sigaction action;
        memset(&action, 0, sizeof(struct sigaction));
        action.sa_handler = Global::terminate;
        sigaction(SIGTERM, &action, NULL);
    }

    void increaseStack(){
        const rlim_t kStackSize = 3 * 4L * 256L * 1024L * 1024L;   // min stack size = 64 Mb
        struct rlimit rl;
        int result;

        result = getrlimit(RLIMIT_STACK, &rl);
        if (result == 0)
        {
            if (rl.rlim_cur < kStackSize)
            {
                rl.rlim_cur = kStackSize;
                result = setrlimit(RLIMIT_STACK, &rl);
                if (result != 0)
                {
                    fprintf(stderr, "setrlimit returned result = %d\n", result);
                }
            }
        }
    }


}