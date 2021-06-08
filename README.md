# pace-2021

**CluES** - heuristic solver for cluster editing problem, written as an entry to the PACE 2021 challenge.

A variety of heuristics are used to find a small cluster editing set of given graph.<br>
The main algorithm works in iterations, in each iteration a new cluster editing set is found (iterations are independent from each other, the longer the solver runs, the more iterations are done and the greater chance of finding a good set).<br>

The **main** branch contains code of CluES submitted to **heuristic track**.<br>
The version of CluES submitted to exact track can be found on **'exact_track'** branch. Please note: there is no guarantee of finding an optimal solution.<!-- - there are some adjustments of parameters done to the heuristic version of CluES to increase probability of finding a good result (hopefully an optimal one on instances from exact track))--><br>
The version of CluES submitted to kernelization track can be found on **'kernelization_track'** branch. Please note: there is no guarantee, that the returned kernel meets the optimality requirements.

<br>


**Requirements**:

CMake VERSION 3.10.2 or higher<br>
c++ 17 or higher

<br>

**Installation**:

Use cmake to obtain a binary file, e.g. in linux in the main directory (pace-2021) you can use the following commands:

mkdir build<br>
cd build<br>
cmake ..<br>
make

After this, the executable file named "CluES" should be in the "build" directory

<br>

**Usage:**

Given a graph in a file example_input.gr, you can run CluES in the following way
 
./CluES < example_input.gr > example_output.out 2>example_logs.err

CluES (version from heuristic track, on main branch) will run until it receives SIGTERM signal (e.g. using "kill -s SIGTERM $pid" command in linux).<br>
The exact and kernelization versions of CluES will terminate for itself and do not need SIGTERM to be sent to the process (but it will terminate earlier if the signal is sent before the compuation is finished).
<br>

**Information on running CluES on kernelization track**:<br><br>
To run the first phase (create an output of kernelization), please run CluES with the following command:<br>
./CluES -\-phase=kernelize < example_input.gr > kernel_output.out<br><br>
In order to lift the solution, please run CluES with command:<br>
./CluES -\-phase=lift-solution < lift_solution_input.in > final_output.out