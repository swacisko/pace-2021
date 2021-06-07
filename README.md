# pace-2020

**CluES** - heuristic solver for cluster editing problem, written as an entry to the PACE 2021 challenge.

A variety of heuristics are used to find a small cluster editing set of given graph.<br>
The main algorithm works in iterations, in each iteration a new cluster editing set is found (iterations are independent from each other, the longer the solver runs, the more iterations are done and the greater chance of finding a good set).<br>


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


CluES will run until it receives SIGTERM signal (e.g. using "kill -s SIGTERM $pid" command in linux).<br>
