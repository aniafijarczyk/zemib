
To test if scripts and programs are working make 3 qsub bash scripts, each executing 1 job (for 3 existing folders)
and test with l-test



### MODEL PREPARATION

* prepare one folder per each model with all necessary files,
each model should have name like this: M* (e.g. M1, M2, M3, M4, ...)
* in example.input file put as many simulations as you need
```
e.g. 4 models each runs 1mln simulations:
use plgrig-long (1000 jobs and 1024 core available)
run each model in 250 jobs, each job runs 4000 simulations
4 x 250 * 4000 = 4 000 000 simulations
```
* ABCtoolbox source should be compiled in zeus (this version works: ABCtoolbox_16_05, compile with: g++ -O3 -o ABCtoolbox *.cpp)
* mstatspop sourcs should be taken from mstatspop source folder
* 1 job uses 1 core



### MULTIPLY FOLDERS WITH EACH MODEL

* use multiply_jobs.py to multiply folders with jobs
* folder with model (M1) is copied to folder aniaf_M1_<int>
* (remove M1, M2 .. folders)
* multiply_jobs.py <int>



### GENERATE QSUB SCRIPTS

* generate a script (qsub_<model>.sh) that has all instructions how to execute the program somewhere
* qsub_generator.py
* use -q switch to choose the queue type (e.g. 1 is plgrid-long), the script reads 'aniaf*' folders


### GENERATE A SCRIPT WITH COMMANDLINES TO RUN ALL JOBS IN PL-GRID

* qsub_cmd.py
* this script reads all qsub_<model>.sh scripts and puts them in one general script called: kolejkuj_<folder>.sh
* the command is: 'qsub qsub_<model>.sh'




