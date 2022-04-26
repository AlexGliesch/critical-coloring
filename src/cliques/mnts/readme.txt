1. To run MN/TS, you can use the following command: 
    mnts.exe   parameter1   parameter2   parameter3   paramter4 
    where parameter1, parameter2, parameter3 and parameter4 are four parameters that we will describe in detail in the following.
2. parameter1 denotes the name of the iuput file (or the name of the instance). The format of the input file should be as follows (for instance): 
   p 1000 449449
   e 2 1
   e 3 1
   e 3 2
   e 4 1
   ...
3. parameter2 denotes the objective value (or Wbest in the original paper) that the MN/TS algorithm tries to reach. 
4. parameter3 is for allocating weight, for vertex i, of i mod parameter3 + 1. In the original paper, parameter3 is set to 200, 10 and 1 (unweighted case) respectively.
5. parameter4 is a very important parameter which denotes the search depth of the tabu search. In the original paper, parameter4 is set as follows. parameter4 = 4000 for the instances of the weighted case (MWCP). For the unweighted case (MCP), parameter4 = 10000 except for the brock and san families (DIMACS) for which parameter4 is equal to 100.
6. Please e-mail to us (wu@info.univ-angers.fr) if You're having trouble running MN/TS or You've found a bug in MN/TS.
7. For academic puposes, you may use the software as you wish. If you wish to use this software for commercial applications, please obtain the prior permission of Qinghua Wu (wu@info.univ-angers.fr). 
 

