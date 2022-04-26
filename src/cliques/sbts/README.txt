Yan Jin and Jin-Kao Hao. General swap-based multiple neighborhood tabu search for finding maximum independent set. Engineering Applications of Artificial Intelligence 37: 20-33, 2015.

Important note: Please make sure that the above paper is cited whenever our code is used in your research.

Technical notes:

1. To run SBTS, you can use the following command:

     SBTS.exe   parameter1   parameter2   parameter3   
  or   
     SBTS.exe   parameter1   parameter2   parameter3   paramter4

     where parameter1, parameter2, parameter3 and parameter4 are four parameters that we will describe in detail in the following.

2. parameter1 denotes the name of the iuput file (or the name of the instance). Notice that the instance should be for the maximum independent set problem. 
   The format of the input file should be as follows (for instance):

   p edge 200 5066
   e    1   2
   e    1   3
   e    1   4
   ...

2. parameter2 denotes the name of the output file.

3. parameter3 denotes the objective value that the SBTS algorithm tries to reach.

4. parameter4 is the tabu tenure of the SBTS.


Contact Us

Please e-mail (jin@info.univ-angers.fr) under any of the following conditions:

    You're having trouble running SBTS
    You've found a bug in SBTS
    You have any other questions, comments or suggestions
