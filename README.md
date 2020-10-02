# 5G_user_scheduling_optimization
This project aims at designing optimal packet schedulers by providing an online solution to the NP-hard integer linear programming problem of network scheduling.

# How to test the code
We provide the command lines to execute the Java scripts in the code directory.
## Preprocessing step
arguments: String [] path, boolean display (= false)

output: Returns the results of preprocessing the file with address path. The boolean display is used to request the display of Xnkm tables after preprocessing.

example: java Preprocessing /home/user/Documents/testfiles/test1.txt true


## Dynamic Programming solution
arguments: String [] path, boolean display (= false)

output: Returns the results after applying dynamic programming algorithms to the path file. The boolean display is used to request the display of Xnkm tables after processing.

example: java DynamicPrograming /home/user/Documents/testfiles/test1.txt true


## Greedy solution
input: String [] path, boolean display (= false)

output: Returns the results after applying the greedy algorithm on the path file. The boolean display is used to request the display of Xnkm tables after processing.

example: java Greedy /home/user/Documents/testfiles/test1.txt true

## First online solution
arguments: int iterations (= 10,000)

output: Returns the results obtained after applying the online algorithm on an 'iteration' number of operations. These results are obtained compared to the greedy algorithm applied to the same problem once all the users have been generated, after the fact.

example: java Online 1000000

## Second online solution
arguments: int iterations (= 10,000)

output: Returns the results obtained after applying the online2 algorithm on a number of 'iteration' of operations. These results are obtained compared to the greedy algorithm applied to the same problem once all the users have been generated, after the fact.

example: java Online2 1000000

## Third online solution
arguments: int iterations (= 10,000)

output: Returns the results obtained after applying the online3 algorithm on a number of 'iteration' operations. These results are obtained compared to the greedy algorithm applied to the same problem once all the users have been generated, after the fact.

example: java Online3 1000000
