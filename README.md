# Bayesian-Suffix-Trees

This repository contains a C++ implementation for the 3 algorithms: CTW, BCT, k-BCT. CTW returns the mean marginal likelihood of our BCT framework as the weighted probability at the root of Tmax. BCT returns the MAP tree model, and k-BCT the top-k a-posteriori most likely trees.


The file "paper.cpp" contains the main function and is ready to be built and executed. The only inputs that are needed are the ones defined as global variables at the top of the file: alphabet size m, maximum depth D, prior hyperparameter beta, and k_max for the top-k trees. Then, the only thing needed is to load the data xn as a short vector (this is done at the beginning of the main function). These are the only inputs that the user needs to add in order to run CTW and k-BCT. 



The header file "node.h" defines structures for representing nodes and trees in C++. The header file "Functions.h" contains all the required functions for the implementation of the main algorithms. Changes to these header files are not at all needed in order to run the algorithms.
