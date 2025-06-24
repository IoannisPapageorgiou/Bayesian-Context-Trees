# Bayesian-Context-Trees

UPDATE: We have now implemented a more user-friendly R package containins this code. Please have a look at the R package BCT, at: https://cran.r-project.org/package=BCT

This repository contains a C++ implementation for the Bayesian Context Tree (BCT) framework described in the paper http://arxiv.org/abs/2007.14900. The main algorithms implemented here are CTW, BCT and k-BCT. CTW returns the prior predictive likelihood (averaged over models and parameters) as the weighted probability at the root of Tmax. BCT returns the MAP tree model, and k-BCT the top-k a-posteriori most likely models. Also, a function calculating the log-loss incurred from sequential prediction is included.


The file "main_bct.cpp" contains the main function and is ready to be built and executed. The only inputs that are needed are the ones defined as global variables at the top of the file: alphabet size m, maximum depth D, prior hyperparameter beta, and k_max for the top-k trees. Then, the only thing needed is to load the data xn as a short vector (this is done at the beginning of the main function).  



The header file "node.h" defines structures for representing nodes and trees in C++. The header file "Functions.h" contains all the required functions for the implementation of the main algorithms. Changes to these header files are not at all needed in order to run the algorithms.
