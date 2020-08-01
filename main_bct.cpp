//#include "stdafx.h"     //This is for Visual Studio
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>  
#include <iomanip>    
#include <vector>
#include <random>
#include <map>
#include <time.h> 

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//IMPORTANT INPUTS                              

const int m = 2; // alphabet size
const int D = 10; // max depth
long double beta = 0.5; // prior hyper-parameter
const long double alpha = pow((1.0 - beta), (1.0 / (m - 1.0)));
const short k_max = 5; //top-k trees for k-bct algorithm

//1.These are the inputs that need to be changed by a user for running the code.
//2.The ONLY remaining thing is to load the data file in the first line of the main function.


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Functions.h"

//MAIN


int main() {

	//1. Read data xn
	string s = "example.txt";                      // After setting the initial inputs at the top of the file,
	vector <short> xn = read(s);                   // this is the ONLY remaining thing needed to run the code for CTW, k-BCT.
	cout << "File was read " << endl;           
	
	// !!! If sequential prediction is required instead of model selection, uncomment the next 2 lines !!!  

	//int train_size = 500;                       // Size of training set                    
	//log_loss(xn, train_size); return 0;         // This evaluates the log-loss and writes it to file "log_loss.txt"


	//2. Initialise tree to store Tmax 
	tree T;             
	init_tree(T); 
	vector <tree> trees(k_max, T);//initialise top-k trees



	//3. Find pre-processing nodes at each depth

	vector <node *> init; // pointers for nodes of pre-processing stage (not needed  for root node)

	if (D > 0) {

		for (short d = 0; d < D; d++) {
			init.push_back(new node); // initialise them
		}
		preproc(init);

	}




	//4. Update for each sequence symbol to build Tmax and calculate estimated probabilities

	for (int i = D; i < xn.size(); i++) {
	
		short s = xn[i];          // current symbol
		vector <short> ct(D);     // current context

		for (int j = 0; j < D; j++) {

			ct[j] = xn[i - j - 1];  // sets context
									
		}
		
		update(T, s, ct);           // Updates sequentially Tmax when symbol "s" follows context "ct" in xn
	}

	cout << "Tree was built" << endl;	


	//5. Run CTW and k-BCT

	vector<double> odds(k_max, 0.0); // stores posterior odds for top-k trees
	kbct(T, trees, init, odds);      // runs CTW followed by k-BCT

	cout << endl << "End of k-BCT" << endl;

}









