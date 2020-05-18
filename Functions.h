//#include "stdafx.h"   // this is of Visual Studio
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

//Important inputs that are set in main

extern const int m; // alphabet size
extern const int D; // max depth
extern long double beta; // prior hyper-parameter
extern const long double alpha;
extern const short k_max; //top-k trees for k-bct algorithm

#include "node.h" //definitions for nodes and trees 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// FUNCTION PROTOTYPES

//1. Functions for initialisation/preprocessing
vector<short> read(string s);                   // reads file and stores xn into vector  
vector<short> read2(string s);                  // reads file and stores xn into vector (one symbol in each line)
void init_tree(tree& T);                        // some operations to initialise tree
void preproc(vector <node *> init);             // does the preprocessing calcultions needed for k-bct



//2. Functions used for building Tmax
void update(tree& T, short s, vector <short> ct);    // updates tree for each new symbol
void occur(node * N, short s);                       // node occurrence updates
void insert(tree& T, vector <short> ct2, short ch);  // inserts node with context ct2 to T, links with existing child




//3. Important-main algorithms
long double ctw(tree& T);                      // CTW algorithm
long double bct(tree& T);                      // finds MAP tree
void ctw_bct(tree& T);                         // runs CTW and BCT together


void kbct_forw(tree &T, vector <node *> init);                      // forward pass of kbct, takes improper and gives improper tree, calculates vectors lm, c for each node
void kbct_back(vector <node *> init, tree T, vector <tree> &trees); // backward pass of kbct
void kbct(tree& T, vector <tree> &trees, vector <node *> init, vector<double> & odds); //runs ctw and kbct together




//4.The remaining secondary functions, which are auxiliary things (but are needed to perform the rest in C++)

int show_leaves(tree T);                       // prints leaves of tree and returns number of leaves
void label(tree& T);                           // sets contexts to the nodes of a tree
tree copytree(tree T);                         // copies only non-deleted nodes to store proper pruned tree


void comb_initial3(int d, vector <node *> init);                        //used in preprocessing calculations
void comb(int d, int k, tree &T, vector <node *> init);                 // combinations for k-bct forward pass
vector<vector<double> > cartesian_prod(const vector<vector<double>>& v);  //combinations for doubles
vector<vector<short> > cartesian_prod_int(const vector<vector<short>>& v);//combinations for shorts


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//FUNCTION DEFINITIONS


void update(tree& T, short s, vector <short> ct) {

	node * temp = T[0][0]; // start from root

	occur(temp, s);


	for (int j = 0; j < D; j++) {

		if (temp->child[ct[j]] != NULL) { // if child node exists in tree

			temp = temp->child[ct[j]];    // move to child
			occur(temp, s);               // child occurred
		}

		else {                            // create children up to depth D

			vector <short> ct2 = ct;        // context of node to be created
			short ch = 0;                   // shows which child exists

			for (int k = 0; k < D - j; k++) {



				//insert node with context ct2 to tree
				insert(T, ct2, ch);

				occur(T[ct2.size()].back(), s); // inserted nodes occurs

				ch = ct2.back();
				ct2.pop_back();

			}

			j = D + 5;  // break out of previous loop if a child doesn't exist
			temp->child[ch] = T[ct2.size() + 1].back();


		}

	}

}

void occur(node * N, short s) {

	N->a[s]++; //update count vectors a_s
	int M = 0;
	for (int i = 0; i < m; i++) { M = M + N->a[i]; }

	N->le = 1.0* N->le + log2(1.0* N->a[s] - 0.5) - log2(0.5*m + 1.0* M - 1.0); //update log-estimated probability

}

void insert(tree& T, vector <short> ct2, short ch) {

	int d = ct2.size();            // node depth
	node * init = new node;        // initialise a node pointer
	T[d].push_back(init);         // add a node to tree at corresponding depth



	if (d == D) {

		T[d].back()->leaf = 1;   // no children if leaf

	}

	else {                     // set address of children given by ch

		T[d].back()->child[ch] = T[d + 1].back();

	}

}




int show_leaves(tree T) { // uses node class with contexts, returns number of leaves

	int n_leaves = 0;
	int maxdepth = 0;

	for (int u = 0; u < D + 1; u++) {
		for (int v = 0; v < T[u].size(); v++)
		{

			if (T[u][v]->leaf == 1) {

				cout << u << v << " node ct is ";
				for (int m = 0; m < T[u][v]->s.size(); m++)
				{
					cout << T[u][v]->s[m];
				}


				cout << endl;
				n_leaves++;

			}

			maxdepth = u;
		}
	}

	cout << "max depth is  " << maxdepth << endl;
	cout << "number of leaves is " << n_leaves << endl;
	return n_leaves;

}






void init_tree(tree& T) { //initialise tree to have D+1 rows (i.e. for depth 0 up to D) and add root node to depth0


	vector <node*> row;
	for (int i = 0; i < D + 1; i++) { T.push_back(row); }

	T[0].push_back(new node); //only add root node which always exists


	if (D == 0) {
		T[0][0]->leaf = 1;
	}

}


long double ctw(tree& T) {                   // algorithm takes improper tree and finds mean marg lik 

	for (int d = D; d > -1; d--) {           // loop over levels

		for (int k = 0; k < T[d].size(); k++) {  // loop over nodes of each level

			if (d == D) {                   // if at max depth, node is a leaf
				T[d][k]->lw = T[d][k]->le;

			}
			else {                         // if node is not a leaf

				long double sum = 0;

				for (int ch = 0; ch < m; ch++) {

					if (T[d][k]->child[ch] != NULL) {           // if child #ch exists
						sum = sum + T[d][k]->child[ch]->lw;     // calculate sum of Le(s)

					}

				}

				//calculate weighted log-prob in two cases as explained in notes for numerical precision

				long double delta = T[d][k]->le - sum + log2(beta) - log2(1.0 - beta);
				if (delta < 30) {

					T[d][k]->lw = log2(1.0 - beta) + sum + log2(1.0 + pow(2.0, delta));

				}
				else {
					T[d][k]->lw = log2(beta) + T[d][k]->le + log2(exp(1))*(pow(2.0, -delta) - pow(2.0, -2.0*delta - 1));

				}

			}
		}
	}

	cout << "Mean marginal likelihood is " << pow(2.0, T[0][0]->lw) << " and has log " << T[0][0]->lw << endl;
	return T[0][0]->lw;              //output value of weighted prob Pw at root
}



long double bct(tree& T) {                   // algorithm takes improper tree finds maximal probability at root and makes tree proper, then finds bct
											 // lw here isused to store the maximal probablity, not the weighted one

											 // First forward pass (leaves to root) to calculate maximal probabilities Pm at every node

	if (D == 0) { // if iid data

		return  T[0][0]->le;              //output value of max prob at root
	}

	for (int d = D; d > -1; d--) {           // loop over levels



		for (int k = 0; k < T[d].size(); k++) {  // loop over initially existing nodes of each level

			if (d == D) {                   // if at max depth, node is a leaf (if and only if for improper tree)
				T[d][k]->lw = T[d][k]->le;

			}

			else {                         // if node is not a leaf

				long double sum = 0;

				for (short ch = 0; ch < m; ch++) {

					if (T[d][k]->child[ch] == NULL) {           // if child #ch does not exist, it is equivalent with

						if (d < D - 1) {
							sum = sum + log2(beta);
						}

					}

					else {                                        // if child ch exists

						sum = sum + T[d][k]->child[ch]->lw;       // sum of log-probs at children

					}

				}

				// calculate maximal log-prob as explained in notes



				if (log2(1.0 - 1.0*beta) + sum > log2(beta) + T[d][k]->le) { // maximum achieved by children term

					T[d][k]->lw = log2(1.0 - 1.0*beta) + sum;                // set max prob of node


				}

				else {                                                        // maximum achived by curent node

					T[d][k]->lw = log2(beta) + T[d][k]->le;                  // set max prob of node and mark to be pruned
					T[d][k]->leaf = 1;

					for (short ch = 0; ch < m; ch++) {   // for child # ch of each node


														 //cout << "prune child " << ch << " of node " << d << k << endl;

						if (T[d][k]->child[ch] != NULL) {                       // if child exists

							T[d][k]->child[ch]->a[0] = -1;                     // mark child to be destructed
							T[d][k]->child[ch] = 0;                            // destruct connection with child

						}

					}
				}


			}
		}
	}

	cout << endl << " end of bct forward pass" << endl;
	// Then backward pass (root to leaves), to prune tree and destroy the required nodes
	// Use a[0] =-1 to mark nodes that need to be pruned

	int length[D + 1] = { 0 };
	for (int d = 0; d < D + 1; d++) {
		length[d] = T[d].size();
	}



	for (int d = 0; d < D + 1; d++) { // root to leaves now

		int check = 0;

		for (int k = 0; k < length[d]; k++) {

			if ((T[d][k]->a[0] == -1)) {  // node was marked to be deleted

				for (short ch = 0; ch < m; ch++) {   // for child # ch of each node

					if (T[d][k]->child[ch] != NULL) {                       // if child exists

						T[d][k]->child[ch]->a[0] = -1;                     // mark child to be destructed

					}
				}

				//delete T[d][k];                                   // destruct node
				//T[d].erase(T[d].begin() + k);                    // destruct the pointer of the node from the tree
				//k--;
				//length[d]--;
				//T[d].shrink_to_fit();                          // releases memory but takes much more time, so it is commented out
			}

			else {

				if (T[d][k]->leaf == 0) {               // if child not a leaf and not deleted

					for (short ch = 0; ch < m; ch++) {

						if (T[d][k]->child[ch] == NULL) {


							node * init = new node;                // insert child to make tree proper
							T[d][k]->child[ch] = init;             // connect child with parent node
							T[d + 1].push_back(init);              // store at appropriate tree depth
							init->leaf = 1;                        // denote it leaf

							if (d < D - 1) {
								init->lw = log2(beta);            // set maximal prob for leaf at depth < D, if leaf is at depth D then logP=0;
							}


						}
					}
				}
			}

			if (T[d][k]->a[0] == -1) { //node doesn't exist in reality
				check++;
			}
		}
		if (check == length[d]) { //then no need to look at higher depths
			cout << " Pm,root is " << pow(2.0, T[0][0]->lw) << " and has log " << T[0][0]->lw << endl;
			return T[0][0]->lw;              // output value of weighted prob Pw at root
		}
	}

	cout << " Pm,root is " << pow(2.0, T[0][0]->lw) << " and has log " << T[0][0]->lw << endl;
	return T[0][0]->lw;              // output value of weighted prob Pw at root
}



void ctw_bct(tree& T) {



	long double pwl = ctw(T);

	cout << "end of CTW" << endl;

	long double pml = bct(T);

	cout << "end of bct" << endl;

	T = copytree(T);

	label(T); //add contexts to nodes



	int n_leaves = show_leaves(T);  // show leaves only


	long double prior = log2(pow(alpha, (n_leaves - 1.0))*pow(beta, (n_leaves - T[D].size())));// log-prior


	cout << "prior is " << pow(2, prior) << " and has log " << prior << endl;

	long double posterior = pml - pwl; // log-posterior

	cout << "Posterior is " << pow(2, posterior) << " and has log " << posterior << endl;




}

vector<short> read(string s) { //could use char

	vector <short> xn = { 0 };
	xn.pop_back();

	ifstream fin;
	fin.open(s);

	char c = fin.get() - '0';

	while (fin.good()) {

		xn.push_back(c);
		c = fin.get() - '0';

	}

	fin.close();


	cout << "Size of Xn is " << xn.size() << endl;


	return xn;
}

vector<short> read2(string s) { //could use char, here for m > 9, in text file: each symbol in new line, after final symbol need one extra newline

	vector <short> xn = { 0 };
	xn.pop_back();

	ifstream fin;
	fin.open(s);

	char c[3] = { 0 };

	fin.getline(c, 3);


	while (fin.good()) {

		xn.push_back(atoi(c));
		fin.getline(c, 3);

	}

	fin.close();


	cout << "Size of Xn is " << xn.size() << endl;


	return xn;
}


void label(tree& T) { // takes as input proper tree with no contexts and writes their context in node.s

	for (int d = 0; d < D + 1; d++) {

		for (int k = 0; k < T[d].size(); k++) {

			if (T[d][k]->leaf == 0) {

				for (short ch = 0; ch < m; ch++) {

					T[d][k]->child[ch]->s = T[d][k]->s;
					T[d][k]->child[ch]->s.push_back(ch);
				}


			}
		}

	}

}




tree copytree(tree T) {   // takes tree from BCT which doesn't delete nodes (only marks them so virtually they dont exist)
						  // copies only non-marked nodes to another  tree

						  // initialise tree of depth D
	tree T2;
	init_tree(T2);
	T2[0].pop_back();

	for (int d = 0; d < D + 1; d++) {

		int check = 0;

		for (int k = 0; k < T[d].size(); k++) {

			if (T[d][k]->a[0] > -1) {

				T2[d].push_back(T[d][k]);
			}

			else { check++; }
		}

		if (check == T[d].size()) {
			return T2;
		}
	}

	return T2;

}



void preproc(vector <node *> init) {   // preprocessing needed for k-bct, gives init[0] for nodes at d=1,...,init[D-1] for d=D
									   // no need to include root here as it is always in the data (obviously)


	init[D - 1]->c.push_back(zeros);       // for d=D , c=0 and lm[0]=0 from construction 

	for (short d = D - 2; d > -1; d--) {

		init[d]->lm[0] = log2(beta);      // for smaller depth first add c=0 with p=logbeta
		init[d]->c.push_back(zeros);


		comb_initial3(d, init);           // then find all combinations and keep the top k of them

	}


}


vector<vector<double> > cartesian_prod(const vector<vector<double>>& v) { // cartesian product of vectors
	vector<vector<double>> s = { {} };                                  // returns a matrix with all combinations
	for (auto& u : v) {
		vector<vector<double>> r;
		for (auto& x : s) {
			for (auto y : u) {
				r.push_back(x);
				r.back().push_back(y);
			}
		}
		s.swap(r);
	}
	return s;
}

vector<vector<short> > cartesian_prod_int(const vector<vector<short>>& v) {   // finds cartesian product of vector of short ints
	vector<vector<short>> s = { {} };                                       // stores in matrix
	for (auto& u : v) {
		vector<vector<short>> r;
		for (auto& x : s) {
			for (auto y : u) {
				r.push_back(x);
				r.back().push_back(y);
			}
		}
		s.swap(r);
	}
	return s;
}

void comb(int d, int k, tree &T, vector <node *> init) {                   // finds combinations of children of node T[d][k]
																		   // computes the products for lm and keeps top k 
	vector<vector<double>> v;  // stores combinations of probabilities
	vector<vector<short>> c;  // stores respective indices


	for (short ch = 0; ch < m; ch++) {


		node * pointer;

		if (T[d][k]->child[ch] == NULL) {

			pointer = init[d];                         // if child does not exists (not occured) use preprocessing node

		}

		else {

			pointer = T[d][k]->child[ch];                 // else use children of node

		}

		v.push_back(pointer->lm);

		vector <short> temp;

		for (short j = 1; j < pointer->lm.size() + 1; j++) { // use cartesian products of {1,2,3}... to find position vectors

			temp.push_back(j);
		}

		c.push_back(temp);
	}

	vector<vector<double> > v2 = cartesian_prod(v);          // cartesian products : all combinations
	c = cartesian_prod_int(c);

	for (int i = 0; i < v2.size(); i++) {                  //loop combiantions, keep what is necessary

		double sum = log2(1.0 - beta);

		for (short j = 0; j < v2[i].size(); j++) {
			sum = sum + v2[i][j];                          // sum of log-lm= prod-lm
		}

		vector <vector <short>> tempc;                    // use temporal variables to keep sort appropriately

		vector <double> temp;

		if (sum > T[d][k]->lm[T[d][k]->lm.size() - 1]) {  // if sum> min of ordered list lm

			int j = T[d][k]->lm.size() - 1;
			bool test = 0;

			while (j > 0) {

				if ((sum > T[d][k]->lm[j]) && (sum <= T[d][k]->lm[j - 1])) {  // if > of prev and < next, inlude there
																			  // stick element below the first one which is equal to sum
					temp.push_back(sum);
					tempc.push_back(c[i]);

					short limit = 0;

					if (T[d][k]->lm.size() < k_max) {
						limit = 1;                       // if the size of lm<k then add without replacement
					}

					for (int q = 0; q < T[d][k]->lm.size() - j - 1 + limit; q++) { //set temp
						temp.push_back(T[d][k]->lm[j + q]);
						tempc.push_back(T[d][k]->c[j + q]);
					}

					if (T[d][k]->lm.size() < k_max) {   // if size<k then add to list (here initialise space)
						T[d][k]->lm.push_back(0);
						T[d][k]->c.push_back({ 0 });
					}

					for (int q = 0; q < temp.size(); q++) {

						T[d][k]->lm[j + q] = temp[q];    // use temp for assignments
						T[d][k]->c[j + q] = tempc[q];
					}

					j = 1;            // used to break out of loop when added this combination
					test = 1;         // test used to add element at the top of the list
				}

				j--;

			}

			if (test == 0) {         // if sum> all lm then add it at the top of the list
									 // code identical with before at index 0

				temp.push_back(sum);
				tempc.push_back(c[i]);

				short limit = 0;

				if (T[d][k]->lm.size() < k_max) {
					limit = 1;
				}

				for (int q = 0; q < T[d][k]->lm.size() - 1 + limit; q++) {
					temp.push_back(T[d][k]->lm[q]);
					tempc.push_back(T[d][k]->c[q]);
				}

				if (T[d][k]->lm.size() < k_max) {
					T[d][k]->lm.push_back(0);
					T[d][k]->c.push_back({ 0 });
				}

				for (int q = 0; q < temp.size(); q++) {

					T[d][k]->lm[q] = temp[q];
					T[d][k]->c[q] = tempc[q];
				}

			}
		}

		else { // if element <= last element of list, only include it at the bottom if size<k

			if (T[d][k]->lm.size() < k_max) {
				T[d][k]->lm.push_back(sum);
				T[d][k]->c.push_back(c[i]);
			}
		}





	}


}



void kbct_forw(tree &T, vector <node *> init) { // forward pass of kbct algorithm

	for (int d = D; d > -1; d--) {               // loop for leaves to root

		for (int k = 0; k < T[d].size(); k++) {

			if (d == D) {

				T[d][k]->lm[0] = T[d][k]->le;    // at leaves set lm[0]=le with position vector c=0
				T[d][k]->c.push_back(zeros);



			}


			else {

				T[d][k]->lm[0] = log2(beta) + T[d][k]->le; // for d<D first add the c=0 lm=le combination
				T[d][k]->c.push_back(zeros);               // if sth is equal to that it will be stuck below it
														   // intuitevely keep c=0 higher to "reward" pruning at ties

				comb(d, k, T, init);                       // find combinations, sort and keep top k of them in list

			}
		}
	}

	for (short p = 0; p < T[0][0]->lm.size(); p++) {   //output at root
		cout << "prob is " << T[0][0]->lm[p] << endl;

	}
}

void kbct_back(vector <node *> init, tree T, vector <tree> &trees) {   // backward loop of kbct
																	   // builds top k trees using the lm and c's

	for (int i = 0; i < k_max; i++) { // loop over k-top trees

		cout << "i is " << i << endl;


		if (T[0][0]->c[i][0] != 0) {  // if the c at the root is not zero, then don't prune there and add nodes at d=1

			node * temp = new node;  // create new node
			*temp = *T[0][0];        // initialise from T (so have children)
			trees[i][0][0] = temp;   // add to new tree

			for (int ch = 0; ch < m; ch++) { // always add m children to have proper tree output

				if (trees[i][0][0]->child[ch] != NULL) { // if child exists, initialise new node like the child and
					node * temp2 = new node;             // add it to the tree at depth 1
					*temp2 = *trees[i][0][0]->child[ch];
					trees[i][1].push_back(temp2);
					trees[i][0][0]->child[ch] = temp2; // connect it to the root of the new tree
				}

				else {  // if child doesn't exist, initialise new node from preprocessing step and add to tree, connect to root
					node * newnode = new node;
					*newnode = *init[0];
					trees[i][1].push_back(newnode);
					trees[i][0][0]->child[ch] = newnode;
				}
			}

			for (int d = 0; d < D - 1; d++) { // after adding the m nodes at d=1, check about pruning iteratively

				for (int k = 0; k < trees[i][d].size(); k++) {

					if (trees[i][d][k]->leaf == 0) {   // if node is not a leaf

						for (int j = 0; j < m; j++) {  // for all its children (all of them exist so next line probably not necessary)

							if (trees[i][d][k]->child[j]->leaf == 0) { // probably not necessary

								int index = 0;
								if (d == 0) { index = i; }
								short t = trees[i][d][k]->c[index][j] - 1; // index t as in notes (-1 because of c++ indexing from 0 and not from 1)


																		   // check if appropriate c is not 0, then no pruning and adding all children

								if (trees[i][d][k]->child[j]->c[t][0] != 0) {


									trees[i][d][k]->child[j]->c[0] = trees[i][d][k]->child[j]->c[t]; // after node examined and children added, take t of next step from examined node
																									 // (node will not be examined again so this is ok)

									for (int ch = 0; ch < m; ch++) {

										if (trees[i][d][k]->child[j]->child[ch] != NULL) {  // if child exists 

											node * temp3 = new node;                        // create new child
											*temp3 = *trees[i][d][k]->child[j]->child[ch];  // initialise from tree T
											trees[i][d + 2].push_back(temp3);               // add it to tree
											trees[i][d][k]->child[j]->child[ch] = temp3;    // connect it to parent
										}

										else { //if child doesn't exist initialise from preprocessing and then as before

											node * newnode2 = new node;
											*newnode2 = *init[d + 1];
											if (d == D - 2) {
												newnode2->leaf = 1; // the adding tree at depth D, so mark it a leaf
											}
											trees[i][d + 2].push_back(newnode2);
											trees[i][d][k]->child[j]->child[ch] = newnode2;
										}
									}
								}

								else { // if pruning at node 

									trees[i][d][k]->child[j]->leaf = 1; // mark it to be a leaf and make all children pointer =0 so they don't exist
									for (int ch = 0; ch < m; ch++) { trees[i][d][k]->child[j]->child[ch] = 0; }

								}
							}
						}
					}
				}
			}
		}

		else { // else if T[0][0] c[i]=0 then keep only root node

			node * newnode3 = new node;
			*newnode3 = *T[0][0];
			newnode3->leaf = 1;
			for (int ch = 0; ch < m; ch++) { newnode3->child[ch] = 0; }
			trees[i][0][0] = newnode3;
		}




	}




}

void kbct(tree& T, vector <tree> &trees, vector <node *> init, vector<double> & odds) { // call all kbct functions together

	long double pwl = ctw(T);

	cout << "end of CTW" << endl;

	kbct_forw(T, init);

	cout << " end of forward pass" << endl;

	kbct_back(init, T, trees);

	cout << " end of backward pass" << endl;


	for (int i = 0; i < k_max; i++) {

		label(trees[i]);



		int n_leaves = show_leaves(trees[i]);  // choose to show leaves only


		long double prior = log2(pow(alpha, (n_leaves - 1.0))*pow(beta, (n_leaves - trees[i][D].size())));// log-prior


		cout << "prior is " << pow(2, prior) << " and has log " << prior << endl;

		cout << " Pm,root is " << pow(2.0, T[0][0]->lm[i]) << " and has log " << T[0][0]->lm[i] << endl;

		long double pml = T[0][0]->lm[i];

		long double posterior = pml - pwl; // log-posterior

		cout << "Posterior is " << pow(2, posterior) << " and has log " << posterior << endl;

		odds[i] = pow(2, T[0][0]->lm[0] - pml);

		cout << "Posterior odds = " << pow(2, T[0][0]->lm[0] - pml) << endl;
	}


}



void comb_initial3(int d, vector <node *> init) {    // finds combinations from preprocessing stage



	vector<vector<double>> v;  // stores combinations of probabilities
	vector<vector<short>> c;  // stores respective indices


	for (short ch = 0; ch < m; ch++) {


		node * pointer;
		pointer = init[d + 1];


		v.push_back(pointer->lm);

		vector <short> temp;

		for (short j = 1; j < pointer->lm.size() + 1; j++) {

			temp.push_back(j);
		}

		c.push_back(temp);
	}

	vector<vector<double> > v2 = cartesian_prod(v);
	c = cartesian_prod_int(c);

	for (int i = 0; i < v2.size(); i++) {

		double sum = log2(1.0 - beta);

		for (short j = 0; j < v2[i].size(); j++) {
			sum = sum + v2[i][j];
		}

		vector <vector <short>> tempc;

		vector <double> temp;

		if (sum > init[d]->lm[init[d]->lm.size() - 1]) {

			int j = init[d]->lm.size() - 1;
			bool test = 0;

			while (j > 0) {

				if ((sum > init[d]->lm[j]) && (sum <= init[d]->lm[j - 1])) {

					temp.push_back(sum);
					tempc.push_back(c[i]);

					short limit = 0;

					if (init[d]->lm.size() < k_max) {
						limit = 1;
					}

					for (int q = 0; q < init[d]->lm.size() - j - 1 + limit; q++) {
						temp.push_back(init[d]->lm[j + q]);
						tempc.push_back(init[d]->c[j + q]);
					}

					if (init[d]->lm.size() < k_max) {
						init[d]->lm.push_back(0);
						init[d]->c.push_back({ 0 });
					}

					for (int q = 0; q < temp.size(); q++) {

						init[d]->lm[j + q] = temp[q];
						init[d]->c[j + q] = tempc[q];
					}

					j = 1;
					test = 1;
				}

				j--;

			}

			if (test == 0) {

				temp.push_back(sum);
				tempc.push_back(c[i]);

				short limit = 0;

				if (init[d]->lm.size() < k_max) {
					limit = 1;
				}

				for (int q = 0; q < init[d]->lm.size() - 1 + limit; q++) {
					temp.push_back(init[d]->lm[q]);
					tempc.push_back(init[d]->c[q]);
				}

				if (init[d]->lm.size() < k_max) {
					init[d]->lm.push_back(0);
					init[d]->c.push_back({ 0 });
				}

				for (int q = 0; q < temp.size(); q++) {

					init[d]->lm[q] = temp[q];
					init[d]->c[q] = tempc[q];
				}

			}
		}

		else {

			if (init[d]->lm.size() < k_max) {
				init[d]->lm.push_back(sum);
				init[d]->c.push_back(c[i]);
			}
		}





	}


}

