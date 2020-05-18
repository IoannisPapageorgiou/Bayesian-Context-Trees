#include <vector>
using namespace std;

extern const int m; // alphabet size , initialised in main file


typedef vector <vector <short > > matrix;

//Introduce a class for nodes

class node 
{
public:

	vector <short> s = {};               // Node context

	int a[m] = { 0 };                    // Occurrences for each j in alphabet

	double le = 0;                       // Logarithm of estimated probability Pe(as) 
	double lw = 0;                       // Logarithm of weighted probability probability
	vector <double> lm = { 0 };          // List  of log-maximal probabilities
	matrix c;                            // Position vectors for k-bct



	bool leaf = 0;                       // Indicates if node is a leaf

	node * child[m] = { 0 };             // Pointers for children nodes

};


// Introduce a structure for a tree: group of node pointers

typedef vector <vector <node *> >  tree; // description: tree[m] stores pointers of nodes at depth m



const vector <short> zeros(m, 0);
