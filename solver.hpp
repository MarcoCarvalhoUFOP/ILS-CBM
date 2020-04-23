#include <random>
#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <list>
#include <queue>
#include <bitset>
using namespace std;

#include <chrono>
using namespace std::chrono;

double TIME_LIMIT;
float ALPHA = 0.2;
float BETA = 0.99;
float GAMMA = 0.99;
int MAX_ITER = 150;
int RUNS  = 20;
#define MAX 1500

class Element{

public:
	Element(){
			blocked = false;
			removed = false;
			sz = 0;
	}

	void setBlocked()
	{
		blocked = true;
	}

	bool isBlocked()
	{
		return blocked;
	}

	void setRemoved()
	{
		removed = true;
	}

	bool isRemoved()
	{
		return removed;
	}

	void incrementSize()
	{
		sz++;
	}

	int size()
	{
		return sz;
	}

private:
	bool blocked; //indicates whether the element has been already sequenced
	bool removed; //indicates whether the element has been removed in the preprocessing
	int sz;    //stores the size of the element, in number of pieces
};

//Node class for use in graph search
class Node{
	friend bool operator<(const Node &, const Node &);
	friend bool operator>(const Node &, const Node &);

public:
	Node()
	{
		degree = 0;
		visited = false;
		done = false;
		id = 0;
	}

	void setDegree(int i)
	{
		degree = i;
	}

	void setVisited()
	{
		visited = true;
	}

	void setDone()
	{
		done = true;
	}

	void setId(int i)
	{
		id = i;
	}

	int getDegree()
	{
		return degree;
	}

	bool wasVisited()
	{
		return visited;
	}

	bool isDone()
	{
		return done;
	}

	int getId()
	{
		return id;
	}

private:
	int id;		 //stores the index of the node
	bool visited; //indicates whether the node has been reached by the BFS
	bool done;
	int degree;	 //stores the degree of the node
};

//GLOBAL variables
std::default_random_engine generator;
vector < bitset <MAX> > bitMatrix;
vector< vector<float> > adjacencyMatrix;
vector< int > pi;
vector< int > permutation;
vector< pair<int, int> > pairs;
vector< Element > element;
vector< vector<int> > inputMatrix;
vector< list <int> > dominator;
vector<Node> node;
int m, n, nOriginal, lowerBound, convergence;
int pp, initialNode, searchRule, shaking, localSearch, mh, verbose;
unsigned seed;
Node pivot;
high_resolution_clock::time_point ILSt1, ILSt2;
duration<double> ILStimeSpan;

bool operator<(const Node &a, const Node &b)
{
	if(searchRule == 1)
		return a.degree > b.degree;
	else if (searchRule == 2)
				return a.degree < b.degree;
			else if (searchRule == 3){
				int i, j, k, l, aux;

				i = pivot.id;
				k = i;
				j = a.id;
				l = a.id;

				if(j > i){
					aux = i;
					i = j;
					j = aux;
				}

				if(l > k){
					aux = k;
					k = l;
					l = aux;
				}

				return adjacencyMatrix[i][j] > adjacencyMatrix[k][l];
			}

			return false;
}

bool operator>(const Node &a, const Node &b)
{
	return a.degree > b.degree;
}

/*
Initializes the structures, vectors and matrices used
*/
void initialization()
{
	int i;

	TIME_LIMIT = n/2;

	inputMatrix.resize(m);

	for (i = 0; i < m; i++) {
			inputMatrix[i].resize(n);
	}

	adjacencyMatrix.resize(m);

	for (i = 0; i < m; i++) {
			adjacencyMatrix[i].resize(m);
	}

	node.resize(m);

	for (i = 0; i < m; i++)
		node[i].setId(i);

	bitMatrix.resize(n+2);
	dominator.resize(n);
	element.resize(n);

	seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);

	convergence = 0;
}

/*
Reads the problem from a file specified by fileName
*/
void readProblem(string fileName)
{
	int i;
	int j;
	int k;

	char nome[256];
	FILE* fpIn = fopen(fileName.c_str(), "r");                      //input file

	fscanf(fpIn, "%d %d", &m, &n);

	nOriginal = n;

	initialization();                                       //initializes all structures, vectors and matrices

	for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++) {
					fscanf(fpIn, "%d", &inputMatrix[i][j]);
					if (inputMatrix[i][j] != 0) {
							inputMatrix[i][j] = 1;
							bitMatrix[j+1][i]=1;
							element[j].incrementSize();
					}
			}
	}

	fclose(fpIn);
}

/*
Builds the graph
*/
void buildGraph()
{
	int i, j, k;

		for (j = 0; j < n; j++) {
			for (i = 0; i < m; i++) {
					if ((!element[j].isRemoved()) && (inputMatrix[i][j] == 1)) {
							for (k = 0; k < m; k++) {
									if ((i != k) && (inputMatrix[k][j] == 1)) {  //updates the adjacency inputMatrix
											if (k < i) {
													adjacencyMatrix[k][i]+=0.5;
											}
											else
											{
													adjacencyMatrix[i][k]+=0.5;
											}
									}
							}
					}
			}
	}

	for (i = 0; i < m; i++) {
			for (j = 0; j < m; j++) {
					if (i < j && adjacencyMatrix[i][j]>0) {
							node[i].setDegree(node[i].getDegree()+1);
					}
					else if (i > j && adjacencyMatrix[j][i]>0) {
							node[i].setDegree(node[i].getDegree()+1);
					}
			}
	}
}

int deltaBitwise(int a, int b, int c)
{
	return (~bitMatrix[a]&bitMatrix[b]&~bitMatrix[c]|bitMatrix[a]&~bitMatrix[b]&bitMatrix[c]).count();
}

int deltaLuisShift (int i, int j)
{
int jLeft, jRight, iRight, iLeft;

if(verbose == 1)
	cout<<"delta pair: "<<i<<" "<<j<<endl;

if(j+1 >= n)
	jRight = nOriginal+1;
else
	jRight = permutation[j+1]+1;

if (j-1 < 0)
		jLeft = 0;
else
		jLeft = permutation[j-1]+1;

		if(i+1 >= n)
				iRight = nOriginal+1;
		else
				iRight = permutation[i+1]+1;

		if (i-1 < 0)
				iLeft = 0;
		else
				iLeft = permutation[i-1]+1;


i = permutation[i];
j = permutation[j];
if(verbose == 1)
{
	cout<<"Delta evaluation started! We're interested in swapping columns "<<i<<" and "<<j<<endl;
	cout<<"Permutation:"<<endl;
	for (int k = 0; k < n; k++)
			cout<<permutation[k]<<" ";
	cout<<"\nBitmatrix:"<<endl;
	for (int k = 0; k < nOriginal+2; k++)
		if (!element[k].isRemoved())
			cout<<bitMatrix[k]<<endl;
}

i++;
j++;

//cout<<"\nparameters -> iLeft:"<<iLeft<<" i:"<<i<<" iRight: "<<iRight<<" jLeft:"<<jLeft<<" j:"<<j<<" jRight:"<<jRight<<endl;

return -deltaBitwise(iLeft, i, iRight) + deltaBitwise(jLeft, i, j);
}

int deltaLuisExchange(int i, int j)
{
int jLeft, jRight, iRight, iLeft;

if(verbose == 1)
	cout<<"delta pair: "<<i<<" "<<j<<endl;

if(j+1 >= n)
	jRight = nOriginal+1;
else
	jRight = permutation[j+1]+1;

if (j-1 < 0)
		jLeft = 0;
else
		jLeft = permutation[j-1]+1;

		if(i+1 >= n)
				iRight = nOriginal+1;
		else
				iRight = permutation[i+1]+1;

		if (i-1 < 0)
				iLeft = 0;
		else
				iLeft = permutation[i-1]+1;


i = permutation[i];
j = permutation[j];
if(verbose == 1)
{
	cout<<"Delta evaluation started! We're interested in swapping columns "<<i<<" and "<<j<<endl;
	cout<<"Permutation:"<<endl;
	for (int k = 0; k < n; k++)
			cout<<permutation[k]<<" ";
	cout<<"\nBitmatrix:"<<endl;
	for (int k = 0; k < nOriginal+2; k++)
		if (!element[k].isRemoved())
			cout<<bitMatrix[k]<<endl;
}

i++;
j++;

//cout<<"\nparameters -> iLeft:"<<iLeft<<" i:"<<i<<" iRight: "<<iRight<<" jLeft:"<<jLeft<<" j:"<<j<<" jRight:"<<jRight<<endl;

return -deltaBitwise(iLeft, i, iRight) - deltaBitwise(jLeft, j, jRight) + deltaBitwise(jLeft, i, jRight) + deltaBitwise(iLeft, j, iRight);
}

int deltaLuisInversion(int i, int j)
{
int jLeft, jRight, iRight, iLeft;

if(verbose == 1)
	cout<<"delta pair: "<<i<<" "<<j<<endl;

if(j+1 >= n)
	jRight = nOriginal+1;
else
	jRight = permutation[j+1]+1;

if (j-1 < 0)
		jLeft = 0;
else
		jLeft = permutation[j-1]+1;

		if(i+1 >= n)
				iRight = nOriginal+1;
		else
				iRight = permutation[i+1]+1;

		if (i-1 < 0)
				iLeft = 0;
		else
				iLeft = permutation[i-1]+1;


i = permutation[i];
j = permutation[j];

if(verbose == 1)
{
	cout<<"Delta evaluation started! We're interested in swapping columns "<<i<<" and "<<j<<endl;
	cout<<"Permutation:"<<endl;
	for (int k = 0; k < n; k++)
			cout<<permutation[k]<<" ";
	cout<<"\nBitmatrix:"<<endl;
	for (int k = 0; k < nOriginal+2; k++)
		if (!element[k].isRemoved())
			cout<<bitMatrix[k]<<endl;
}

i++;
j++;

//cout<<"\nparameters -> iLeft:"<<iLeft<<" i:"<<i<<" iRight: "<<iRight<<" jLeft:"<<jLeft<<" j:"<<j<<" jRight:"<<jRight<<endl;

return -deltaBitwise(iLeft, i, iRight) - deltaBitwise(jLeft, j, jRight) + deltaBitwise(iRight, i, jRight) + deltaBitwise(iLeft, j, jLeft);
}

/*
Evaluates the current solution
*/
int evaluation()
{
	int value = 0;

	for (int i = 0; i < m; i++){
		int cont=0;
		for (int j = 0; j < permutation.size(); j++)
		{
				if(inputMatrix[i][permutation[j]] == 1)
					cont++;
	  		if (inputMatrix[i][permutation[j]] == 0){
					if(cont >0){
						value++;
						cont = 0;
				  }
		    }
    }

		if(cont > 0)
		  value ++;
	}

	return value;						//returns the maximum number of open stacks
}

/*
 Pre-processing procedure
 */
void preProcessingRepeated()
{
	int i, j, k, index, index2;
	int flag = 1;
	int counter = 0;

	if(verbose == 1)
		cout<<"Preprocessing started, removing duplicated columns"<<endl;

	for (i = 0; i < n; i++) {
			if (!element[i].isRemoved()) {
					for (j = i + 1; j < n; j++) {
							if (!element[j].isRemoved()) {
								if (element[i].size() == element[j].size()) {               //if elements i and j have the same size
										if(verbose == 1)
											cout<<"Columns "<<i+1<<" and "<<j+1<<" have the same size";
										index = j;
										index2 = i;
										for (k = 0; k < m; k++) {
												if (inputMatrix[k][i] != inputMatrix[k][j]) {             //and they differ in any position
													if(verbose == 1)
														cout<<"... but are different at row "<<k+1<<endl;

														flag = 0;                                   //then no one is dominated
														break;
												}
										}
										if (flag == 1) {                                                //if the a dominance condition is met
											if(verbose == 1)
												cout<<"... and are equal"<<endl;
											element[index].setRemoved();                                //the element is removed
											if(verbose == 1)
												cout<<"Column "<<index+1<<" removed"<<endl;
											counter++;

											dominator[index2].push_back(index);
											if(verbose == 1)
												cout<<"Column "<<index+1<<" added to the dominated list of column "<<index2+1<<endl;

											if (dominator[index].size() > 0)      {                     //if the dominated element dominates other elements
													dominator[index2].splice(dominator[index2].end(), dominator[index]); //VERIFICAR
													//transferDominated(&dominator[index2], &dominator[index]);       //the list is transferred
													if(verbose == 1)
														cout<<"Transferring the columns dominated by "<<index2+1<<endl;
											}
											if (index == i)
												break;
						   		}
								}
								flag = 1;
							}
						}
				}
		}
		nOriginal = n;                                                  //stores the original number of elements
		n -= counter;                                                           //updates the number of elements

		if(verbose == 1)
			cout<<counter<<" columns removed"<<endl;
}

bool compareColumns(int i, int j)
{
	return element[i].size()>element[j].size();
}

void preProcessingDominated()
{
	int i, j, k, index, index2;
	int flag = 1;
	int counter = 0;
	vector<int> dominanceOrdering;

	for(i = 0; i < nOriginal; i++)
		if (!element[i].isRemoved())
			dominanceOrdering.push_back(i);

	sort(dominanceOrdering.begin(), dominanceOrdering.end(), compareColumns);

	if(verbose == 1)
		cout<<"Preprocessing started, removing dominated columns"<<endl;

	for (i = 0; i < dominanceOrdering.size(); i++) {
			if (!element[dominanceOrdering[i]].isRemoved()) {
					for (j = i + 1; j < dominanceOrdering.size(); j++) {
							if (!element[dominanceOrdering[j]].isRemoved()) {
								if ((element[dominanceOrdering[i]].size() > element[dominanceOrdering[j]].size())) {      //if element i is larger than element j
									if(verbose == 1)
										cout<<"Column "<<dominanceOrdering[i]+1<<" is larger than column "<<dominanceOrdering[j]+1;
									index = dominanceOrdering[j];
									index2 = dominanceOrdering[i];
									for (k = 0; k < m; k++) {
										if ((inputMatrix[k][dominanceOrdering[i]] == 0) && (inputMatrix[k][dominanceOrdering[j]] == 1)) { //and j has a piece that i doesn't
											if(verbose == 1)
												cout<<"... both are different at row "<<k+1<<endl;
											flag = 0;                                   						//then j is not dominated by i
											break;
										}
									}
								}
								else flag = 0;
								if (flag == 1) {                                                			//if the a dominance condition is met
									if(verbose == 1)
											cout<<"... and dominates it"<<endl;

									element[index].setRemoved();                                			//the element is removed
									if(verbose == 1)
										cout<<"Column "<<index+1<<" removed"<<endl;
									counter++;

									dominator[index2].push_back(index);

									if(verbose == 1)
										cout<<"Column "<<index+1<<" added to the dominated list of column "<<index2+1<<endl;

									if (dominator[index].size() > 0)      {                     //if the dominated element dominates other elements
											dominator[index2].splice(dominator[index2].end(), dominator[index]); //VERIFICAR

											if(verbose == 1)
												cout<<"Transferring the columns dominated by "<<index2+1<<endl;
									}
									if (index == i)
										break;
							}
						  flag = 1;
						}
				}
			}
		}

		if(pp == 2)
			nOriginal = n;                                                  //stores the original number of elements
		n -= counter;                                                           //updates the number of elements

		if(verbose == 1)
			cout<<counter<<" columns removed"<<endl;
}

/*
Terminates all data structures.
*/
void termination()
{
	int i;

	for (i = 0; i < m; i++) {
			inputMatrix[i].clear();
	}

	inputMatrix.clear();

	for (i = 0; i < m; i++) {
			adjacencyMatrix[i].clear();
	}

	adjacencyMatrix.clear();

	node.clear();

	for (i = 0; i < n; i++) {
			dominator[i].clear();
	}

	for (i = 0; i < n; i++)
		bitMatrix[i].reset();

	dominator.clear();
	element.clear();
	permutation.clear();
	pi.clear();
	pairs.clear();

	m=0;
	n=0;
	nOriginal=0;
	lowerBound=0;
}

/*
Prints the solution information to the output file
*/
void printSolution(string inputFileName, int initialSolutionValue, int finalSolutionValue, double initialRunningTime, double finalRunningTime, int run)
{
		int i, j;

    string outputFileName = "Solution_run_"+to_string(run)+"_" + inputFileName;
		string permutationFileName = "Permutation_run_"+to_string(run)+"_" + inputFileName;

		ofstream fpPermutation(permutationFileName);

		fpPermutation <<inputFileName<<endl;

		for(int i = 0; i < permutation.size(); i++){
				fpPermutation << permutation[i] << " ";
		}

		fpPermutation.close();

		ofstream fpSolution(outputFileName);				//file that contains the information about the solution of a problem instance

		fpSolution << "Algorithm settings:"<< endl;
		if(verbose == 1)
			fpSolution << "Verbose mode on"<< endl;

		fpSolution<<"Alpha: "<<ALPHA<<endl;
		fpSolution<<"Beta: "<<BETA<<endl;
		fpSolution<<"Gamma: "<<GAMMA<<endl;
		fpSolution<<"Max iterations: "<<MAX_ITER<<endl;

		switch(pp){
		 	case 1:		fpSolution<<"Preprocessing repeated columns"<<endl; break;
			case 2:		fpSolution<<"Preprocessing dominated columns"<<endl;  break;
			case 3: 	fpSolution<<"Preprocessing repeated and dominated columns"<<endl;  break;
			default: 	fpSolution<<"No preprocessing"<<endl;
		}
		switch(searchRule){
			case 1:		fpSolution<<"Search guided by minimum degree"<<endl; break;
			case 2:		fpSolution<<"Search guided by maximum degree"<<endl;  break;
			case 3: 	fpSolution<<"Search guided by maximum edge weight"<<endl;  break;
			default: 	fpSolution<<"No search guidance"<<endl;
		}
		switch(initialNode){
			case 1:		fpSolution<<"BFS started by minimum degree node"<<endl; break;
			case 2:		fpSolution<<"BFS started by maximum degree node"<<endl;  break;
			case 3: 	fpSolution<<"BFS started by random node"<<endl;  break;
		 	default: 	fpSolution<<"No choice of initial node"<<endl;
	 	}
		switch(shaking){
			case 1:		fpSolution<<"Shaking: 2-opt"<<endl; break;
			case 2:		fpSolution<<"Shaking: 2-swap"<<endl; break;
			default: 	fpSolution<<"No choice of shaking"<<endl;
	  }
	  switch(localSearch){
			case 1:		fpSolution<<"Local Search: 2-opt"<<endl; break;
			case 2:		fpSolution<<"Local Search: 2-swap"<<endl; break;
			case 3:  	fpSolution<<"Local Search: 1-block grouping"<<endl; break;
			case 4:  	fpSolution<<"Local Search: 2-opt and 1-block grouping"<<endl; break;
			case 5:  	fpSolution<<"Local Search: 2-swap and 1-block grouping"<<endl; break;
			default: 	fpSolution<<"No choice of local search"<<endl;
	  }
	  switch(mh){
			 case 1:		fpSolution<<"Metaheuristic Search: ILS"<<endl; break;
			 case 2:		fpSolution<<"Metaheuristic Search: SA"<<endl; break;
			 default: 	fpSolution<<"No choice of Metaheuristic Search"<<endl;
    }

		fpSolution << endl;

		fpSolution << "Input File: "<< inputFileName << endl << endl;
		fpSolution << "Instance properties: \n" << "# Rows: " << m << "\n# Columns: "  << nOriginal << endl <<endl;
		fpSolution << "# Consecutive blocks: \n"<< "BFS: "<< initialSolutionValue<< "\nILS: "<<finalSolutionValue << endl<<endl;
		fpSolution << "Running time (s): \n" << "BFS: "<< initialRunningTime<< "\nILS: "<<finalRunningTime << endl;
		fpSolution << "Run: "<<run<< "/"<<RUNS << endl;
		fpSolution <<"Metaheuristic convergence: "<< convergence <<endl;
		fpSolution << "\nSolution: ";

		for(int i = 0; i < permutation.size(); i++){
				fpSolution << permutation[i] << " ";
		}

		fpSolution << endl << endl;
		fpSolution << "Original matrix: " << endl;

		for(int i = 0; i < m; i ++){
				for(int j = 0; j < nOriginal; j++){
						fpSolution <<(int) inputMatrix[i][j] << " ";
				}
				fpSolution << endl;
		}

		fpSolution << endl << endl;
		fpSolution << "Solution matrix: " << endl;

		for(int i = 0; i < m; i++){
				for(int j = 0; j < nOriginal; j++){
						fpSolution << (int)inputMatrix[i][permutation[j]] << " ";
				}
				fpSolution << endl;
		}

		fpSolution.close();
}

/*
Sets the lowerbound
*/
void setLowerBound()
{
	lowerBound = m;
}

/*
Gives the column sequencing correspondent to the row sequencing, original version
*/
void columnSequencing()
{
    int i, j;

		if(verbose == 1)
			cout<<"Column sequencing started"<<endl;

    for (i = m - 1; permutation.size()<n; i--) {                    				//traverses the pieces sequence in reversal order
			if(verbose == 1)
				cout<<"Checking which columns have a nonzero at row "<<pi[i]+1<<endl;

			for (j = 0; (j < nOriginal && permutation.size()<n); j++) { 			//finds the element that include the current piece
						if(verbose == 1){
							cout<<"Checking column "<<j+1<<endl;
							cout<<"n: "<<n<<" i:"<<i+1<<" permutation.size():"<<permutation.size()<<endl;
						}
            if ((!element[j].isBlocked()) && (!element[j].isRemoved())) { //if the element wasn't already sequenced and wasn't removed by the preprocessing
                if (inputMatrix[pi[i]][j] == 1) {
									if(verbose == 1)
										cout<<"--- Column "<<j+1<<endl;
                  element[j].setBlocked();
									permutation.push_back(j);
                }
            }
        }
    }
}

void addDominatedColumns()
{
	int i;
	vector<int>::iterator it;
	vector <int> aux(permutation);

	for(i=0; i<aux.size(); i++)
	{
			if (!dominator[aux[i]].empty()){
				if(verbose == 1){
					cout<<"\nSolution before:"<<endl;
					for (int j=0; j<permutation.size();j++)
						cout<<permutation[j]<<" ";
					cout<<endl;
				}

				it = find(permutation.begin(), permutation.end(), aux[i]);
				permutation.insert(it+1, dominator[aux[i]].begin(), dominator[aux[i]].end());
				dominator[aux[i]].clear();

				if(verbose == 1){
					cout<<"Columns dominated by column "<<aux[i]<<" added to the solution"<<endl;
					for (int j=0; j<permutation.size();j++)
						cout<<permutation[j]<<" ";
					cout<<endl;
				}
			}
	}
}

int findMinimumDegree()
{
    int index  = -1;
    int i;
    int smaller = m * m + 1;

    for (i = 0; i < m; i++) {
        if ((node[i].getDegree() < smaller) && (!node[i].isDone()) && (!node[i].wasVisited())) {
            smaller = node[i].getDegree();                                                 //selects the node of smaller degree not sequenced yet
            index = node[i].getId();                                                                  //stores its index
        }
    }

    return index;
}

int findMaximumDegree()
{
    int index  = -1;
    int i;
    int larger = -1;

    for (i = 0; i < m; i++) {
        if ((node[i].getDegree() > larger) && (!node[i].isDone()) && (!node[i].wasVisited())) {
            larger = node[i].getDegree();                                                 //selects the vertex of smaller degree not sequenced yet
            index = node[i].getId();                                                                  //stores its index
        }
    }

    return index;
}

int getRandomNode()
{
		std::uniform_int_distribution<int> distribution(1,m-1);

		return distribution(generator);
}

void breadthFirstSearch()
{
    int i, index;
    priority_queue <Node> preQ;
		list <Node> Q;

		if(verbose == 1)
			cout<<"BFS started"<<endl;

		if(initialNode == 1)
			Q.push_back(node[findMinimumDegree()]);
		else
			Q.push_back(node[findMaximumDegree()]);

		if(verbose == 1)
			cout<<"Initial node: "<<Q.front().getId()+1<<" degree: "<<Q.front().getDegree()<<endl;
    do {
						if(verbose == 1)
							cout<<"Nodes and degrees in the priority queue"<<endl;
						while(!preQ.empty())
						{
								Q.push_back(preQ.top());
								if(verbose == 1)
									cout<<Q.back().getId()+1<<" ->"<<	 Q.back().getDegree()<<endl;
								preQ.pop();
						}


				if (preQ.size()>0 && verbose) {
						cout<<"Priority queue tranferred to queue"<<endl;
        }

				if(Q.size() == 0)
				{
					if(searchRule == 2)
						Q.push_back(node[findMaximumDegree()]);                                        //finds the node with maximum degree and enqueue it
					else
						Q.push_back(node[findMinimumDegree()]);                                        //finds the node with maximum degree and enqueue it
					if(verbose == 1)
						cout<<"Empty queue. The new starting node is "<<Q.front().getId()+1<<" with degree "<<node[Q.front().getId()].getDegree()<<endl;
				}

				pi.push_back(Q.front().getId());                                               				//dequeue and inserts into the pieces sequence
				if(verbose == 1)
					cout<<"Node inserted in the sequence:  "<<pi.back()+1<<endl;
				index = Q.front().getId();
				if(searchRule == 3)
					pivot = node[index];
				Q.pop_front();

				if(verbose == 1)
					cout<<"Exploring the neighborhood of node "<<index+1<< " with degree "<<node[index].getDegree()<<endl;
        for (i = 0; i < index && 0 < node[index].getDegree(); i++) { 								//searches as long as the degree of the node allow
            if (adjacencyMatrix[i][index] > 0) {                                   //searches for its neighbours
								if(verbose == 1)
									cout<<"-- Node "<<i+1<<" is a neighbor of node "<<index+1<<endl;
                if ((!node[i].wasVisited()) && (!node[i].isDone())) {               //if they are not in queue or pieces sequence
										if(verbose == 1)
											cout<<"-- Node "<<i+1<<" marked as visited"<<endl;
                    node[i].setVisited();                                           //blocks the neighbour
										preQ.push(node[i]);

										if(verbose == 1)
											cout<<"Node "<<i+1<<" added to the priority queue"<<endl;
                }
            }
        }

        for (i = index + 1; i < m && 0 < node[index].getDegree(); i++) { //same as the loop above, but the indexes of the adjacency matrix are reversed (remember it is an upper diagonal matrix)
            if (adjacencyMatrix[index][i] > 0) {
								if(verbose == 1)
									cout<<"-- Node "<<i+1<<" is a neighbor of node "<<index+1<<endl;
                if ((!node[i].wasVisited()) && (!node[i].isDone())) {
										if(verbose == 1)
											cout<<"-- Node "<<i+1<<" marked as visited"<<endl;
                    node[i].setVisited();
										preQ.push(node[i]);
										if(verbose == 1)
											cout<<"Node "<<i+1<<" added to the priority queue"<<endl;
                }
            }
        }
        node[index].setDone();                                                  //removes the dequeued nodes from the problem

				if(verbose == 1){
					cout<<"Node "<<index+1<<" removed from the search"<<endl;
					cout<<"pi.size() "<<pi.size()<<endl;
					cout<<"pi contents: ";
					for(i = 0; i<pi.size(); i++)
						cout<<pi[i]<<" ";
					cout<<endl<<endl;
				}
    }
		while (pi.size() < m);                                                  //until all nodes are in the sequence

		Q.clear();

		while(!preQ.empty())
			preQ.pop();
}

bool checkPermutation()
{
		sort(permutation.begin(), permutation.end());
		for (int i= 0; i<permutation.size(); i++)
			if(permutation[i]!=i)
				return false;

		return true;
}

void generatePairs(){

	string bitmask(2, 1); 				// K leading 1's
	bitmask.resize(n, 0); 				// N-K trailing 0's
	pair <int,int> par;
	vector<int> auxiliar;

	do {
				for (int i = 0; i < n; ++i) // [0..N-1] integers
				{
					if (bitmask[i]){
						auxiliar.push_back(i);
					}
				}

				par = make_pair (auxiliar[0],auxiliar[1]);
				pairs.push_back(par);
				auxiliar.clear();
			} while (prev_permutation(bitmask.begin(), bitmask.end()));

			bitmask.clear();
}

void twoSwapShake()
{
	shuffle(pairs.begin(), pairs.end(), default_random_engine(seed));

	if(verbose == 1)
		cout<<"\ntwoSwapShake started"<<endl;

	for(int i=0; i<pairs.size()*ALPHA; i++){
		if(verbose == 1){
			for (int j = 0; j<permutation.size(); j++)
				cout<<permutation[j]<<" ";
			cout<<endl;
		}

		if(verbose == 1)
			cout<<"Pair "<<pairs[i].first<<" "<<pairs[i].second<<" swapped"<<endl;

		swap(permutation[pairs[i].first], permutation[pairs[i].second]);

		if(verbose == 1){
			for (int j = 0; j<permutation.size(); j++)
				cout<<permutation[j]<<" ";
			cout<<endl;
		}
	}
}

void twoSwapDescent()
{
	shuffle(pairs.begin(), pairs.end(), default_random_engine(seed));

	if(verbose == 1){
		cout<<"\ntwoSwapDescent started"<<endl;
	}

	for(int i=0; i<pairs.size()*BETA; i++){

		ILSt2 = high_resolution_clock::now();
		ILStimeSpan = duration_cast<duration<double> >(ILSt2 - ILSt1);
		if (ILStimeSpan.count() >= TIME_LIMIT)
			return;

		if(verbose == 1){
			for (int j = 0; j<permutation.size(); j++)
				cout<<permutation[j]<<" ";
			cout<<endl;
			for(int i = 0; i < m; i++){
					for(int j = 0; j < n; j++){
							cout << (int)inputMatrix[i][permutation[j]] << " ";
					}
					cout << endl;
			}
		}

		if(pairs[i].first != pairs[i].second-1)
		{
			if(deltaLuisExchange(pairs[i].first, pairs[i].second) < 0)
			{
				if(verbose == 1)
					cout<<"Pair "<<pairs[i].first<<" ("<<*(permutation.begin()+pairs[i].first)<<") "<<pairs[i].second<<" ("<<*(permutation.begin()+pairs[i].second)<<") swapped"<<endl;

				swap(permutation[pairs[i].first], permutation[pairs[i].second]);

				if(verbose == 1){
					for (int j = 0; j<permutation.size(); j++)
						cout<<permutation[j]<<" ";
					cout<<endl;

					for(int i = 0; i < m; i++){
							for(int j = 0; j < n; j++){
									cout << (int)inputMatrix[i][permutation[j]] << " ";
							}
							cout << endl;
					}
				}
				shuffle(pairs.begin(), pairs.end(), default_random_engine(seed));
				i=-1;
			}
		}
	}
}

void insertionDescent()
{
	shuffle(pairs.begin(), pairs.end(), default_random_engine(seed));

	if(verbose == 6)
		cout<<"\ninsertionDescent started"<<endl;

	for(int i=0; i<pairs.size()*BETA; i++){

		ILSt2 = high_resolution_clock::now();
		ILStimeSpan = duration_cast<duration<double> >(ILSt2 - ILSt1);
		if (ILStimeSpan.count() >= TIME_LIMIT)
			return;

		if(pairs[i].first != pairs[i].second-1)
			if(deltaLuisShift(pairs[i].first, pairs[i].second) < 0)
			{
				if(verbose == 6){
					cout<<" delta: "<<deltaLuisShift(pairs[i].first, pairs[i].second)<<endl;
					cout<<"\nNegative delta! "<<pairs[i].first<<" "<<pairs[i].second<<endl;
				}

				int whereFrom, columnIndex;
				vector<int>::iterator whereTo;

				columnIndex = permutation[pairs[i].first];
				whereFrom = pairs[i].first;
				whereTo = permutation.begin()+pairs[i].second+1;

				if(verbose == 6){
					for (int j = 0; j<permutation.size(); j++)
						cout<<permutation[j]<<" ";
					cout<<endl;
				}

				permutation.insert(whereTo, columnIndex);

				if(verbose == 6){
					for (int j = 0; j<permutation.size(); j++)
						cout<<permutation[j]<<" ";
					cout<<endl;
				}

				permutation.erase(permutation.begin() + whereFrom);

				if(verbose == 6){
					for (int j = 0; j<permutation.size(); j++)
						cout<<permutation[j]<<" ";
					cout<<endl;
					getchar();
				}

				shuffle(pairs.begin(), pairs.end(), default_random_engine(seed));
				i=-1;
			}
	}
}

void twoOptShake()
{
	shuffle(pairs.begin(), pairs.end(), default_random_engine(seed));

	for(int i=0; i<pairs.size()*ALPHA; i++)
		if(pairs[i].second-pairs[i].first==1)
			swap(permutation[pairs[i].first], permutation[pairs[i].second]);
		else
			reverse(permutation.begin()+pairs[i].first, permutation.begin()+pairs[i].second);
}

void twoOptDescent()
{

	int d2;




	int bestSolution = evaluation();
	int currentSolution = 0;

	shuffle(pairs.begin(), pairs.end(), default_random_engine(seed));

	if(verbose == 1){
		cout<<"\ntwoOptDescent started"<<endl;
		cout<<"Initial solution: "<<bestSolution<<endl;
		getchar();
	}

	for(int i=0; i<pairs.size()*BETA; i++){
		/*
		if(verbose == 1){
			for (int j = 0; j<permutation.size(); j++)
				cout<<permutation[j]<<" ";
			cout<<endl;
		}
		*/
		ILSt2 = high_resolution_clock::now();
		ILStimeSpan = duration_cast<duration<double> >(ILSt2 - ILSt1);
		if (ILStimeSpan.count() >= TIME_LIMIT)
			return;

		if(verbose == 1)
			cout<<"Pair "<<pairs[i].first<<" ("<<*(permutation.begin()+pairs[i].first)<<") "<<pairs[i].second<<" ("<<*(permutation.begin()+pairs[i].second)<<") reversed"<<endl;

			if(verbose == 10)
			{
					d2 = deltaLuisInversion(pairs[i].first+1, pairs[i].second+1);
			}


		if(pairs[i].second-pairs[i].first==1)
			swap(permutation[pairs[i].first], permutation[pairs[i].second]);
		else
			reverse(permutation.begin()+pairs[i].first, permutation.begin()+pairs[i].second);

		/*if(verbose == 1){
			for (int j = 0; j<permutation.size(); j++)
				cout<<permutation[j]<<" ";
			cout<<endl;
		}
		*/
		currentSolution = evaluation();

		if(verbose == 10)
			cout<<"Current solution: "<<currentSolution<<" - best solution: "<<bestSolution<<endl;

		if(currentSolution < bestSolution)
		{
			if(verbose == 10)
			{
				cout<<"delta Luis: "<<d2<<endl;
				getchar();
			}

			shuffle(pairs.begin(), pairs.end(), default_random_engine(seed));
			i=-1;
			bestSolution = currentSolution;

			if(verbose == 1){
				cout<<"New best solution: "<<bestSolution<<endl;
				//getchar();
			}
		}
		else if(currentSolution > bestSolution)
				 {
					 if(pairs[i].second-pairs[i].first==1)
					 	swap(permutation[pairs[i].first], permutation[pairs[i].second]);
					 else
					 	reverse(permutation.begin()+pairs[i].first, permutation.begin()+pairs[i].second);
						//reverse(permutation.begin()+pairs[i].first, permutation.begin()+pairs[i].second);
				 }
	}
}

void oneBlockGrouping()
{
	 vector<int> rows;
	 pair<int, int> firstBlock, secondBlock;

 	 for(int i = 0; i < m; i++)
	 		rows.push_back(i);

	 if(verbose == 1){
		 cout<<"\noneBlockGrouping started"<<endl;
		 getchar();
		 for(int i = 0; i < m; i++){
		 		for(int j = 0; j < n; j++)
					cout<<inputMatrix[i][permutation[j]]<<" ";
	 			cout<<endl;
		 }
	 }

	 shuffle (rows.begin(), rows.end(), default_random_engine(seed));

	 if(verbose == 1)
 		cout<<"rows shuffled"<<endl;

	 int k, l, o;

		for(int i = 0; i < m*GAMMA; i++){
			for (int j= 0; j < permutation.size(); j++){
				 if(inputMatrix[rows[i]][permutation[j]] > 0){
					   if(verbose == 4)
								cout<<"Random row: "<<rows[i]+1<<endl;

						 k = j;
						 while(k+1 < permutation.size() && inputMatrix[rows[i]][permutation[k+1]] > 0)
								 k++;

						 //salva os dados do primeiro block
						 firstBlock = make_pair(j, k);

						 if(verbose == 4)
	 					 		cout<<"Found the first 1-block between columns "<<firstBlock.first+1<<" and "<<firstBlock.second+1<<endl;

							//a partir do primeiro block encontra o segundo
						 for(l = k+1; l < permutation.size(); l++){
								 if(inputMatrix[rows[i]][permutation[l]] >= 1  && l < permutation.size()){
										 o = l;
										 while(o+1 < permutation.size() && inputMatrix[rows[i]][permutation[o+1]] >= 1 )
												 o++;
										 break;
								 }
						 }

						 //se existir o segundo block, salva suas informacoes
						 if(l != permutation.size()){

								 secondBlock = make_pair(l, o);
								 if(verbose == 4)
										cout<<"Found the second 1-block between columns "<<secondBlock.first+1<<" and "<<secondBlock.second+1<<endl;
						 }
						 else
						 {
								 secondBlock = make_pair(-1, -1);
								 if(verbose == 1)
										cout<<"Did not find a second 1-block"<<endl;
						 }

						 //neste momento os dois one-blocks foram identificados.
						 //se dois block foram encontrados chamar a funcao de troca
						 if(secondBlock.first != -1){
							 	int whereFrom, columnIndex;
								vector<int>::iterator whereTo;

								whereFrom = firstBlock.first;

 						 		for(int inner = firstBlock.first; inner <=firstBlock.second; inner++)
								{
									if(verbose == 4){
										 cout<<"Now trying column in position"<<whereFrom+1<<endl;
									}

									ILSt2 = high_resolution_clock::now();
									ILStimeSpan = duration_cast<duration<double> >(ILSt2 - ILSt1);
									if (ILStimeSpan.count() >= TIME_LIMIT)
										return;

									int columnIndex = permutation[whereFrom];

									std::uniform_real_distribution<double> realDistribution(0.0,1.0);

									int offset;

									int d1, d2, delta;

									delta = 0;
									d1 = deltaLuisShift (whereFrom, secondBlock.first);
									d2 = deltaLuisShift (whereFrom, secondBlock.second+1);

									if(d1 < d2)
									{
										whereTo = permutation.begin()+secondBlock.first;			//if the insertion is BEFORE the second 1-block
										offset = secondBlock.first;
										delta = d1;
									}
									else if (d2 < d1)
									{
										whereTo = permutation.begin()+secondBlock.second+1;		//if the insertion is AFTER the second 1-block
										offset = secondBlock.second+1;
										delta = d2;
									}
									if(verbose == 4)
									{
											cout<<"Trying to move the "<<whereFrom+1<<"st/nd/rd/th column ("<<columnIndex<<") to near the element "<<*whereTo<<endl;
											for(int k = 0; k<permutation.size(); k++)
												cout<<permutation[k]<<" ";
											cout<<endl;
									}

									if (delta<0)
									{
									 	permutation.insert(whereTo, columnIndex);
										permutation.erase(permutation.begin() + whereFrom);

										if(verbose == 4)
										{
												cout<<"New solution"<<endl;
												for(int k = 0; k<permutation.size(); k++)
													cout<<permutation[k]<<" ";
												cout<<endl;
										}

										j = secondBlock.first-1;
										while(inputMatrix[rows[i]][permutation[j]] > 0)
												 j--;
							    }
								}
						 }
						 else
						 {
							 	 //caso nao exista o contador recebe o final do primeiro block
							 	 if(verbose == 1){
										cout<<"Moving forward..."<<endl;
							 	 }
								 j=k;
						 }
					}
			 }
	 }
}

int simulatedAnnealing()
{
	int i;
	int bestSolution = evaluation();
	int currentSolution = bestSolution;
	int newSolution;
	float T = -0.60*bestSolution/log(0.5);
	float Tend = -0.45*bestSolution/log(0.5);
	float cooling = 0.995;
	convergence = 0;

	const double e = 2.718281828;

	vector<int> backupBest(permutation);
	vector<int> backupCurrent(permutation);

	std::uniform_real_distribution<double> realDistribution(0.0,1.0);

	generatePairs();

	if(verbose == 1){
		cout<<"SA started"<<endl;
		cout<<"Initial solution: "<<bestSolution<<endl;
	}

	for(i=0; i<MAX_ITER && T > 0.01; i++)
	{
		backupCurrent = permutation;
		if(shaking == 1)
			twoOptShake();
		else
				twoSwapShake();

		switch(localSearch){
			case 1: twoOptDescent(); break;
			case 2: twoSwapDescent(); break;
			case 3: oneBlockGrouping(); break;
			case 4: twoOptDescent(); oneBlockGrouping(); break;
			case 5: twoSwapDescent(); oneBlockGrouping(); break;
			default: twoSwapDescent(); oneBlockGrouping(); break;
		}

		newSolution = evaluation();

		if(verbose == 1){
			cout<<"New solution: "<<newSolution<<" - current solution: "<<currentSolution<<" - best solution: "<<bestSolution<<endl;
		}

		if(newSolution < currentSolution)
		{
			currentSolution = newSolution;
			backupCurrent = permutation;

			if(verbose == 1){
				cout<<"New current solution: "<<currentSolution<<endl;
			}

			if(newSolution < bestSolution)
			{
				bestSolution = newSolution;
				backupBest = permutation;
				convergence = i;

				if(verbose == 1){
					cout<<"New best solution: "<<bestSolution<<endl;
				}
			}

			if (bestSolution == lowerBound)
				break;
		}
		else
		{		if(realDistribution(generator) <= pow(e, - (newSolution - bestSolution) / T))
				{
					if(verbose == 1){
						cout<<"New current solution: "<<currentSolution<<endl;
					}
				}
				else
					permutation = backupCurrent;
		}

		T=T*pow(-0.41*bestSolution/log(0.5)/-0.45*bestSolution/log(0.5), 1/i);
	}

	permutation = backupBest;
	return bestSolution;
}

int iteratedLocalSearch()
{
	int i;
	int bestSolution = evaluation();
	int currentSolution = 0;
	convergence = 0;
	vector<int> backup(permutation);

	generatePairs();

	if(verbose == 1){
		cout<<"ILS started"<<endl;
		cout<<"Initial solution: "<<bestSolution<<endl;
		cout<<"solution size: "<<n<<endl;
		getchar();
	}

	ILSt1 = high_resolution_clock::now();

	for(i=0; i<MAX_ITER; i++)
	{
		backup = permutation;
		if(shaking == 1)
			twoOptShake();
		else
				twoSwapShake();

		switch(localSearch){
			case 1: twoOptDescent(); break;
			case 2: twoSwapDescent(); break;
			case 3: oneBlockGrouping(); break;
			case 4: twoOptDescent(); oneBlockGrouping(); break;
			case 5: twoSwapDescent(); oneBlockGrouping();insertionDescent(); break;
			default: oneBlockGrouping(); break;
		}

		currentSolution = evaluation();

		if(verbose == 1){
			cout<<"Current solution: "<<currentSolution<<" - best solution: "<<bestSolution<<endl;
		}

		if(currentSolution > bestSolution)
		{
			permutation = backup;
		}
		else
		{
			bestSolution = currentSolution;
			convergence = i;
			if(verbose == 1){
				cout<<"New best solution: "<<bestSolution<<endl;
		 }
		 if (bestSolution == lowerBound)
				break;
		}

		ILSt2 = high_resolution_clock::now();

		ILStimeSpan = duration_cast<duration<double> >(ILSt2 - ILSt1);
		if (ILStimeSpan.count() >= TIME_LIMIT)
			return bestSolution;
	}
	return bestSolution;
}

/*
Main procedure
*/

void singleRun(int& initialSolutionValue, int& finalSolutionValue, double& initialRunningTime, double& finalRunningTime, string inputFileName, int run)
{
	int i, j;
	int FinalSolutionValue = INT_MAX;
	initialSolutionValue = 0;
	finalSolutionValue = 0;
	initialRunningTime = 0.0;
	finalRunningTime = 0.0;

	if(verbose == 1)
		cout<<"\nProcedure started"<<endl;

	readProblem(inputFileName);

	if(verbose == 1)
	{
		cout<<"Problem read"<<endl;
		cout<<"m: "<<m<<" n: "<<n<<endl;
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	setLowerBound();

	if(verbose == 1)
		cout<<"Lower bound set"<<endl;

	if(pp==1 || pp==3)
		preProcessingRepeated();

	if(pp==2 || pp==3)
		preProcessingDominated();

	if(verbose == 1)
		cout<<"Preprocessing performed"<<endl;

	high_resolution_clock::time_point t3 = high_resolution_clock::now();
	buildGraph();

	if(verbose == 1)
		cout<<"Built graph"<<endl;

	breadthFirstSearch();

	if(verbose == 1)
		cout<<"BFS performed"<<endl;

	columnSequencing();

	if(verbose == 1)
		cout<<"Columns sequenced"<<endl;

	initialSolutionValue = evaluation();

	if(verbose == 1)
		cout<<"Initial solution evaluated"<<endl;

	high_resolution_clock::time_point t4 = high_resolution_clock::now();

	if(mh == 1)
		finalSolutionValue = iteratedLocalSearch();
	else
		finalSolutionValue = simulatedAnnealing();

	if(verbose == 1)
		cout<<"Iterated Local Search performed"<<endl;

	if(pp!=0){
		addDominatedColumns();
		if(verbose == 1)
			cout<<"Dominated columns added to the solution"<<endl;
	}

  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  duration<double> finalTimeSpan = duration_cast<duration<double> >(t2 - t1);
	duration<double> initialTimeSpan = duration_cast<duration<double> >(t4 - t3);

	initialRunningTime = initialTimeSpan.count();
	finalRunningTime =  finalTimeSpan.count();

	printSolution(inputFileName, initialSolutionValue, finalSolutionValue, initialRunningTime, finalRunningTime, run);

	if (!checkPermutation())
	{
		cout<<"\n\nERROR! Solution is not a permutation of "<<nOriginal<<" elements"<<endl;
		getchar();
	}

	termination();
}
