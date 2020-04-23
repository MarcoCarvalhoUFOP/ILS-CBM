#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <cstdio>
#include <algorithm>
using namespace std;

int n, m;
vector< vector<int> > inputMatrix;

//int myints[] = {2, 87, 125, 35, 124, 181, 73, 17, 165, 25, 108, 99, 174, 4, 88, 195, 140, 154, 89, 65, 44, 122, 34, 137, 20, 153, 8, 53, 189, 12, 151, 40, 67, 118, 192, 51, 196, 186, 150, 39, 185, 143, 32, 14, 162, 82, 113, 36, 11, 24, 27, 141, 120, 159, 167, 98, 111, 83, 112, 135, 145, 166, 47, 37, 101, 26, 85, 62, 0, 31, 121, 90, 134, 190, 68, 117, 63, 3, 128, 22, 76, 5, 75, 50, 149, 123, 86, 23, 199, 79, 102, 180, 92, 115, 28, 146, 104, 97, 93, 59, 45, 139, 164, 61, 173, 21, 15, 84, 18, 136, 177, 147, 74, 46, 116, 130, 94, 126, 163, 169, 43, 58, 198, 78, 72, 129, 152, 41, 29, 138, 184, 38, 70, 172, 182, 77, 6, 100, 193, 56, 170, 178, 188, 175, 95, 144, 197, 183, 142, 160, 156, 55, 155, 1, 119, 9, 161, 171, 81, 148, 168, 30, 158, 109, 110, 157, 91, 33, 103, 10, 131, 187, 191, 176, 16, 57, 107, 114, 194, 71, 69, 179, 106, 19, 133, 49, 42, 7, 105, 64, 13, 127, 52, 60, 54, 80, 48, 96, 132, 66};

vector< int > permutation;//(myints, myints + sizeof(myints) / sizeof(int) );
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
	if (fpIn == NULL) perror ("Error opening file");

	fscanf(fpIn, "%d %d", &m, &n);

	inputMatrix.resize(m);

	for (int i = 0; i < m; i++) {
			inputMatrix[i].resize(n);
	}

	for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
					fscanf(fpIn, "%d", &inputMatrix[i][j]);
					if (inputMatrix[i][j] != 0) {
							inputMatrix[i][j] = 1;
					}
			}
	}
	fclose(fpIn);
}

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

	cout<<"Solution value: "<<value<<endl;
	return value;
}

bool checkPermutation()
{
		sort(permutation.begin(), permutation.end());
		for (int i= 0; i<permutation.size(); i++)
			if(permutation[i]!=i)
				return false;

		return true;
}


int main(int argc, char* argv[])
{
	int element;
	string inputFileName, instance;

	ifstream permFile;

	cout<<"Enter file name: "<<endl;
	cin>>inputFileName;

	permFile.open(inputFileName, ifstream::in);

	if(!permFile.is_open())
		cout<<"\n\nERROR! File not opened!"<<endl;

	permFile>>instance;

	readProblem(instance);

	for(int i=0; i<n; i++)
	{
		permFile>>element;
		permutation.push_back(element);
	}

	evaluation();

	if (!checkPermutation())
		cout<<"\n\nERROR! Solution is not a permutation of "<<n<<" elements"<<endl;

	permFile.close();

	return 0;
}
