#include <cstdio>
#include <climits>
using namespace std;
#include "solver.hpp"

int main(int argc, char* argv[])
{
	double averageInitialSolutionValue;
	double averageFinalSolutionValue;
	double bestSolutionValue = INT_MAX;
	double averageRunningTime;
	double averageConvergence;

	int finalSolutionValue;
	int initialSolutionValue;
	int i = 0, j;

	double initialRunningTime;
	double finalRunningTime;

	string inputFileName;

	verbose = 0;
	pp = 0;
	searchRule = 3;
	initialNode = 3;
	localSearch = 4;
	shaking = 2;
	mh = 1;

	vector<string> arguments(argv + 1, argv + argc);
	string::size_type sz;

	for(int i=0; i<arguments.size(); i+=2)
	{
		if(arguments[i]== "-v")
					verbose = stoi(arguments[i+1], &sz);								//tracing level (1) all (2) all solution values (3) only best values
			else if(arguments[i]== "-pp")
							pp = stoi(arguments[i+1], &sz);									//preprocessing: (0) none (1) equal sizes (2) smaller and larger (3) both
				else if(arguments[i]== "-s")
							searchRule = stoi(arguments[i+1], &sz); 				//search rule (1) minimum degree (2) maximum degree (3) maximum edge weight
					else if(arguments[i]== "-i")
								initialNode = stoi(arguments[i+1], &sz);			//BFS initial node (1) minimum degree (2) maximum degree (3) random
								else if(arguments[i]== "-alpha")
											ALPHA = stod(arguments[i+1], &sz);					//percentage of ILS perturbation
											else if(arguments[i]== "-beta")
														BETA = stod(arguments[i+1], &sz);			//percentage of ILS local search
														else if(arguments[i]== "-it")
																	MAX_ITER = stoi(arguments[i+1], &sz);	//maximum number of ILS iterations
																		else if(arguments[i]== "-gamma")
																						GAMMA = stoi(arguments[i+1], &sz);	//percentage of rows analyzied during 1-block local search
																				 else if(arguments[i]== "-ls")
																								localSearch = stoi(arguments[i+1], &sz);	//choice of local search (1) 2-opt (2) 2-swap (3) 1-block grouping (4) 1+3 (5) 2+3
																							else if(arguments[i]== "-p")
			 																								shaking = stoi(arguments[i+1], &sz);	//choice of shaking (1) 2-opt (2) 2-swap
																									 else if(arguments[i]== "-mh")
							 																								mh = stoi(arguments[i+1], &sz);	//choice of metaheuristic (1) ILS (2) SA
																												else if(arguments[i]== "-r")
	 							 																								RUNS = stoi(arguments[i+1], &sz);	//number of runs per instance
	}

	ifstream fpIndex("index.txt");
	ofstream fpOut("RESULTS_SUMMARY.txt");

	if(verbose == 1)
	{
			cout<<"Algorithm settings:"<<endl;
			cout<<"Verbose mode on"<<endl;
			switch(pp){
			 case 1:	cout<<"Preprocessing repeated columns"<<endl; break;
			 case 2:	cout<<"Preprocessing dominated columns"<<endl;  break;
			 case 3: 	cout<<"Preprocessing repeated and dominated columns"<<endl;  break;
			 default: cout<<"No preprocessing"<<endl;
		  }
			switch(searchRule){
				case 1:	cout<<"Search guided by minimum degree"<<endl; break;
				case 2:	cout<<"Search guided by maximum degree"<<endl;  break;
				case 3: cout<<"Search guided by maximum edge weight"<<endl;  break;
				default: cout<<"No search guidance"<<endl;
		  }
	    switch(initialNode){
			 case 1:	cout<<"BFS started by minimum degree node"<<endl; break;
			 case 2:	cout<<"BFS started by maximum degree node"<<endl;  break;
			 case 3: 	cout<<"BFS started by random node"<<endl;  break;
			 default: cout<<"No choice of initial node"<<endl;
	    }
	    switch(shaking){
			 case 1:	cout<<"Shaking: 2-opt"<<endl; break;
		 	 case 2:	cout<<"Shaking: 2-swap"<<endl; break;
			 default: cout<<"No choice of shaking"<<endl;
	    }
	 	  switch(localSearch){
			 case 1:	cout<<"Local Search: 2-opt"<<endl; break;
			 case 2:	cout<<"Local Search: 2-swap"<<endl; break;
			 case 3:  cout<<"Local Search: 1-block grouping"<<endl; break;
			 case 4:  cout<<"Local Search: 2-opt and 1-block grouping"<<endl; break;
			 case 5:  cout<<"Local Search: 2-swap and 1-block grouping"<<endl; break;
			 default: cout<<"No choice of local search"<<endl;
		 }
		 switch(mh){
			case 1:		cout<<"Metaheuristic Search: ILS"<<endl; break;
			case 2:		cout<<"Metaheuristic Search: SA"<<endl; break;
			default: 	cout<<"No choice of Metaheuristic Search"<<endl;
		 }

		 cout<<"Alpha: "<<ALPHA<<endl;
		 cout<<"Beta: "<<BETA<<endl;
		 cout<<"Gamma: "<<GAMMA<<endl;
		 cout<<"Max iterations: "<<MAX_ITER<<endl;
	}

	while(fpIndex>>inputFileName)
	{
		averageInitialSolutionValue = 0.0;
		averageFinalSolutionValue = 0.0;
		bestSolutionValue = INT_MAX;
		averageRunningTime = 0.0;
		averageConvergence = 0.0;

		i++;

		for(j=0; j<RUNS; j++)
		{
			singleRun(initialSolutionValue, finalSolutionValue, initialRunningTime, finalRunningTime, inputFileName, j+1);

			if(verbose == 2)
				printf("RUN %d/%d - PROBLEM %d: %s %d %d %lf %lf\n", j+1, RUNS, i, inputFileName.c_str(), initialSolutionValue, finalSolutionValue, initialRunningTime, finalRunningTime);
			fpOut<<inputFileName<<" "<<j+1<<" "<<initialSolutionValue<<" "<<initialRunningTime<<" "<<finalSolutionValue<<" "<<finalRunningTime<<endl;

			averageInitialSolutionValue += initialSolutionValue;
			averageFinalSolutionValue += finalSolutionValue;
			averageRunningTime += finalRunningTime;
			averageConvergence += convergence;

			if(finalSolutionValue < bestSolutionValue)
				bestSolutionValue = finalSolutionValue;
		}
		if(RUNS > 1)
		{
			averageInitialSolutionValue/=RUNS;
			averageFinalSolutionValue/=RUNS;
			averageRunningTime/=RUNS;
			averageConvergence/=RUNS;

			if(verbose == 2)
			{
				cout<<"\nAVERAGE RESULTS: Initial Solution: "<<averageInitialSolutionValue<<" Final Solution: "<< averageFinalSolutionValue<<" Running Time: "<<averageRunningTime<<" Convergence: "<<averageConvergence<<endl;
				cout<<"Best solution value: "<<bestSolutionValue<<"\n"<<endl;
			}
			else if(verbose == 3)
							cout<<"PROBLEM " <<i<<" "<<inputFileName.c_str()<<" best solution value: "<<bestSolutionValue<<endl;

			fpOut<<"\nAVERAGE RESULTS: Initial Solution: "<<averageInitialSolutionValue<<" Final Solution: "<< averageFinalSolutionValue<<" Running Time: "<<averageRunningTime<<" Convergence: "<<averageConvergence<<endl;
			fpOut<<"Best solution value: "<<bestSolutionValue<<endl;
			fpOut<<endl;
		}
	}
}
