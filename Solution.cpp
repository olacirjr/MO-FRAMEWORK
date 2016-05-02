//class solution
class Solution{
public:
	Problem problem;
	
	Solution();//constructor
	~Solution();//destructor
	Solution(const Solution &source);//copy
	Solution& operator= (const Solution &source);//assignment
	
	double *decisionVector;
	double *objectiveVector;
	
	double crowdingDistance;
	double weightedDistance;// stores the weightedDistance relative to the swarm
	bool dominated;
	bool evalSolution;//only evaluate this solution if this flag is true

	//compare two Solutions
	//param - the solution to be compared to
	bool isEqual(const Solution &sol) const;
	
	double localOptimizer(double* weight);

	//evaluate this solution
	void evaluate();
	
	//prints the decison and objective vectors of a solution
	void print();
};

// // // constructor randomly initialize the variables
Solution::Solution(){
	crowdingDistance = -1; //if the crowding distance is negative, it is not set
	dominated=true;
	evalSolution=true;

	decisionVector = new double[decisionNumber];
	objectiveVector = new double[objectiveNumber];
	
	if(!strcmp(algorithm, "cmaes-mopso")){//if it is CMAES-MOPSO, start the solution with a gaussian
		for(int i=0;i<decisionNumber;i++){
			double mean=(rand()/(double)RAND_MAX);
			double sigma=0.3;//as in the CMA-ES code
			double normalizedValue= mean+sigma*cmaes_random_Gauss();
			while(normalizedValue > 1 || normalizedValue < 0)//resample until is within the bounds
				normalizedValue= mean+sigma*cmaes_random_Gauss();
			
// 			normalizedValue=std::min(normalizedValue,1.0);//limit to the maximum value of 1
// 			normalizedValue=std::max(normalizedValue,0.0);//limit to the minimum value of 0
			
			decisionVector[i]=inferiorPositionLimit[i]+ (normalizedValue* (superiorPositionLimit[i]-inferiorPositionLimit[i]) );
		}
	}else{
		for(int i=0;i<decisionNumber;i++)
			decisionVector[i]=inferiorPositionLimit[i]+ ((rand()/(double)RAND_MAX)* (superiorPositionLimit[i]-inferiorPositionLimit[i]) );
	}
// 	problem.evaluate(decisionVector, objectiveVector);
}

//Destructor
Solution::~Solution(){
	delete[] objectiveVector;
	delete[] decisionVector;
}

//Copy constructor
Solution::Solution(const Solution &source){
// 	if (this != &source) {
// 		if (objectiveVector != NULL) {
// // 			delete[] objectiveVector;
// 		}
// 	}
	objectiveVector = new double[objectiveNumber];
	decisionVector = new double[decisionNumber];
	
	crowdingDistance=source.crowdingDistance;
	weightedDistance=source.weightedDistance;
	dominated=source.dominated;
	evalSolution=source.evalSolution;
	
	memcpy(objectiveVector, source.objectiveVector, sizeof(double)*objectiveNumber);
	memcpy(decisionVector, source.decisionVector, sizeof(double)*decisionNumber);
}

//Assignment operator
Solution& Solution::operator= (const Solution &source){
	crowdingDistance=source.crowdingDistance;
	weightedDistance=source.weightedDistance;
	dominated=source.dominated;
	evalSolution=source.evalSolution;
	
	memcpy(objectiveVector, source.objectiveVector, sizeof(double)*objectiveNumber);
	memcpy(decisionVector, source.decisionVector, sizeof(double)*decisionNumber);
}

//evaluate the solution according to a predetermined problem
void Solution::evaluate(){
	if(evalSolution){
		problem.evaluate(decisionVector, objectiveVector);
		crowdingDistance= -1;
		weightedDistance=-1;
		for(int o=0;o<objectiveNumber;o++){
			maxObjectives[o]=std::max(maxObjectives[o], objectiveVector[o]);
			minObjectives[o]=std::min(minObjectives[o], objectiveVector[o]);
		}
	}
}

//comparator if two solutions are equal
//param - the solution to be compared to
bool Solution::isEqual(const Solution &sol) const{
	for(int i=0;i<objectiveNumber;i++){
		//if( sol.objectiveVector[i] != objectiveVector[i] ){
		if( round_5(sol.objectiveVector[i]) != round_5(objectiveVector[i]) ){
			return false;
		}
	}
	return true;
}

void Solution::print(){
	printVector(objectiveVector,objectiveNumber);
	printf("-> ");
	printVector(decisionVector,decisionNumber);
	printf("\n");
}

double Solution::localOptimizer(double* weight){
	double  rhobeg = 0.5, rhoend = 0.0001;
	int iprint=0, maxfun = 1000;
	example_state state;
	state.nprob = 0;
	
// 	if(objectiveVector[0]<2000){
		globalTmpWeight=weight;
		int rc = cobyla(decisionNumber, 2*decisionNumber, decisionVector, rhobeg, rhoend, iprint, &maxfun, localOptimizerObjFunc, &state);
// 	}
}