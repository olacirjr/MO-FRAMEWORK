class Particle{
public:

	Solution solution; //solution contained by this particle
	double *velocity; //velocity of this particle
	Solution localBest; //local best solution
	Solution globalBest; //global best solution
	char leaderType[50];
	
	double* weight; //pointer to the weight vector of the repository
	
	Particle();
	~Particle();
	
// 	operator Solution(){return solution; }

	//initialize the particle with random values
	void initialize();
	//special initialization of particles with the local best being initialized with a ramdom solution from a given set
	void initialize(Solution* candidateSet, int candidateNumber);
	//compute the speed of the particle
	void computeSpeed();
	//update the position according to the velocity
	void updatePosition();
	//choose the global leader
	//params - a set of solutions representing the actual pareto front
	//		 - the code of the strategy for selecting the leader (set in the main.cpp)
	//		 - the number of solutions composing the given front
	//		 - the current status of the random number generator
	void chooseGlobalLeader(Solution* repository, int repActualSize);
	//strategy to choose the global leader - chooses the leader randomly
	//params - a set of solutions representing the actual pareto front
	//		 - the number of solutions composing the given front
	//		 - the current status of the random number generator
	void randomLeader(Solution* repository, int repActualSize);
	//strategy to choose the global leader - chooses the leader based in the crowding distance
	//params - a set of solutions representing the actual pareto front
	//		 - the number of solutions composing the given front
	//		 - the current status of the random number generator
	void crowdingDistanceLeader(Solution* repository, int repActualSize);
	//The method NWSum of leader selection
	//params - a set of solutions representing the actual pareto front
	//		 - the number of solutions composing the given front
	void NWSumLeader(Solution* repository, int repActualSize);
	//The method Sigma of leader selection
	//params - a set of solutions representing the actual pareto front
	//		 - the number of solutions composing the given front
	void SigmaLeader(Solution* repository, int repActualSize);
	//update the local best solution
	void updateLocalLeader();
	//reduce the speed if the particle exceed the bounds
	//void reduceSpeed();
	//the turbulence operator
	void turbulence();
	//keep the position in the swarm specified bounds in IMulti algorithm
	int truncatePositionIMulti(double* centroid, double range, double* smallerDecision, double* largerDecision);
	
	//void truncatePositon();

private:

	double c1;
	double c2;

	//initialize the velocity randomly
	void initializeVelocity();
	//update omega value randomly
	double getOmega();
	//update phi value randomly
	double getPhi();
	//update c1 value randomly
	double getC1();
	//update C2 value randomly
	double getC2();
	//calculate Fi value
	double getFi();
	
	Particle(const Particle &source){}//copy
	Particle& operator= (const Particle &source){}//assignment
};

//Constructor
Particle::Particle(){
	memset(solution.objectiveVector, 0, sizeof(double)*objectiveNumber);
	velocity = new double[decisionNumber];
}

//destructor
Particle::~Particle(){
	delete[] velocity;
}

//initialize the particle with random values
void Particle::initialize(){
	solution.evaluate();
	
// 	globalBest=localBest=solution;//this does not work, do not know why
	globalBest=solution;
	localBest=solution;

	initializeVelocity(); //initialize the velocity with random values in 0...0
}
//initialize the particle with random values
void Particle::initialize(Solution* candidateSet, int candidateNumber){
	int rnd=rand() % candidateNumber;
	
// 	globalBest=localBest=solution; //this does not work, do not know why
	globalBest=solution;
	localBest=solution;
	
	initializeVelocity(); //initialize the velocity with random values in 0...0
}

double velocityConstriction(double v, int var){
	double result=v;
	double dmax = deltaMax[var];
	double dmin = deltaMin[var];
    
    
    if (v > dmax) {
      result = dmax;
    }
    if (v < dmin) {
      result = dmin;
    }
    return result;
  } // velocityConstriction
  
    // constriction coefficient (M. Clerc)
double constrictionCoefficient(double c1, double c2) {
    double rho = c1 + c2;
    //rho = 1.0 ;
    if (rho <= 4) {
      return 1.0;
    } else {
      return 2.0 / (2.0 - rho - sqrt( (rho*rho) - 4.0 * rho));
    }
  } // constrictionCoefficient
  
//compute the speed of the particle
void Particle::computeSpeed(){
	double r1 = getPhi();
	double r2 = getPhi();
	double C1 = getC1();
	double C2 = getC2();
	//	
	for (int var = 0; var < decisionNumber; var++) {
	//Computing the velocity of this particle
		velocity[var]=velocityConstriction( constrictionCoefficient(C1, C2) * (INERTIA * velocity[var] +
		C1 * r1 * (localBest.decisionVector[var] - solution.decisionVector[var]) +
		C2 * r2 * (globalBest.decisionVector[var]- solution.decisionVector[var])), var);
	}
}

//update the particle position according to the velocity
void Particle::updatePosition(){
	for(int i=0;i<decisionNumber;i++){
		solution.decisionVector[i]+=velocity[i];
		
		//original SMPSO truncation
// 		if (solution.decisionVector[i] < inferiorPositionLimit[i]) {
// 			solution.decisionVector[i]=inferiorPositionLimit[i];
// 			velocity[i] = velocity[i] * -1; //
// 		}
// 		if (solution.decisionVector[i] > superiorPositionLimit[i]) {
// 			solution.decisionVector[i]=superiorPositionLimit[i];
// 			velocity[i] = velocity[i] * -1; //
// 		}
	}
}

//initialize the velocity
void Particle::initializeVelocity(){
	for(int i=0;i<decisionNumber;i++)
		//velocity[i]=(rand()/(double)RAND_MAX);
		velocity[i]=0.0; //JMetal
}

// //update omega value randomly
// double Particle::getOmega(){
// 	//int idx = blockIdx.x*blockDim.x + threadIdx.x;
// 	//curandState localState = devStates[idx];
// 	return (double)((rand()/(double)RAND_MAX)*MAX_OMEGA);
// 	//double value= (double)(curand_uniform(&localState)*MAX_OMEGA);
// 	//devStates[idx] = localState;
// 	//return value;
// }

//update phi value randomly
double Particle::getPhi(){
	//int idx = blockIdx.x*blockDim.x + threadIdx.x;
	//curandState localState = devStates[idx];
	//return fmod((rand()/(double)RAND_MAX),MAX_PHI);
	//double value= (double)(curand_uniform(&localState)*MAX_PHI);
	//devStates[idx] = localState;
	//return value;
	return (double)((rand()/(double)RAND_MAX)*MAX_PHI);
}

//update c1 value randomly
double Particle::getC1(){
	//int idx = blockIdx.x*blockDim.x + threadIdx.x;
	//curandState localState = devStates[idx];
	//return c1=(rand()/(double)RAND_MAX)+1.5;
	//c1=(curand_uniform( &localState)+1.5);
	//devStates[idx] = localState;
	//return c1;
	return (rand()/(double)RAND_MAX)+1.5;
}

//update c2 value randomly
double Particle::getC2(){
// 	int idx = blockIdx.x*blockDim.x + threadIdx.x;
// 	curandState localState = devStates[idx];
// 	//return c2=(rand()/(double)RAND_MAX)+1.5;
// 	c2=(curand_uniform( &localState)+1.5);
// 	devStates[idx] = localState;
// 	return c2;
	return (rand()/(double)RAND_MAX)+1.5;
}

//calculate Fi value
double Particle::getFi(){
	if(c1+c2 > 4)
		return c1+c2;
	else
		return 1;
}
//choose the global leader randomly
//params - a set of solutions representing the actual pareto front
//		 - the number of solutions composing the given front
//		 - the current status of the random number generator
void Particle::randomLeader(Solution* repository, int repActualSize){
	int idx=(rand() % repActualSize);
	
	globalBest=repository[idx];
// 	memcpy(&globalBest, &repository[idx], sizeof(Solution));
}
//choose the global leader according to the crowding distance
//params - a set of solutions representing the actual pareto front
//		 - the number of solutions composing the given front
//		 - the current status of the random number generator
void Particle::crowdingDistanceLeader(Solution* repository, int repActualSize){
	//updateCrowdingDistances(repository, repActualSize);

	//int idx = blockIdx.x*blockDim.x + threadIdx.x;
	//curandState localState = devStates[idx];
	
	int rnd1=(rand() % repActualSize);
	int rnd2=(rand() % repActualSize);
	
	//int rnd1=(int)(curand_uniform( &localState)*repActualSize);
	//int rnd2=(int)(curand_uniform( &localState)*repActualSize);

	//devStates[idx] = localState;

	if(repository[rnd1].crowdingDistance <0 || repository[rnd2].crowdingDistance<0)
		fprintf(stderr,"\nCROWDING DISTANCE LEADER ERROR!!\n");

	if(repository[rnd1].crowdingDistance > repository[rnd2].crowdingDistance){//the solution in the less crowded region (bigger crowding distance) is selected as leader
		globalBest=repository[rnd1];
// 		memcpy(&globalBest, &repository[rnd1], sizeof(Solution));
// 		for(int o=0;o<objectiveNumber;o++)
// 			globalBest.objectiveVector[o]=repository[rnd1].objectiveVector[o];
// 		for(int d=0;d<decisionNumber;d++)
// 			globalBest.decisionVector[d]=repository[rnd1].decisionVector[d];
		
	}else{
		globalBest=repository[rnd2];
// 		memcpy(&globalBest, &repository[rnd2], sizeof(Solution));
// 		for(int o=0;o<objectiveNumber;o++)
// 			globalBest.objectiveVector[o]=repository[rnd2].objectiveVector[o];
// 		for(int d=0;d<decisionNumber;d++)
// 			globalBest.decisionVector[d]=repository[rnd2].decisionVector[d];
	}
}
//The method NWSum of leader selection
//params - a set of solutions representing the actual pareto front
//		 - the number of solutions composing the given front
void Particle::NWSumLeader(Solution* repository, int repActualSize){
	double bestNWSum = -1;
// 	Solution *GBest;
	int indBest=-1;
	
	for(int p=0;p<repActualSize;p++) {
		//calculate the NWSum metric according to the actual solution and the solutions from the repository
		double NWSum= calculateNWSum(solution, repository[p]); //in util.cpp
		//The higher the NWSum the better
		//if this NWSum is higher than the previous, replace
		if(NWSum>bestNWSum){
			bestNWSum= NWSum;
// 			GBest = &repository[p];
			indBest=p;
		}
	}
	if(indBest == -1){
		fprintf(stderr, "Problem on NWSum leader %f\n", bestNWSum);
		exit(1);
	}
	
	globalBest=repository[indBest];
// 	memcpy(&globalBest, GBest, sizeof(Solution));
// 	for(int o=0;o<objectiveNumber;o++)
// 		globalBest.objectiveVector[o]=repository[indBest].objectiveVector[o];
// 	for(int d=0;d<decisionNumber;d++)
// 		globalBest.decisionVector[d]=repository[indBest].decisionVector[d];
}
//The method Sigma of leader selection
//params - a set of solutions representing the actual pareto front
//		 - the number of solutions composing the given front
void Particle::SigmaLeader(Solution* repository, int repActualSize){
	double bestSigma = MAXDOUBLE;
	int indBest=-1;
	int sigmaSize=(int)combination(objectiveNumber, 2);
	double sol[objectiveNumber];
	double solSigmaVector[sigmaSize];
	double repSol[objectiveNumber];
	double repSolSigmaVector[sigmaSize];
	
	normalizeObjectives(sol, solution.objectiveVector);
	calculateSigmaVector(solSigmaVector, sol);
	
	for(int p=0;p<repActualSize;p++) {
		normalizeObjectives(repSol, repository[p].objectiveVector);
		calculateSigmaVector(repSolSigmaVector, repSol);
		
		double sigmaDistance = calculateEuclideanDistance(solSigmaVector, repSolSigmaVector, sigmaSize);
		
		//The smaller the SigmaDistance the better
		if(sigmaDistance<bestSigma){
			bestSigma= sigmaDistance;
			indBest=p;
		}
	}
	if(indBest == -1){
		fprintf(stderr, "Problem on Sigma leader\n");
		exit(1);
	}
	globalBest=repository[indBest];
// 	memcpy(&globalBest, GBest, sizeof(Solution));
// 	for(int o=0;o<objectiveNumber;o++)
// 		globalBest.objectiveVector[o]=repository[indBest].objectiveVector[o];
// 	for(int d=0;d<decisionNumber;d++)
// 		globalBest.decisionVector[d]=repository[indBest].decisionVector[d];
}
//update the local leader
void Particle::updateLocalLeader(){
	if(!decomposition){//if it is not decomposition, update the local best using dominance
		//verify the dominance relation between
		//return 1 if sol1 dominates sol2, -1 if sol2 dominates sol1, 0 if they do not dominate each other, 2 if they are equal
		int tmp=dominance(solution.objectiveVector, localBest.objectiveVector, objectiveNumber);
		if(tmp == 1){//if the new position is better then the old
// 			memcpy(&localBest, &solution, sizeof(Solution));
			localBest=solution;
// 			for(int o=0;o<objectiveNumber;o++)
// 				localBest.objectiveVector[o]=solution.objectiveVector[o];
// 			for(int d=0;d<decisionNumber;d++)
// 				localBest.decisionVector[d]=solution.decisionVector[d];
		}else{
			if(tmp==0){ // if the solutions non dominate each other
				//if usado no algoritmo do andré
				//if(curand_uniform( &localState) > 0.5){ //the new solution have 50% chance to replace the old
					//localBest.deepCopy(solution); //this solution is the best of the particle since now
					localBest=solution;
// 					memcpy(&localBest, &solution, sizeof(Solution));
// 					for(int o=0;o<objectiveNumber;o++)
// 						localBest.objectiveVector[o]=solution.objectiveVector[o];
// 					for(int d=0;d<decisionNumber;d++)
// 						localBest.decisionVector[d]=solution.decisionVector[d];
				//}
			}
		}
	}else{//if it is decomposition
		//Using the weights in the local leader assignment
		
// 		if(PBI(newSol, weight) < PBI(oldSol, weight)){
		if(getScalarValue(solution.objectiveVector, weight) < getScalarValue(localBest.objectiveVector, weight))
			localBest=solution;
		
	}
}

//choose the global leader
//params - a set of solutions representing the actual pareto front
//		 - the code of the strategy for selecting the leader (set in the main.cpp)
//		 - the number of solutions composing the given front
//		 - the current status of the random number generator
void Particle::chooseGlobalLeader(Solution* repository, int repActualSize){
	if(!strcmp(leaderType, "rnd"))
		randomLeader(repository, repActualSize);
	if(!strcmp(leaderType, "cd"))
		crowdingDistanceLeader(repository, repActualSize);
	if(!strcmp(leaderType, "nwsum"))
		NWSumLeader(repository, repActualSize);
	if(!strcmp(leaderType, "sigma"))
		SigmaLeader(repository, repActualSize);
	
	if(strcmp(leaderType, "rnd") && strcmp(leaderType, "cd") && strcmp(leaderType, "nwsum") && strcmp(leaderType, "sigma")){
		fprintf(stderr,"INVALID LEADER SELECTION METHOD! (%s)\n", leaderType);
		exit(1);
	}
}
// //reduce the speed of the particles when it position exceeds the bounds
// void Particle::reduceSpeed(){
// 	for(int i=0;i<decisionNumber;i++)
// 		velocity[i] *= VELOCITY_REDUCTION;
// }
// //apply the turbulence operator - Andre
// void Particle::turbulence(){
// 	int idx = blockIdx.x*blockDim.x + threadIdx.x;
// 	curandState localState = devStates[idx];
// 	for(int i=0;i<decisionNumber;i++){
// 		double pos = solution.decisionVector[i];
// 		//double prob=(rand()/(double)RAND_MAX);
// 		double prob =curand_uniform( &localState);
// 		double delta;
// 		double mutation_prob=(1.0/decisionNumber);
// 		if(prob<mutation_prob){
// 			//double u =(rand()/(double)RAND_MAX);
// 			double u =curand_uniform( &localState);
// 			if(u < 0.5)
// 				delta=pow(2.0*u,1.0/(decisionNumber+1))-1;
// 			else
// 				delta=1-pow(2.0*(1-u),1.0/(decisionNumber+1));
// 		}else
// 			delta=0;
// 		solution.decisionVector[i]=(pos+delta*MAX_MUT);
// 	}
// 	devStates[idx] = localState;
// }

//apply the turbulence operator - JMetal
void Particle::turbulence(){
	double probability=1.0/decisionNumber;
	//int idx = blockIdx.x*blockDim.x + threadIdx.x;
	//curandState localState = devStates[idx];
	
	double rnd, delta1, delta2, mut_pow, deltaq;
	double y, yl, yu, val, xy;
	for (int var=0; var < decisionNumber; var++) {
		//if (curand_uniform( &localState) <= probability)
		if ((rand()/(double)RAND_MAX) <= probability){
			y      = solution.decisionVector[var];
			yl     = inferiorPositionLimit[var];//lower bound            
			yu     = superiorPositionLimit[var];//upper bound
			delta1 = (y-yl)/(yu-yl);
			delta2 = (yu-y)/(yu-yl);
			//rnd = curand_uniform( &localState);
			rnd=(rand()/(double)RAND_MAX);
			mut_pow = 1.0/(20.0+1.0);
			if (rnd <= 0.5){
				xy     = 1.0-delta1;
				val    = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(20.0+1.0)));
				deltaq =  pow(val,mut_pow) - 1.0;
			}
			else{
				xy = 1.0-delta2;
				val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(20.0+1.0)));
				deltaq = 1.0 - (pow(val,mut_pow));
			}
			y = y + deltaq*(yu-yl);
			if (y<yl)
				y = yl;
			if (y>yu)
				y = yu;
			solution.decisionVector[var]= y;
		}
	} // for
	//devStates[idx] = localState;
}

// //keep the position in the swarm specified bounds in IMulti algorithm
// bool Particle::truncatePositionIMulti(double* centroid, double range, int decisionNumber, double* inferiorLimit, double* superiorLimit){ //according to JMetal
// 	bool truncate=false;
// 	
// 	for(int i=0;i<decisionNumber;i++){
// 		double inf = max(inferiorLimit[i], (centroid[i]-range)*(superiorLimit[i]-inferiorLimit[i]));
// 		double sup = min(superiorLimit[i], (centroid[i]+range)*(superiorLimit[i]-inferiorLimit[i]));
// 		if(solution.decisionVector[i]<inf){
// 			solution.decisionVector[i]=inf;
// 			velocity[i] = velocity[i] * -1;
// 			truncate=true;
// 		}
// 		if(solution.decisionVector[i]>sup){
// 			solution.decisionVector[i]=sup;
// 			velocity[i] = velocity[i] * -1;
// 			truncate=true;
// 		}
// 	}
// 	return truncate;
// }
//keep the position in the swarm specified bounds in IMulti algorithm
int Particle::truncatePositionIMulti(double* centroid, double range, double* smallerDecision, double* largerDecision){ //original do André
	int truncate=0;
	
	for(int i=0;i<decisionNumber;i++){
		double var_i=solution.decisionVector[i];
		double inf = std::max(inferiorPositionLimit[i], centroid[i]-(range*(superiorPositionLimit[i]-inferiorPositionLimit[i])) );
		double sup = std::min(superiorPositionLimit[i], centroid[i]+(range*(superiorPositionLimit[i]-inferiorPositionLimit[i])) );
		
// 		printf("%f < %f < %f (%f +- %f)\n", inf, var_i, sup, centroid[i], (range*(superiorPositionLimit[i]-inferiorPositionLimit[i])));
// 		printf("inf: %f inferior: %f centroid: %f range: %f (%f)\n",inf, inferiorPositionLimit[i], centroid[i],range,superiorPositionLimit[i]-inferiorPositionLimit[i]);
		double min=std::max(smallerDecision[i], inf);
		double max=std::min(largerDecision[i], sup);
		
		if(var_i != var_i){//nan detection
			if((rand()/(double)RAND_MAX) < 0.5){//set to outside upper or lower limit randomly
				fprintf(stderr,"WARNING! Invalid solution truncated to lower limit\n");
				var_i=inferiorPositionLimit[i]-1;
			}else{
				fprintf(stderr,"WARNING! Invalid solution truncated to upper limit\n");
				var_i=superiorPositionLimit[i]+1;
			}
		}
		
		if(var_i<inf){
			if(!strcmp(truncType, "imulti")){
				double part1=(sup-inf);
				double part2=fabs(inf - var_i);
				double part3=fmod(part2,part1);
				double part4=sup-part3;
				solution.decisionVector[i]=part4;
			}
			
			if(!strcmp(truncType, "rdmInSwarm"))
				solution.decisionVector[i]=((rand()/(double)RAND_MAX)*(max-min))+min; //new truncation sugested by Roberto
				
			if(!strcmp(truncType, "extremes"))
				solution.decisionVector[i]=inf;
			
			if(!strcmp(truncType, "random"))
				solution.decisionVector[i]=inf+ ((rand()/(double)RAND_MAX)* (sup-inf) );
			
// 			solution.decisionVector[i]=((rand()/(double)RAND_MAX)*(sup-inf))+inf; //new truncation by myself
// 			printf("inf: %f < (%f -> %f) < %f\n", inf, var_i, solution.decisionVector[i], sup);
			
			if(inf != inferiorPositionLimit[i])
				truncate++;
		}
		if(var_i>sup){
			if(!strcmp(truncType, "imulti")){
				double part1=(sup-inf);
				double part2=fabs(var_i-sup);
				double part3=fmod(part2,part1);
				double part4=part3+inf;
				solution.decisionVector[i]=part4;
			}
			
			if(!strcmp(truncType, "rdmInSwarm"))
				solution.decisionVector[i]=((rand()/(double)RAND_MAX)*(max-min))+min; //new truncation sugested by Roberto
				
			if(!strcmp(truncType, "extremes"))
				solution.decisionVector[i]=sup;
  			
//   			solution.decisionVector[i]=((rand()/(double)RAND_MAX)*(sup-inf))+inf; //new truncation by myself
// 			printf("sup: %f < (%f -> %f) < %f\n", inf, var_i, solution.decisionVector[i], sup);
			
			if(!strcmp(truncType, "random"))
				solution.decisionVector[i]=inf+ ((rand()/(double)RAND_MAX)* (sup-inf) );
			
			if(sup != superiorPositionLimit[i])
				truncate++;
		}
// 		if(solution.decisionVector[i] > superiorPositionLimit[i] || superiorPositionLimit[i] < inferiorPositionLimit[i]){
// 			printf("\ntruncation problem -- %f < %f < %f\n", inferiorPositionLimit[i], solution.decisionVector[i], superiorPositionLimit[i]);
// 			exit(1);
// 		}
	}
	return truncate;
}



// void Particle::truncatePositon(){
// 	//JMetal
// 	for (int var = 0; var < decisionNumber; var++) {
// 		if (solution.decisionVector[var] < solution.inferiorPositionLimit) {
// 			solution.decisionVector[var]=solution.inferiorPositionLimit;
// 			velocity[var] = velocity[var] * -1; //
// 		}
// 		if (solution.decisionVector[var] > solution.superiorPositionLimit) {
// 			solution.decisionVector[var]=solution.superiorPositionLimit;
// 			velocity[var] = velocity[var] * -1; //
// 		}
// 	}
// }
