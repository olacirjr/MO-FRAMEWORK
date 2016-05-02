struct Neighbor{
	int index;
	double distance;
};
void r_BOA(Swarm &sw, int swarm_edaSize, bool init);
void cma_es(Swarm &sw);
class Swarm{
	public:
		double modelQuality;//stores the quality of the model
		
		Particle* particles;
		
		Neighbor* neighborhood; //neighborhood for MOEA/D like approach
		int neighborhoodSize;
		
		Repository repository; //repository of the swarm
		double *centroid; //centroid of the swarm in I-Multi like approaches
		
		Swarm();
		~Swarm();
		int getSize();
		void setSize(int size);
		
		//initialize the particles of the population
		void initializeParticles();
		//special method to initialize the particles of the population setting as initial local leader a ramdom solution from a given set
		void initializeParticles(Solution* candidateSet, int candidateNumber);
		//initialize the repository with the first non-dominated solutions
		void initializeRepository();
		//choose the global leaders according to the predefined strategy (PSO)
		void chooseGlobalLeaders(Repository &rep);
		//calculate the velocity of the particles (PSO)
		void calculateVelocity();
		//update the position of the particles based on its velocity, leaders and current positions (PSO)
		void updatePosition();
		//apply the turbulence factor in a predefined percentage of the particles (PSO)
		void turbulence();
		//evaluate the particles according to its current position (PSO)
		void evaluation();
		//update the repository with the best new non-dominated solutions
		void updateRepository();
		//update the local leader of the particles
		void updateParticlesMemory();
		
		void PSO();
		
		void rBOA();
		
		void CMAES();
		
		//****************  CMA-ES parameters  ********//
		double sigma, det;
		double **C;
		double **B;
		double *pC, *pSigma, *mean, *D;
		bool init;
		int gen;
		int solsKept;//counter of solutions that are not injected from other swarms, i.e. kept from the previous iteration
		//********************************************//
	private:
		int swarmSize; //population
		Swarm(const Swarm &source){}//copy
		Swarm& operator= (const Swarm &source){}//assignment
};

Swarm::Swarm(){
	swarmSize=0;
	centroid = new double[decisionNumber+objectiveNumber];
	repository.outerSwarm=this;
	particles=NULL;
	neighborhood=NULL;
	
	if(!strcmp(algorithm, "cmaes-mopso")){
		init=true;
		gen=0;
		pC = new double[decisionNumber];
		pSigma = new double[decisionNumber];
		mean = new double[decisionNumber];
		D = new double[decisionNumber];
		C = new double*[decisionNumber];
		B = new double*[decisionNumber];
		for(int i=0; i<decisionNumber;i++){
			mean[i]=(rand()/(double)RAND_MAX);//start the mean in [0,1] since the decision vectors are normalized for CMA-ES
			B[i] = new double[decisionNumber];
			C[i] = new double[decisionNumber];
		}
	}
}

Swarm::~Swarm(){
	delete[] centroid;
	repository.outerSwarm=NULL;
	if(particles != NULL){
		delete[] particles;
		particles=NULL;
	}
	if(neighborhood != NULL){
		delete[] neighborhood;
		neighborhood=NULL;
	}
	
	if(!strcmp(algorithm, "cmaes-mopso")){
		delete[] pC;
		delete[] pSigma;
		delete[] mean;
		delete[] D;
		for(int i=0; i<decisionNumber;i++){
			delete[] B[i];
			delete[] C[i];
		}
		delete[] B;
		delete[] C;
	}
}

int Swarm::getSize(){
	return swarmSize;
}

void Swarm::setSize(int size){
	if(swarmSize == size)//if both are the same size, skip
		return;

	if(particles != NULL)
		delete[] particles;

	swarmSize=size;
	particles = new Particle[swarmSize];
}

void Swarm::initializeParticles(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	funcEvals=0;

	for(int p=0;p<swarmSize;p++){
		sprintf(particles[p].leaderType, "%s", leader);
		particles[p].initialize();
		particles[p].weight=repository.weight;//particle has a pointer pointing to the weight vector of its repository
	}
	if(repository.getActualSize() > 0)
		repository.clear();
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Particles initialized in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::initializeParticles(Solution* candidateSet, int candidateNumber){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	for(int p=0;p<swarmSize;p++){
		sprintf(particles[p].leaderType, "%s", leader);
		particles[p].initialize(candidateSet, candidateNumber);
		particles[p].weight=repository.weight;//particle has a pointer pointing to the weight vector of its repository
	}
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Particles initialized in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::initializeRepository(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	repository.initialize(archiver, repositorySize);
// 	for(int p=0;p<swarmSize;p++){
// 		repository.add(particles[p].solution); //tries to insert the solutions in the repository
// 	}
// 	repository.organize(); //finishes with an organized reposiory
	updateRepository();
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Repository initialized in ms =\t%03.2f\n\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::chooseGlobalLeaders(Repository &rep){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	if(rep.getActualSize()==0){
		fprintf(stderr,"\nERROR ON CHOOSING GLOBAL LEADER! EMPTY REPOSITORY\n");
		exit(1);
	}
	
	if(!strcmp(particles[0].leaderType, "cd") || !strcmp(algorithm, "hmopso")) //if the selection method is crowding distance or hmopso, which uses several leader methods, then update the crowding distances
		updateCrowdingDistances(rep.getSolutions(), rep.getActualSize());
	
	for(int p=0;p<swarmSize;p++)
		particles[p].chooseGlobalLeader(rep.getSolutions(), rep.getActualSize());
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Global leaders chosen in in ms=\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::calculateVelocity(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	for(int p=0;p<swarmSize;p++)
		particles[p].computeSpeed();
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Velocities computed in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::updatePosition(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	for(int p=0;p<swarmSize;p++){
		particles[p].updatePosition();
	}
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Positions updated in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::turbulence(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	for(int p=0;p<swarmSize;p++){
		if(p % 6 == 0){
			particles[p].turbulence();
		}
	}
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Turbulence applied in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::evaluation(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
// 	for(int o=0;o<objectiveNumber;o++){//reset the largest and smallest values before evaluation
// 		maxObjectives[o]=MAXDOUBLE*-1;
// 		minObjectives[o]=MAXDOUBLE;
// 	}
// 	updateMaxMinObjs();
	if(!strcmp(truncType, "rdmInSwarm"))
		updateLargestSmallestDecVectors(*this);
	for(int p=0;p<swarmSize;p++){
		if(!strcmp(algorithm, "imulti") || !strcmp(algorithm, "cmulti") || !strcmp(algorithm, "cmaes-mopso")){
			//tmp_trunc[s]+=(double)swarm[s].particles[p].truncatePositionIMulti(swarm[s].centroid, range, swarm[s].repository.smallerDecision, swarm[s].repository.largerDecision); //to count how many solutions were truncated
			particles[p].truncatePositionIMulti(centroid, range, repository.smallestDecision, repository.largestDecision);
		}else
			particles[p].truncatePositionIMulti(particles[p].solution.decisionVector, range, repository.smallestDecision, repository.largestDecision);
		particles[p].solution.evaluate();
		
		// 		if(!diversityPhase){
		// 			//truncDims/=swarmSize;
		// 			tmp_trunc[s]/=swarmSize;
		// 			//printf("%f ",truncDims);
		// 		}
	}
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Solutions evaluated in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::updateParticlesMemory(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	for(int p=0;p<swarmSize;p++)
		particles[p].updateLocalLeader();
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Local leaders updated in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::updateRepository(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	//use local optimizer for tandem problems on decomposition
	if(!strncmp(problem, "tandem", 6) && decomposition){
		for(int p=0;p<swarmSize;p++)
			particles[p].solution.localOptimizer(repository.weight);
	}
	
	repository.particlesEntered=0; //initialize the counters of numbers of particles that enters in the repository per iteration
	for(int p=0;p<swarmSize;p++){
		int updatedSolutions=0;
		if(repository.add(particles[p].solution)){ //tries to insert the solutions in the repository
			updatedSolutions++;
		}
		if(decomposition && swarmNumber > 1){
			repGlobal->add(particles[p].solution);//tries to add the particle to the global repository
			
			int solsUsed=-1;
			if( (rand()/(double)RAND_MAX) < delta )//uses only the neighborhood
				solsUsed=neighborhoodSize;
			else
				solsUsed=swarmNumber;
			
			//shuffles the neighborhood so when updating a limited number of solutions, there is no bias
			int neighborhoodShuffled[solsUsed];
			for(int i=0;i<solsUsed;i++)
				neighborhoodShuffled[i]=neighborhood[i].index;//shuffles the neighborhood
			std::random_shuffle( neighborhoodShuffled+1, neighborhoodShuffled+solsUsed ); //do not include the first neighbor (itself) on the shuffling
			
			for(int i=1;i<solsUsed;i++){ //tries to insert the particle in all swarms, except the first, which is this one and is already inserted
// 				int neighbor=neighborhood[i].index;
				int neighbor=neighborhoodShuffled[i];
				if(updatedSolutions < maxReplacements){
					if(swarms[neighbor].repository.add(particles[p].solution)){
						updatedSolutions++;
// 						printf("swarm %d updated neighbor %d -- %d\n", this->neighborhood[0].index, neighbor, updatedSolutions);
						swarms[neighbor].repository.organize();//organize the repository if the particle entered
					}
				}
			}
		}
	}
	repository.organize();
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Repository updated in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::PSO(){
	
	updateParticlesMemory();
	
	if(decomposition){
		Repository repTemp;
		mergeNeighboringRepositories(repTemp, *this);
		chooseGlobalLeaders(repTemp);
	}else
		chooseGlobalLeaders(repository);
	
	calculateVelocity();
	
	updatePosition();
	
	turbulence();
}

void Swarm::rBOA(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	for(int p=0;p<swarmSize;p++){
		if(vectorNan(particles[p].solution.decisionVector, decisionNumber)){
			fprintf(stderr,"\n ERROR! NaN on decision vector before EDA application\n");
			exit(1);
		}
		if(vectorNan(particles[p].solution.objectiveVector, objectiveNumber)){
			fprintf(stderr,"\n ERROR! NaN on objective vector before EDA application\n");
			exit(1);
		}
	}
	
	if(repository.getActualSize() > 2 || decomposition){ //if there are enough particles, the solutions are updated normally (not overriten)
		r_BOA(*this, swarmSize, init);
		if(init) init=false;
	}
	
	for(int p=0;p<swarmSize;p++){
		if(vectorNan(particles[p].solution.objectiveVector, objectiveNumber)){
			fprintf(stderr,"\n ERROR! NaN on objective vector after EDA application\n");
			exit(1);
		}
		if(vectorNan(particles[p].solution.decisionVector, decisionNumber)){
			fprintf(stderr,"\n ERROR! NaN on decision vector after EDA application\n");
			exit(1);
		}
	}
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "rBOA processed in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}
	
void Swarm::CMAES(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	for(int p=0;p<swarmSize;p++){
		if(vectorNan(particles[p].solution.decisionVector, decisionNumber)){
			fprintf(stderr,"\n ERROR! NaN on decision vector before EDA application\n");
			exit(1);
		}
		if(vectorNan(particles[p].solution.objectiveVector, objectiveNumber)){
			fprintf(stderr,"\n ERROR! NaN on objective vector before EDA application\n");
			exit(1);
		}
	}
	
// 	if(iteration >0){
// 		for(int s=0;s<swarmNumber;s++)
// 			likelihood_log[s]=logLikelihood(swarms[s].repository.getSolution(0).decisionVector, swarms[sw].mean, swarms[sw].C, swarms[sw].B, swarms[sw].D, swarms[sw].det);
// 		printVectorToFile(likelihood_log, swarmNumber, _matCovB);
// 	}
// 
// 	if(sw==0){//print the cov matrix
// 		printf("\n\nCov mat bef on it: %d \n", iteration);
// 		for (int i = 0; i < decisionNumber; ++i){
// 			printVector(swarms[0].C[i], decisionNumber);
// 			printVectorToFile(swarms[0].C[i], decisionNumber, _matCovB);
// 			printf("\n");
// 		}
// 	}
	
	cma_es(*this);
	
// 	if(sw==50){//print the cov matrix
// // 		printf("\n\nCov mat aft on it: %d \n", iteration);
// 		for (int i = 0; i < decisionNumber; ++i){
// 			for(int j=0;j<decisionNumber;j++){
// 				if(i>j )
// 					printf("%f ",swarms[sw].C[i][j]);
// 				else
// 					printf("%f ",swarms[sw].C[j][i]);
// 			}
// 			printf("\n");
// 		}
// 		printf("\n\n\n");
// 	}
// 
// 	if(iteration >0){
// 		for(int s=0;s<swarmNumber;s++)
// 			likelihood_log[s]=logLikelihood(swarms[s].repository.getSolution(0).decisionVector, swarms[sw].mean, swarms[sw].C, swarms[sw].B, swarms[sw].D, swarms[sw].det);
// 		printVectorToFile(likelihood_log, swarmNumber, _matCovA);
// 	}
	
	for(int p=0;p<swarmSize;p++){
		if(vectorNan(particles[p].solution.objectiveVector, objectiveNumber)){
			fprintf(stderr,"\n ERROR! NaN on objective vector after EDA application\n");
			exit(1);
		}
		if(vectorNan(particles[p].solution.decisionVector, decisionNumber)){
			fprintf(stderr,"\n ERROR! NaN on decision vector after EDA application\n");
			
// 		if(!strcmp(algorithm, "cmaes-mopso")){
// 			for(int s=0;s<swarmNumber;s++){
// 				printf("\n------------------Final covariance matrix-------------------");
// 				for(int i=0;i<decisionNumber;i++){
// 					printf("\n");
// 					printVector(swarms[s].C[i], decisionNumber);
// 				}
// 				printf("\n-------------------------------------\n");
// 			}
// 		}
			
			exit(1);
		}
	}
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "CMA-ES processed in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}