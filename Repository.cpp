class Swarm; //prototype of class Swarm to allow the declaration of the outerSwarm pointer
//class repository
class Repository{
	public:

		Swarm *outerSwarm;// used to obtain parameters from the swarm that contains this repository
		
		double *smallestDecision; //smaller decision values found during the search
		double *largestDecision; //larger decision values found during the search
		
		double *weight; //weight vector associated to this repository, if using weights
		
		char archiverType[50]; //name of the archiver used
		bool archiverUsed; //flag to register if the repository was ever full
		int particlesEntered;//counts solutions that enters in the repository in each generation

		Repository();
		~Repository();
		
		void initialize(const char* At, int Rs);
		//tries to add a solution in the repository if it is non-dominated
		//param - the candidate solution to be inserted in the repository
		bool add(Solution &candidate);
		//tries to add a set of candidate solutions the repository if it is non-dominated
		//param - the candidate solutions to be inserted in the repository
		void add(Repository &rep);
		//after the update the repository gets messy (full of holes), this method make it continuous again.
		void organize();
		//safely get a solution from the repository
		//param - the index of the desired solution
		const Solution& getSolution(int index);
		//safely get the solutions from the repository
		Solution* getSolutions();
		//switch between the implemented archivers
		//param - the solution to be submitted to the archiver
		bool archiver(Solution &candidate);
		//archiver wich remove a random solution of the repository
		//param - the solution to be submitted to the archiver
		bool randomArchiver(Solution &candidate);
		//archiver wich consider the crowding distance to keep/exclude a solution
		//param - the solution to be submitted to the archiver
		bool crowdingDistanceArchiver(Solution &candidate);
		//archiver proposed by Carvalho and Pozo //Removes the solution with worst distance to the ideal point
		//param - the solution to be submitted to the archiver
		bool idealArchiver(Solution &candidate);
		//Removes the solution with highest weighted sum. Needs a weight vector
		//param - the solution to be submitted to the archiver
		bool pbiArchiver(Solution &candidate);
		//Removes the solution with highest weighted sum. Needs a weight vector
		//param - the solution to be submitted to the archiver
		bool tchArchiver(Solution &candidate);
		//Removes the solution with highest weighted sum. Needs a weight vector
		//param - the solution to be submitted to the archiver
		bool WCPArchiver(Solution &candidate);
		//Removes the solution with highest weighted sum. Needs a weight vector
		//param - the solution to be submitted to the archiver
		bool WSumArchiver(Solution &candidate);
		//Removes the solution with highest weighted sum. Needs a weight vector
		//param - the solution to be submitted to the archiver
		bool r_idealArchiver(Solution &candidate);
		 //Removes the solution with highest log-likelihood
		//param - the solution to be submitted to the archiver
		bool largestLikelihoodArchiver(Solution &candidate); //located in the util.cpp file due to cross referencing issues on compilation time
		 //Removes the solution with smallest log-likelihood
		//param - the solution to be submitted to the archiver
		bool smallestLikelihoodArchiver(Solution &candidate); //located in the util.cpp file due to cross referencing issues on compilation time
		 //Removes the solution with highest log-likelihood
		//param - the solution to be submitted to the archiver
		bool largestDistanceArchiver(Solution &candidate); //located in the util.cpp file due to cross referencing issues on compilation time
		//Removes the solution with smallest log-likelihood
		//param - the solution to be submitted to the archiver
		bool smallestDistanceArchiver(Solution &candidate); //located in the util.cpp file due to cross referencing issues on compilation time
		//archiver proposed by Laumanns and Zenklusken
		//param - the solution to be submitted to the archiver
		bool MGAArchiver(Solution &candidate);
		//param - the solution to be submitted to the archiver
		bool R2Archiver(Solution &candidate);
		//param - the solution to be submitted to the archiver
		bool HVArchiver(Solution &candidate);
		//get the actual size of the repository
		int getActualSize(){return actualSize;}
		//clear the repository
		void clear();
		//exclude a solution of the repository using its real index
		//param - the index of the solution to be excluded from the repository
		void exclude(int index);
		//add a solution ignoring whether the repository is full or not
		void forceAdd(Solution &candidate);

	private:
		Repository(const Repository &source){}//copy
		Repository& operator= (const Repository &source){}//assignment
		
		int repositorySize;//maximum particles stored in this repository before the filter is called
// 		Solution solutions[maxRepositorySize+1];
		Solution* solutions;
		int actualSize; //actual number of solutions in the repository
// 		bool controlSolutions[maxRepositorySize+1];//control if one position have or not solutions associated
		bool* controlSolutions;
		//insert a solution in the repository with no verification
		//param - the candidate solution to be inserted in the repository (no verification)
		void insert(const Solution &solution);
};
//constructor initialize the variables
// Repository::Repository(int At, int Rs){
Repository::Repository(){
	repositorySize=-1;
	actualSize=-1;
	solutions=NULL;
	controlSolutions=NULL;
	smallestDecision = new double[decisionNumber]; //smallest decision values found during the search
	largestDecision = new double[decisionNumber]; //largest decision values found during the search
	for(int i=0;i<decisionNumber;i++){//initialize with the worst values possible
		smallestDecision[i]=superiorPositionLimit[i];
		largestDecision[i]=inferiorPositionLimit[i];
	}
	weight = new double[objectiveNumber]; //weight vector associated to this repository, if using weights
}
Repository::~Repository(){
	delete[] smallestDecision;
	delete[] largestDecision;
	delete[] weight;
	if(solutions != NULL){
		delete[] solutions;
		delete[] controlSolutions;
	}
}

//initialize the main variables of the class
void Repository::initialize(const char* At, int Rs){
	if(repositorySize != Rs){
		if(solutions != NULL){ //if there is already memory allocated, re-alocate because the size may have changed
			delete[] solutions;
			delete[] controlSolutions;
		}
		solutions = new Solution[Rs+1];
		controlSolutions = new bool[Rs+1];
	}
	
	memset(controlSolutions, false, sizeof(bool)*(Rs+1));
	
	sprintf(archiverType, "%s",At);
	repositorySize=Rs; //maximum particles stored in this repository before the filter is called
	actualSize=0; //the repository initialize with 0 solutions
	archiverUsed=false;//initialize the flag stating if the filter was ever used
}

//tries to add a solution in the repository if it is non-dominated
//param - the candidate solution to be inserted in the repository
 bool Repository::add(Solution &candidate){
	bool isDominated=false;
	bool equal=false;
	int dom;
	bool enteredArchive=false;
	if(actualSize==0){ //if the repository is empty, insert the solution
		insert(candidate);
		enteredArchive=true;
	}else{
		for(int s=0;s<repositorySize+1;s++){
			if(controlSolutions[s]){//if this solution is valid
				if(!candidate.isEqual(getSolution(s))){ //if the solutions are not equal
					//verify the dominance relation between two vectors
					//return 1 if sol1 dominates sol2, -1 if sol2 dominates sol1, 0 if they do not dominate each other, 2 if they are equal
					if(!strcmp(archiverType, "pbi") || !strcmp(archiverType, "tch")  || !strcmp(archiverType, "wcp") || !strcmp(archiverType, "wsum") || !strcmp(archiverType, "r-ideal"))//repositories based on decomposition, if decomposition do not check dominance
						dom=0;
					else
						dom=dominance(candidate.objectiveVector, getSolution(s).objectiveVector, objectiveNumber);

					if(dom == 1){//if the candidate dominates the solution in the repository
						exclude(s);
					}else{
						if(dom == -1){//if a solution in the repository dominates the candidate
							isDominated=true;
							if(vectorZero(getSolution(s).objectiveVector, objectiveNumber)){
								fprintf(stderr, "\nERROR! Trying to insert in the repository a solution whose objectives are all 0\n");
								exit(1);
							}
							break;
						}
					}
				}else{ //if the solutions are equal, discard the candidate
					equal=true;
					break;
				}
			}
		}

		if(!isDominated && !equal){ //if the solution is non-dominated
			candidate.dominated=false;
			if(actualSize<repositorySize){//if the repository is not empty nor full
				insert(candidate);//insert the solution
				enteredArchive=true;
			}else{ //if the repository is full
				enteredArchive=archiver(candidate);
			}
		}else{
			candidate.dominated=true;
			enteredArchive=false;
		}
	}
	return enteredArchive;
}
//tries to add a set of candidate solutions the repository if it is non-dominated
//param - the candidate solutions to be inserted in the repository
void Repository::add(Repository &rep){
	if(rep.getActualSize()+actualSize > repositorySize){
		fprintf(stderr,"\nWARNING repository already contains %d solutions, you are trying to add more %d. The total size is %d, the maximum size is %d\n", actualSize, rep.getActualSize(), actualSize+rep.getActualSize(), repositorySize);
	}
	
	for(int s=0;s<rep.getActualSize();s++){
		add(rep.getSolutions()[s]);
	}
}
//insert a solution in the repository with no verification
//param - the candidate solution to be inserted in the repository (no verification)
void Repository::insert(const Solution &solution){
	for(int s=0;s<repositorySize+1;s++){
		if(!controlSolutions[s]){//if the solution is invalid
			//solutions[s].deepCopy(solution);
// 			memcpy(&solutions[s], &solution, sizeof(Solution));
// 			for(int o=0;o<objectiveNumber;o++)
// 				solutions[s].objectiveVector[o]=solution.objectiveVector[o];
// 			for(int d=0;d<decisionNumber;d++)
// 				solutions[s].decisionVector[d]=solution.decisionVector[d];
			solutions[s]=solution;
			
			controlSolutions[s]=true; //the solution is now valid
			actualSize++;
			particlesEntered++;
			//break;
			return;
		}
	}
	fprintf(stderr,"\nERROR INSERTING SOLUTION!!!");
	exit(1);
}
//exclude a solution of the repository using its real index
//param - the index of the solution to be excluded from the repository
void Repository::exclude(int index){
	if(controlSolutions[index]){//if the solution is valid
		controlSolutions[index]=false; //the solution can be overriten
		actualSize--;
	}else{
		fprintf(stderr,"\nERROR! TRYING TO REMOVE AN INVALID SOLUTION (%d).\n", index);
		exit(1);
	}
}
//after the update the repository gets messy (full of holes), this method make it continuous again.
void Repository::organize(){
	for(int i=0;i<actualSize;i++){ // verify the valid solutions
		if(!controlSolutions[i]){//if the solution is invalid
			for(int s=actualSize;s<repositorySize+1;s++){ // verify the solutions in the positions they shouldnt be
				if(controlSolutions[s]){//if the solution is valid
// 					solutions[i]=solutions[s];
// 					memcpy(&solutions[i],&solutions[s],sizeof(Solution));
// 					for(int o=0;o<objectiveNumber;o++)
// 						solutions[i].objectiveVector[o]=solutions[s].objectiveVector[o];
// 					for(int d=0;d<decisionNumber;d++)
// 						solutions[i].decisionVector[d]=solutions[s].decisionVector[d];
					solutions[i]=solutions[s];
					
					controlSolutions[i]=true;
					controlSolutions[s]=false;
					break;
				}
			}
		}
	}
}
//safely get a solution from the repository
//param - the index of the desired solution
const Solution& Repository::getSolution(int index){
	if(!controlSolutions[index]){
		fprintf(stderr,"\nERROR, TRYING TO GET AN INVALID SOLUTION FROM THE REPOSITORY (sol:%d repSize: %d)\n", index, actualSize);
		exit(1);
	}
// 	Solution sol = solutions[index];
// 	return sol;
	return solutions[index];
}

//safely the solutions from the repository
Solution* Repository::getSolutions(){
	Solution* sol= solutions;
	return sol;
}
//switch between the implemented archivers
//param - the solution to be submitted to the archiver
 bool Repository::archiver(Solution &candidate){
	archiverUsed=true;
	if(!strcmp(archiverType, "rnd"))
		return randomArchiver(candidate);
	if(!strcmp(archiverType, "cd"))
		return crowdingDistanceArchiver(candidate);
	if(!strcmp(archiverType, "ideal"))
		return idealArchiver(candidate);
	if(!strcmp(archiverType, "pbi"))
		return pbiArchiver(candidate);
	if(!strcmp(archiverType, "tch"))
		return tchArchiver(candidate);
	if(!strcmp(archiverType, "wcp"))
		return WCPArchiver(candidate);
	if(!strcmp(archiverType, "wsum"))
		return WSumArchiver(candidate);
	if(!strcmp(archiverType, "r-ideal"))
		return r_idealArchiver(candidate);
	if(!strcmp(archiverType, "larger-likelihood"))
		return largestLikelihoodArchiver(candidate);
	if(!strcmp(archiverType, "smaller-likelihood"))
		return smallestLikelihoodArchiver(candidate);
	if(!strcmp(archiverType, "larger-distance"))
		return largestDistanceArchiver(candidate);
	if(!strcmp(archiverType, "smaller-distance"))
		return smallestDistanceArchiver(candidate);
	if(!strcmp(archiverType, "mga"))
		return MGAArchiver(candidate);
	if(!strcmp(archiverType, "r2"))
		return R2Archiver(candidate);
	if(!strcmp(archiverType, "hv"))
		return HVArchiver(candidate);
	
	if(strcmp(archiverType, "rnd") && strcmp(archiverType, "cd") && strcmp(archiverType, "ideal") && strcmp(archiverType, "mga") && strcmp(archiverType, "pbi") && strcmp(archiverType, "tch") && strcmp(archiverType, "wcp") && strcmp(archiverType, "wsum") && strcmp(archiverType, "r-ideal") && strcmp(archiverType, "larger-likelihood") && strcmp(archiverType, "smaller-likelihood") && strcmp(archiverType, "larger-distance") && strcmp(archiverType, "smaller-distance") && strcmp(archiverType, "r2") && strcmp(archiverType, "hv") ){
		fprintf(stderr, "INVALID ARCHIVER! (%s)\n", archiverType);
		exit(1);
	}
}
//archiver wich remove a random solution of the repository
//param - the solution to be submitted to the archiver
bool Repository::randomArchiver(Solution &candidate){
	while(actualSize >= repositorySize)
		exclude(rand() % repositorySize);
	insert(candidate);
	return true; //always insert the new solution
}
//archiver proposed by Carvalho and Pozo //Removes the solutio with worst distance to the ideal point
//param - the solution to be submitted to the archiver
bool Repository::idealArchiver(Solution &candidate){
	insert(candidate);
	organize();
	
	double solNorm[objectiveNumber];
	double idealNorm[objectiveNumber];
	double higherDistance=-1;
	double distance=0;
	int index=-1;
	
	memset(idealNorm, 0, sizeof(double)*objectiveNumber); //since the points are normalized between 0 and 1, the ideal point will be (0,...,0)
		
	for(int i=0;i<actualSize;i++){
		normalizeObjectives(solNorm, solutions[i].objectiveVector);
		
		distance=calculateEuclideanDistance(solNorm, idealNorm, objectiveNumber);
		
		if(distance > higherDistance){
			higherDistance=distance;
			index=i;
		}
	}
	
	if(index==-1){
		fprintf(stderr,"\nIDEAL ARCHIVER ERROR %f\n", higherDistance);
		exit(1);
	}
	
	bool ret=true;
	if(candidate.isEqual(solutions[index]))
		ret=false;
	
	exclude(index);
	
	return ret;
	
// 	/*//For each solution on the front, it calculates the distance to the ideal point
// 	double smallerDistance[actualSize];
// 	
// 	for(int i=0;i<actualSize;i++){
// 		double norm[objectiveNumber];
// 		double sp_norm[objectiveNumber];
// 		for(int j=0;j<objectiveNumber;j++){
// 			norm[j]=normalize(solutions[i].objectiveVector[j],smallerObjs[j],largerObjs[j]);
// 			sp_norm[j]=normalize(smallerObjs[j],smallerObjs[j],largerObjs[j]);
// 		}
// 		smallerDistance[i] = calculateEuclideanDistance(sp_norm, norm, objectiveNumber);
// 		//smallerDistance[i] = calculateEuclideanDistance(smallerObjs, solutions[i].objectiveVector, objectiveNumber);
// 	}
// 
// 	double highDistanceValue = -1.0;
// 	int index = -1;
// 	for (int i = 0; i<actualSize; i++) {
// 		if(smallerDistance[i] > highDistanceValue){
// 			highDistanceValue =smallerDistance[i];
// 			index = i;
// 		}
// 	}*/
}

bool Repository::pbiArchiver(Solution &candidate){
	//http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7020200
	//double weight[]={0.1,0.9};
	insert(candidate);
	organize();
	
	double highDistanceValue=MAXDOUBLE*-1;
	int higherIndex=-1;
	double solNorm[objectiveNumber];
// 	double *solNorm;
	
	for(int i=0;i<actualSize;i++){
		double norm;
		double pbi=0;
		normalizeObjectives(solNorm, solutions[i].objectiveVector);
// 		solNorm=solutions[i].objectiveVector;
		
		pbi=PBI(solNorm, weight);
		if(pbi > highDistanceValue){
			highDistanceValue=pbi;
			higherIndex=i;
		}
	}

	if(higherIndex==-1){
		fprintf(stderr,"\nPBI ARCHIVER ERROR %f %f\n", PBI(solNorm, weight), highDistanceValue);
		candidate.print();
		printf("\nweight: ");
		printVector(weight, objectiveNumber);
		printf("\n");
		exit(1);
	}
	
	bool ret=true;
	if(candidate.isEqual(solutions[higherIndex]))
		ret=false;
	
	exclude(higherIndex);
	
	return ret;
}

bool Repository::tchArchiver(Solution &candidate){
	//http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7020200
	insert(candidate);
	organize();
	
	double highDistanceValue=MAXDOUBLE*-1;
	int higherIndex=-1;
	double solNorm[objectiveNumber];
// 	double *solNorm;
	
	for(int i=0;i<actualSize;i++){
		double norm;
		double tch=0;
		normalizeObjectives(solNorm, solutions[i].objectiveVector);
// 		solNorm=solutions[i].objectiveVector;
		
		tch=TCH(solNorm, weight);
		if(tch > highDistanceValue){
			highDistanceValue=tch;
			higherIndex=i;
		}
	}
	
	if(higherIndex==-1){
		fprintf(stderr,"\nTCH ARCHIVER ERROR %f\n", highDistanceValue);
		printVector(solNorm, objectiveNumber);
		fprintf(stderr,"\n");
		exit(1);
	}
	
	bool ret=true;
	if(candidate.isEqual(solutions[higherIndex]))
		ret=false;
	
	exclude(higherIndex);
	
	return ret;
}
bool Repository::WCPArchiver(Solution &candidate){
	//http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7020200
	//double weight[]={0.1,0.9};
	insert(candidate);
	organize();
	
	double highDistanceValue=MAXDOUBLE*-1;
	int higherIndex=-1;
	
	for(int i=0;i<actualSize;i++){
		double norm[objectiveNumber];
		double wsum=0;
		normalizeObjectives(norm, solutions[i].objectiveVector);
		
		for(int j=0;j<objectiveNumber;j++){
			wsum+=(weight[j]*norm[j])*(weight[j]*norm[j])*(weight[j]*norm[j]); //x^3
		}
		// 		wsum=PBI(solutions[i].objectiveVector, weight, smallerObjs, largerObjs);
		if(wsum > highDistanceValue){
			highDistanceValue=wsum;
			higherIndex=i;
		}
	}
	if(higherIndex==-1){
		fprintf(stderr,"\nWCP ARCHIVER ERROR %f\n", highDistanceValue);
		exit(1);
	}
	bool ret=true;
	if(candidate.isEqual(solutions[higherIndex]))
		ret=false;
	
	exclude(higherIndex);
	
	return ret;
}

bool Repository::WSumArchiver(Solution &candidate){
	//http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7020200
	//double weight[]={0.1,0.9};
	insert(candidate);
	organize();
	
	double highDistanceValue=MAXDOUBLE*-1;
	int higherIndex=-1;
	
	for(int i=0;i<actualSize;i++){
		double norm[objectiveNumber];
		double wsum=0;
		normalizeObjectives(norm, solutions[i].objectiveVector);
		
		for(int j=0;j<objectiveNumber;j++){
			wsum+=(weight[j]*norm[j]);
		}
		if(wsum > highDistanceValue){
			highDistanceValue=wsum;
			higherIndex=i;
		}
	}
	if(higherIndex==-1){
		fprintf(stderr,"\nWSum ARCHIVER ERROR %f\n", highDistanceValue);
		exit(1);
	}
	bool ret=true;
	if(candidate.isEqual(solutions[higherIndex]))
		ret=false;
	
	exclude(higherIndex);
	
	return ret;
}
bool Repository::r_idealArchiver(Solution &candidate){
	insert(candidate);
	organize();
	
	double highDistanceValue=MAXDOUBLE*-1;
	int higherIndex=-1;
	
	for(int i=0;i<actualSize;i++){
		double norm[objectiveNumber];
		double dist=0;
		normalizeObjectives(norm, solutions[i].objectiveVector);

		dist=calculateEuclideanDistance(norm, weight, objectiveNumber);
		if(dist > highDistanceValue){
			highDistanceValue=dist;
			higherIndex=i;
		}
	}
	if(higherIndex==-1){
		fprintf(stderr,"\nR-IDEAL ARCHIVER ERROR %f\n", highDistanceValue);
		exit(1);
	}
	bool ret=true;
	if(candidate.isEqual(solutions[higherIndex]))
		ret=false;
	
	exclude(higherIndex);
	
	return ret;
}
// void Repository::largestLikelihoodArchiver(Solution candidate){} ////located in the util.cpp file due to cross referencing issues on compilation time
// void Repository::smallestLikelihoodArchiver(Solution candidate){} ////located in the util.cpp file due to cross referencing issues on compilation time
// void Repository::largestDistanceArchiver(Solution candidate){} ////located in the util.cpp file due to cross referencing issues on compilation time
// void Repository::smallestDistanceArchiver(Solution candidate){} ////located in the util.cpp file due to cross referencing issues on compilation time

//archiver proposed by Laumanns and Zenklusken
//param - the solution to be submitted to the archiver
bool Repository::MGAArchiver(Solution &candidate){
	insert(candidate);
	organize();
	int b = compute_b_mga(solutions, actualSize);
	int index_removed = -1;
	while(index_removed==-1){
		for(int i = actualSize-1; i>=0;i--){
			double box_i[objectiveNumber];
			box_mga(solutions[i], box_i, b, objectiveNumber);
			for(int j = actualSize-1; j>=0;j--){
				if(i!=j){
					double box_j[objectiveNumber];
					box_mga(solutions[j], box_j, b, objectiveNumber);
					int comparation = dominance(box_i, box_j, objectiveNumber);
					if(comparation == -1 || isEqual(box_i, box_j, objectiveNumber) ){
						index_removed = i;
						break;
					}
				}
			}
			//if(index_removed!=-1){
			//	break;
			//}
		}
		b--;
	}
	
	bool ret=true;
	if(candidate.isEqual(solutions[index_removed]))
		ret=false;
	
	exclude(index_removed);
	
	return ret;
}

//archiver wich consider the contributing r2 to keep/exclude a solution
//param - the solution to be submitted to the archiver
bool Repository::R2Archiver(Solution &candidate){
	double smallestR2=MAXDOUBLE;
	int idxSmallestR2=-1;
	
	insert(candidate);
	
	for(int i=0;i<actualSize;i++)
		getSolutions()[i].crowdingDistance=-1;
	
	organize();
	updateContributingR2(*this); //update the contributing R2 of all solutions
// 	updateContributingHypervolume(*this);

	
	for(int i=0;i<actualSize;i++){//find the valid solution with smallest contributing r2 (stored on the crowding distance field)
		if(controlSolutions[i] && getSolution(i).crowdingDistance<=smallestR2){
			smallestR2=getSolution(i).crowdingDistance;
			idxSmallestR2=i;
		}
	}
	
	if(idxSmallestR2 == -1){
		fprintf(stderr,"\nR2 archiver error!\n The R2 of all particles in the repository is:\n");
		for(int i=0;i<actualSize;i++)
			fprintf(stderr,"%f\n", getSolution(i).crowdingDistance);
		exit(1);
	}
	
// 	fprintf(stderr,"%d --> %f  (%d)\n",idxSmallerR2, smallerR2, actualSize);
// 	if(actualSize > repositorySize)//test of removing all solutions with no contribution
	
	bool ret=true;
	if(candidate.isEqual(solutions[idxSmallestR2]))
		ret=false;
	
	exclude(idxSmallestR2); //remove the solution with the smallest contribution
	
	return ret;
}

//archiver wich consider the contributing HV to keep/exclude a solution
//param - the solution to be submitted to the archiver
bool Repository::HVArchiver(Solution &candidate){
	double smallestHV=MAXDOUBLE;
	int idxSmallestHV=-1;
	
	insert(candidate);
	
	for(int i=0;i<actualSize;i++)
		getSolutions()[i].crowdingDistance=-1;
	
	organize();
	updateContributingHypervolume(*this); //update the contributing HV of all solutions
	// 	updateContributingHypervolume(*this);
	
	
	for(int i=0;i<actualSize;i++){//find the valid solution with smallest contributing HV (stored on the crowding distance field)
		if(controlSolutions[i] && getSolution(i).crowdingDistance<=smallestHV){
			smallestHV=getSolution(i).crowdingDistance;
			idxSmallestHV=i;
		}
	}
	
	if(idxSmallestHV == -1){
		fprintf(stderr,"\nHV archiver error!\n The HV of all particles in the repository is:\n");
		for(int i=0;i<actualSize;i++)
			fprintf(stderr,"%f\n", getSolution(i).crowdingDistance);
		exit(1);
	}
	
	// 	fprintf(stderr,"%d --> %f  (%d)\n",idxSmallerR2, smallerR2, actualSize);
	// 	if(actualSize > repositorySize)//test of removing all solutions with no contribution
	
	bool ret=true;
	if(candidate.isEqual(solutions[idxSmallestHV]))
		ret=false;
	
	exclude(idxSmallestHV); //remove the solution with the smallest contribution
	
	return ret;
}

//archiver wich consider the crowding distance to keep/exclude a solution
//param - the solution to be submitted to the archiver
bool Repository::crowdingDistanceArchiver(Solution &candidate){
	//Solution temp[repositorySize+1]; //create a temporary repository to contain all the solutions plus the candidate
	double smallerCrowdingDistance=MAXDOUBLE;
	int idxSmallerCrowdingDistance=-1;

	//sync(temp); //sincronize the new repository with the solutions already found
//	for(int i=0;i<actualSize;i++)
//		temp[i]=solutions[i];
	//memcpy(temp, solutions, sizeof(Solution)*actualSize);

	//solutions[actualSize]=candidate;//insert the new one
	insert(candidate);
		
	organize();
	updateCrowdingDistances(solutions, actualSize); //update the crowing distances

	for(int i=0;i<actualSize;i++){//find the valid solution with smallest crowding distance
		if(controlSolutions[i] && getSolution(i).crowdingDistance<=smallerCrowdingDistance){
			smallerCrowdingDistance=getSolution(i).crowdingDistance;
			idxSmallerCrowdingDistance=i;
		}
	}
	
	if(idxSmallerCrowdingDistance == -1){
		fprintf(stderr,"\nCrowding Distance archiver error!\n The crowding distances of all particles in the repository is:\n");
		for(int i=0;i<actualSize;i++)
			fprintf(stderr,"%f\n", getSolution(i).crowdingDistance);
		exit(1);
	}
	
	
//	if(!solutions[idxSmallerCrowdingDistance].isEqual(candidate)){ //if the candidate is not the solution with smallest crowding distance
	bool ret=true;
	if(candidate.isEqual(solutions[idxSmallerCrowdingDistance]))
		ret=false;
	
	exclude(idxSmallerCrowdingDistance); //remove the solution with the smallest crowding distance
	
	return ret;
	//insert(candidate);//insert the new solution
//	}

	//free(temp);
}

//clear the repository
void Repository::clear(){
	actualSize=0; //the repository initialize with 0 solutions
	archiverUsed=false;
	memset(controlSolutions, false, sizeof(bool)*(repositorySize+1));
}

//add a solution ignoring whether the repository is full or not
void Repository::forceAdd(Solution &candidate){
	bool isDominated=false;
	bool equal=false;
	if(actualSize==0){ //if the repository is empty, insert the solution
		insert(candidate);
	}else{
		for(int s=0;s<repositorySize+1;s++){
			if(controlSolutions[s]){//if this solution is valid
				if(!candidate.isEqual(getSolutions()[s])){ //if the solutions are not equal
					//verify the dominance relation between two vectors
					//return 1 if sol1 dominates sol2, -1 if sol2 dominates sol1, 0 if they do not dominate each other, 2 if they are equal
					if(dominance(candidate.objectiveVector, getSolution(s).objectiveVector, objectiveNumber) == 1){//if the candidate dominates the solution in the repository
						exclude(s);
					}else{
						if(dominance(candidate.objectiveVector, getSolution(s).objectiveVector, objectiveNumber) == -1){//if a solution in the repository dominates the candidate
							isDominated=true;
							break;
						}
					}
				}else{ //if the solutions are equal, discard the candidate
					equal=true;
					break;
				}
			}
		}

		if(!isDominated && !equal){ //if the solution is non-dominated
			if(actualSize<repositorySize+1){//if there is memory left
				insert(candidate);//insert the solution
			}else{ //if there is not memory left
				fprintf(stderr,"REPOSITORY MEMORY UNAVAILABLE, INCREASE THE REPOSITORY MAXIMUM SIZE\n");
				//exit(1);
			}
		}
	}
}