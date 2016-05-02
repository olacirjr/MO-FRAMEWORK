void myLearn(Swarm &sw);
bool mySample (Swarm &sw);
double cSigma, dSigma;

// #include "cmaes.c"
// cmaes_t evo; /* an CMA-ES type struct or "object" */
// double maxDifference=0;

double cmaes_random_Gauss();
static void QLalgo2 (int n, double *d, double *e, double **V);
static void Householder2(int n, double **V, double *d, double *e);
static void Eigen( int N,  double **C, double *diag, double **Q, double *rgtmp);
static double rgdouMax( const double *rgd, int len);
static double rgdouMin( const double *rgd, int len);

static bool logCMAES=false;

void myLearn(Swarm &sw){
	int N=decisionNumber, hsig;
	double oldMean[N], initialStds[N], BDz[N];
	double muEff, cc, muCov, cCov, trace=0.0, chiN;
	cSigma=0; dSigma=0;
	//double sigma, C[N][N], pC[N], pSigma[N], D[N], mean[N], B[N][N];
	
	//******************  treating the input data ********************//
	Repository repTemp;
	
	//gathering the solutions to be used to learn and putting them on repTemp
	if(decomposition){
		mergeNeighboringRepositories(repTemp, sw);
	}else{
		repTemp.initialize(archSubSwarms, sw.getSize()+sw.repository.getActualSize());
		if(!strcmp(solSet, "population") || !strcmp(solSet, "both"))
			for(int i=0;i<sw.getSize();i++)
				if(!sw.particles[i].solution.dominated)
					repTemp.add(sw.particles[i].solution);
				
		if(!strcmp(solSet, "repository") || repTemp.getActualSize()==0 || !strcmp(solSet, "both"))
			repTemp.add(sw.repository);
	}
	
	repTemp.organize();
	
	
// 	if(!strcmp(archSubSwarms, "p-ideal")  || !strcmp(archSubSwarms, "w-ideal") || !strcmp(archSubSwarms, "r-ideal")){
	if(decomposition){
		updateWeightedDistances(repTemp.getSolutions(), repTemp.getActualSize(), sw.repository.weight);
		std::sort(repTemp.getSolutions(), repTemp.getSolutions()+repTemp.getActualSize(), weightedDistanceComparatorSol);
	}else{
		if(!strcmp(cmaRanking, "cd")){
			updateCrowdingDistances(repTemp.getSolutions(), repTemp.getActualSize());
			std::sort(repTemp.getSolutions(), repTemp.getSolutions()+repTemp.getActualSize(), crowdingComparatorSol);
		}
		if(!strcmp(cmaRanking, "hv")){
			updateContributingHypervolume(repTemp); //stores in the crowding distance field, since both are the higher the contribution, the better, there is no problem
			std::sort(repTemp.getSolutions(), repTemp.getSolutions()+repTemp.getActualSize(), crowdingComparatorSol);
		}
		if(!strcmp(cmaRanking, "r2")){
			if(refSize==-1){
				perror("ERROR ON CMA-ES R2 RANKING! Reference file not set.");
			}
			updateContributingR2(repTemp); //stores in the crowding distance field, stores in the crowding distance field //if the value of the r2 increase after I take this solution, it is important, so the higher is the contribution, the better
			std::sort(repTemp.getSolutions(), repTemp.getSolutions()+repTemp.getActualSize(), crowdingComparatorSol);
		}
	}
	
	//normalizing non dominated solutions
	double nonDominated[repTemp.getActualSize()][N];
	int nonDom=0;
// 	for(int i=0;i<repTemp.getActualSize()/2;i++){
	for(int i=0;i<repTemp.getActualSize();i++){
// 			printVector(repTemp.getSolution(i).decisionVector, decisionNumber);
// 		printVector(repTemp.getSolution(i).objectiveVector, objectiveNumber);
// 		printf("sol:%d: %e\n", i, repTemp.getSolution(i).crowdingDistance);
		for(int j=0;j<N;j++)
			nonDominated[nonDom][j]=normalize(repTemp.getSolution(i).decisionVector[j], inferiorPositionLimit[j], superiorPositionLimit[j]);
		nonDom++; //mu
	}
		
	int mu = nonDom;
	double weights[mu];
	
	//**************************** end of treatment of input data***************//	
	
	/***************setting default parameters according to****************/
	//https://www.lri.fr/~hansen/hansenedacomparing.pdf
	// in my strategy lambda = mu = solSize

	double weightSquare=0, weightSum=0, value=0;

	double smaller=MAXDOUBLE, larger=-MAXDOUBLE;
	if(mu > 2){
		if(!strcmp(weightDistribution, "metric")){ //uses the value of the ranking metric as weight
			for(int i=mu-1;i>=0;i--){// find the smaller > 0
				if(decomposition){
					value=repTemp.getSolution(i).weightedDistance;
					if(value < MAXDOUBLE){
						larger=value;
						break;
					}
				}
				else{
					value=repTemp.getSolution(i).crowdingDistance;
					if(value > 0){ //a little higher than 0
						smaller=value;
						break;
					}
				}
			}
			for(int i=0;i<mu;i++){ //find the larger < MAXDOUBLE
				if(decomposition){
					value=repTemp.getSolution(i).weightedDistance;
					if(value > 0){ //a little higher than 0
						smaller=value;
						break;
					}
				}
				else{
					value=repTemp.getSolution(i).crowdingDistance;
					if(value < MAXDOUBLE){
						larger=value;
						break;
					}
				}
			}
			if(larger == -MAXDOUBLE || smaller == MAXDOUBLE){
				fprintf(stderr, "\nERROR ON ASSIGNING WEIGHTS FOR CMA-ES!!\n");
				fprintf(stderr, "\nsm: %f lg: %f mu: %d\n", smaller, larger, mu);
				exit(1);
			}
		}
	}else{
		larger=1;
		smaller=0;
	}
		
	for(int i=0;i<mu;i++){
		if(!strcmp(weightDistribution, "equal"))
			weights[i]=1.0; //equal weight distribution
		if(!strcmp(weightDistribution, "log"))
			weights[i]=log(mu+1)-log(i+1.0); //according to code
		if(!strcmp(weightDistribution, "linear"))
			weights[i]=mu-i; //according to code
		if(!strcmp(weightDistribution, "metric")){ //uses the value of the ranking metric as weight
			if(decomposition){
				value=repTemp.getSolution(i).weightedDistance;
				weights[i]=1-(normalize(value, (smaller-(smaller/100.0)), larger));
			}else{
				value=repTemp.getSolution(i).crowdingDistance;
				weights[i]=(normalize(value, (smaller-(smaller/100.0)), larger));
			}
// 			weights[i]=exp(normalize(repTemp.getSolution(i).crowdingDistance, smaller, larger));
// 			if(repTemp.getSolution(i).crowdingDistance<=0)
// 				weights[i]=exp(0);
// 			if(repTemp.getSolution(i).crowdingDistance>=MAXDOUBLE)
// 				weights[i]=exp(1);
// 			weights[i]=(normalize(value, (smaller-(smaller/100.0)), larger));
			if(weights[i]<= 0)
				weights[i]=normalize(smaller, (smaller-(smaller/100.0)), larger);//weights[i-1]/2.0;
			if(weights[i]>1)
				weights[i]=1;
			
// 			weights[i]=log(1+weights[i]);
		}
		
		weightSum+=weights[i]; //sum of the weights for posterior normalization
	}
	
	for(int i=0;i<mu;i++){ //normalization of the weights
		weightSquare+=weights[i]*weights[i];
		weights[i]/=weightSum;
// 		printf("\n %f - %f",repTemp.getSolution(i).weightedDistance, weights[i]);
	}
	
// 	muEff=1.0/weightSquare;//paper
	muEff=weightSum*weightSum/weightSquare;//code

	cSigma= (muEff+2.0)/(N+muEff+3.0); //code
// 	cSigma= (muEff+2.0)/(N+muEff+5.0); //paper
	dSigma= (1.0 + 2.0*std::max(0.0, sqrt((muEff-1.0)/(N+1.0)) - 1.0)) * std::max(0.3, 1.0 -(double)N / (1e-6+maxIterations)) + cSigma;//code
// 	dSigma= 1.0+ 2.0*( std::max(0.0, sqrt( (muEff-1.0)/(N+1.0) )-1.0 )) +cSigma; //paper
	cc= 4.0/(N+4.0); //code
// 	cc= (4.0+(muEff/N))/ (N+4.0+(2.0*muEff/N)); //paper
	muCov=muEff;
// 	cCov=((1.0/muCov)*(2.0/( (N+sqrt(2.0))*(N+sqrt(2.0)) )))+((1.0- (1.0/muCov))*std::min(1.0, ((2.0*muEff)-1.0)/( ((N+2.0)*(N+2.0))+muEff )  )); //paper (probably more precise)
	double t1 = 2.0 / ((N+1.4142)*(N+1.4142));
	double t2 = (2.0*muEff-1.) / ((N+2.0)*(N+2.0)+muEff);
	t2 = (t2 > 1) ? 1 : t2;
	t2 = (1.0/muCov) * t1 + (1.0-1.0/muCov) * t2;
	cCov = t2; //code
	
	chiN = sqrt((double) N) * (1.0 - 1.0/(4.0*N) + 1./(21.0*N*N));
	/***************setting default parameters ****************/
	if(sw.init){
		for (int i = 0; i < N; ++i){  //avoid memory garbage
			for (int j = 0; j < N; j++)
			sw.B[i][j]=0.0;										//need to keep
			sw.B[i][i]=1.0;										//need to keep
			initialStds[i] = 0.3;									//does not change //code

			trace += initialStds[i]*initialStds[i];						//does not change
		}
		
		sw.sigma=sqrt(trace/N);										//need to keep
		for (int i = 0; i < N; ++i){  //avoid memory garbage
			for (int j = 0; j < N; j++)
				sw.C[i][j]=0.0;								//need to keep
			sw.pC[i] = sw.pSigma[i] = 0.0;						//need to keep
// 			sw.mean[i]=0.5;									//need to keep
// 			sw.mean[i]=normalize(sw.centroid[i], inferiorPositionLimit[i], superiorPositionLimit[i]);
			if(sw.gen>0)
				sw.mean[i]=(rand()/(double)RAND_MAX);//already initialized within the swarm, when restart does not need to restart, does it?
		}
		for (int i = 0; i < N; ++i) {
			sw.C[i][i] = sw.D[i] = initialStds[i] * sqrt(N / trace);	//need to keep
			sw.C[i][i] *= sw.C[i][i];							//need to keep
		}
		sw.gen=0;
		sw.init=false;
		if(logCMAES)
			printf("\n INIT \n");
		
// 		cmaes_init(&evo, decisionNumber, NULL, NULL, 0, mu, "cmaes_initials.par");
// 		printf("\nlambda: %d\n",evo.sp.lambda);
	}else
		sw.gen++;
		//exit(1);
		
	//******************************************//
	
	//http://arxiv.org/pdf/1110.4181v1.pdf  https://hal.inria.fr/hal-01146738/document //injecting external solutions into moead
	sw.solsKept=0;
	double clipped[mu][N], in[N], t[N];
	for (int i = 0; i < N; i++)
		oldMean[i]=sw.mean[i];
	for(int s=0;s<mu;s++){ //for all solutions
		bool injected=true;
		for(int r=0;r<sw.getSize();r++){//check if the current solution is in the population from the previous generation, otherwise it was injected
			for(int d=0;d<N;d++)
				t[d]=normalize(sw.particles[r].solution.decisionVector[d], inferiorPositionLimit[d], superiorPositionLimit[d]);
			
			if(isEqual(nonDominated[s],t,N)){//if this solution is in the population, it is not injected
				injected=false;
				sw.solsKept++;
			}
		}
		
		for (int i = 0; i < N; i++)// eq(2) turning x_i into y_i
			in[i]=(nonDominated[s][i]-oldMean[i])/sw.sigma;
		
		if(injected){//first step eq(3) : y_i <-- clip(cy,||C^{-1/2}y_i||)
			double sum, tmp[N], sum2=0;
			for (int i = 0; i < N; i++) {
				sum=0;
				for (int j = 0; j < N; ++j)
					sum += sw.B[j][i]*in[i];//B^T*y_i
				tmp[i] = sum / sw.D[i]; //B^T*y_i*D^-1
			}
			for (int i = 0; i < N; i++) {
				sum=0;
				for (int j = 0; j < N; ++j)
					sum += sw.B[i][j] * tmp[j];//B^T*y_i*B*D^-1
				sum2+=sum*sum;
			}
			sum2=sqrt(sum2);
			double cy=sqrt(N)+((2*N)/N+2);
			double clip=cy/sum2;	//clip(c,x) = min(1,c/x)
			clip=std::min(1.0,clip);
			
// 			double * in = nonDominated[s];//x_i
// 			double sum2=0;//C^{-1/2} * y_i
// 			for (int i = 0; i < N; ++i) {
// 				double sum = 0.0;
// 				for (int j = 0; j < N; ++j)
// 					sum += sw.B[i][j] * sw.D[i] * ((in[i]-oldMean[i])/sw.sigma ); //y_i
// 				sum2+=sum*sum;
// 			}
// 			sum2=sqrt(sum2);
// 			double cy=sqrt(N)+((2*N)/N+2);
// 			double clip=cy/sum2;	//clip(c,x) = min(1,c/x)
// 			clip=std::min(1.0,clip);
// 			printf(" -- %f %f %f\n----------------------\n",sum2, cy, clip);
			for (int i = 0; i < N; ++i)
				clipped[s][i]=in[i]*clip;
		}else
			for (int i = 0; i < N; ++i)
				clipped[s][i]=in[i];
	}//end of first step, clipping y_i
	//second step: \Delta m
	//As far as I understood, the differentiation in eq (4) is optional, to be used only if we want a strong impact with an injection, it is not used in the MOEA/D-CMA-ES paper, and is not going to be used here.
	double deltaM[N];
	for (int i = 0; i < N; ++i){ //update of the delta mean using the second part of eq (4) 
		deltaM[i]=0;
		for(int j=0;j<mu;++j)
			deltaM[i]+=weights[j]*clipped[j][i];
	}//end of delta m update
	//mean update eq(5)
	double cm=1;
	for (int i = 0; i < N; ++i)
		sw.mean[i]=oldMean[i]+cm*sw.sigma*deltaM[i];
	//end of mean update
// // 	//delta mean clip eq(6), but only if strong injection of mean shift, we do not use this
// // 	if(injected){
// // 		double * in = deltaM;//deltaM
// // 		double sum, tmp[N], sum2=0;
// // 		for (int i = 0; i < N; i++) {
// // 			sum=0;
// // 			for (int j = 0; j < N; ++j)
// // 				sum += sw.B[j][i]*in[i];//B^T*deltaM
// // 				tmp[i] = sum / sw.D[i]; //B^T*deltaM*D^-1
// // 		}
// // 		for (int i = 0; i < N; i++) {
// // 			sum=0;
// // 			for (int j = 0; j < N; ++j)
// // 				sum += sw.B[i][j] * tmp[j];//B^T*deltaM*B*D^-1
// // 				sum2+=sum*sum;
// // 		}
// // 		sum2=sqrt(sum2);
// // 		double cym=sqrt(2*N)+((2*N)/N+2);
// // 		double clip=cym/(sqrt(muEff)*sum2);	//clip(c,x) = min(1,c/x)
// // 		clip=std::min(1.0,clip);
// // 		for (int i = 0; i < N; ++i)
// // 			deltaM[i]*=clip;//if not injected, deltaM keeps the same
// // 	}
// // 	//end of the delta mean clip
//Equations 7 and 9 were updated by changing the BDz calculation, now it is sqrt(muEff)*deltaM
//Equation 8 is a little different, but there is no apparent reason for it, so we did not change it yet
//Equation 10 is changed in the cov matrix
//Equation 11 is changed in the sigma update
//------------------------------------------------------------//
	
// 	for (int i = 0; i < N; ++i) {//update of the mean and BDz //original CMA-ES update
// 		oldMean[i]=sw.mean[i];
// 		sw.mean[i]=0;
// 		for(int j=0;j<mu;++j)
// 			sw.mean[i]+=weights[j]*nonDominated[j][i]; // in the paper sw.mean = m ;; y_w = (mean-oldMean)/sigma
// 		BDz[i] = sqrt(muEff)*(sw.mean[i] - oldMean[i])/sw.sigma; // BDz * sqrt(muEff) = y_w * sqrt(muEff)
// 	}
	for (int i = 0; i < N; ++i)//BDz with deltaM, as in the paper of injecting solutions into CMA-ES
		BDz[i]=sqrt(muEff)*deltaM[i];
	
// 	if(decomposition && sw.neighborhood[0].index==swarmNumber/2){//if the first neighbor of it (itself) is zero
// 		double sum=0;
// 		for(int i=0;i<N;i++)
// 			sum+=abs(sw.mean[i]-oldMean[i]);
// 		printf("%d -- %f\n",repTemp.getActualSize(), sw.sigma);
// 		for(int i=0;i<mu;i++){
// // 			printf("%f -- ", repTemp.getSolution(i).weightedDistance);
// 			for(int d=0;d<N;d++){
// 				printf("%.4f ", nonDominated[i][d]);
// // 				printf("%.4f ", repTemp.getSolution(i).decisionVector[d]);
// 			}
// 			printf(" (%.4f)\n", getScalarValue(repTemp.getSolution(i).objectiveVector, sw.repository.weight));
// 		}
// 		printf("-----------------------------------------\ncovMat:\n");
// 		for(int i=0;i<N;i++){
// 			for(int j=0;j<N;j++)
// 				printf("%.4f ",sw.C[i][j]);
// 			printf("\n");
// 		}
// 		printf("---------------------------------------\n");
// 		printf("mean: ");
// 		for(int d=0;d<N;d++){
// 			// 			printf("%.4f ", sw.mean[i]);
// 			printf("%.4f ", sw.mean[d]);
// 		}
// 		printf("\nbest: ");
// 		for(int d=0;d<N;d++){
// 			printf("%.4f ", normalize(sw.repository.getSolution(0).decisionVector[d], inferiorPositionLimit[d], superiorPositionLimit[d]) );
// 		}
// 		printf("\n---------------------------------------------------\n");
// 	}
	
	double sum, tmp[N], psxps=0.0;
	//calculate z := D^(-1) * B^(-1) * BDz into tmp // orignal comment
	for (int i = 0; i < N; i++) {
		sum=0;
		for (int j = 0; j < N; ++j)
			sum += sw.B[j][i] * BDz[j];
		tmp[i] = sum / sw.D[i]; //tmp = z = D^(-1) * B' * BDz
	}
	
	//cumulation for sigma (ps) using B*z //original comment
	for (int i = 0; i < N; i++) {
		sum=0;
		for (int j = 0; j < N; ++j)
			sum += sw.B[i][j] * tmp[j]; //sum = sqrt(muEff)*C^(1/2)*y_w = sqrt(muEff)*BD^−1B'*y_w
		
		sw.pSigma[i] = (1.0 - cSigma) * sw.pSigma[i] + sqrt(cSigma * (2.0 - cSigma)) * sum;
		
		/* calculate norm(ps)^2 */
		psxps += sw.pSigma[i] * sw.pSigma[i];
	}

	/* cumulation for covariance matrix (pc) using B*D*z~N(0,C) */
	hsig = sqrt(psxps) / sqrt(1.0 - pow(1.0-cSigma, 2*(sw.gen+1))) / chiN < 1.4 + 2.0/(N+1);
	
	for (int i = 0; i < N; ++i)
		sw.pC[i] = (1.0 - cc) * sw.pC[i] + hsig * sqrt(cc * (2.0 - cc)) * BDz[i];
	
	//************  Adapt_C2  *********// Update covariance matrix
	double ccov1 = std::min(cCov * (1.0/muCov) * 1.0, 1.0); //code
// 	double ccov1 = 2.0/(((N+1.3)*(N+1.3))+muEff); //paper
	double ccovmu = std::min(cCov * (1-1.0/muCov)* 1.0, 1.0-ccov1); //code
// 	double ccovmu = std::min(1.0-ccov1, 2 * ( ( (muEff-2.0)+(1.0/muEff) )/( (N+2.0)*(N+2.0) + (2.0*muEff)/2.0 ) )); //paper
	double sigmasquare = sw.sigma * sw.sigma;
	
	/* update covariance matrix */
	for (int i = 0; i < N; ++i)
		for (int j = 0; j <= i; ++j) {
			sw.C[i][j] = (1.0 - ccov1 - ccovmu) * sw.C[i][j] + ccov1* (sw.pC[i] * sw.pC[j] + (1.0-hsig)*cc*(2.0-cc) * sw.C[i][j]);
			for (int k = 0; k < mu; ++k) { /* additional rank mu update */
// 				sw.C[i][j] += ccovmu * weights[k]* (nonDominated[k][i] - oldMean[i]) * (nonDominated[k][j] - oldMean[j])/ sigmasquare;//original equation
				sw.C[i][j] += ccovmu * weights[k] * clipped[k][i] * clipped[k][j];//changed on eq (10) of the injecting solutions into cmaes paper
			}
		}
	//**************************************//
	
// 	dSigma= (1.0 + 2.0*std::max(0.0, sqrt((muEff-1.0)/(N+1.0)) - 1.0)) * std::max(0.3, 1.0 -(double)N / (1e-6+std::min(evo.sp.stopMaxIter, evo.sp.stopMaxFunEvals/evo.sp.lambda))) + cSigma;//code (used to compare our algorithm to the implemented)
	
	/* update of sigma */
// 	sw.sigma *= exp(((sqrt(psxps)/chiN)-1.0)*cSigma/dSigma);//original equation
	double deltaSigmaMax=1;
	sw.sigma *= exp(std::min(deltaSigmaMax, ((sqrt(psxps)/chiN)-1.0)*cSigma/dSigma ) );//Changed on eq (11) of the injecting solutions into cmaes paper
	
	if(logCMAES){
		printf("\n\nmueff %f, cSigma %f, dSigma %f, ccumcov %f, sigma %f", muEff, cSigma, dSigma, cc, sw.sigma);
		printf("\nhsig %d, psxps %f, chiN %f, gen %d \n", hsig, psxps, chiN, sw.gen);
	// 	printf("\n\nweights - mu: %d : ", mu);
	// 	printVector(weights, mu);
		
		printf("\navg_old : ");
		printVector(oldMean, N);
		printf("\navg :     ");
		printVector(sw.mean, N);

// 		printf("\n\nD mat:");
// 		printVector(sw.D, N);
		printf("\n\nBDz: ");
		printVector(BDz, N);
		printf("\ntmp : ");
		printVector(tmp, N);
		printf("\nps : ");
		printVector(sw.pSigma, N);
		printf("\npc : ");
		printVector(sw.pC, N);
	}
// 		t->rgBDz[i] = sqrt(t->sp.mueff)*(t->rgxmean[i] - t->rgxold[i])/t->sigma; 
	
// 	printf("\n\nrgD (new): ");
// 	printVector(D, N);
// 	printf("\n\nB (new): \n");
// 	for (int i = 0; i < N; ++i){
// 		printVector(B[i], N);
// 		printf("\n");
// 	}
// 	printf("\n\nCov mat: \n");
// 	for (int i = 0; i < N; ++i){
// 		printVector(sw.C[i], N);
// 		printf("\n");
// 	}

//****************comparing the implemented CMA-ES with our version**********************//
// evo.sp.mu=mu;
// evo.sp.lambda=mu;
// 
// for(int i=0;i<mu;i++)
// 	for(int j=0;j<N;j++)
// 		evo.rgrgx[i][j]=nonDominated[i][j];
// cmaes_SamplePopulation(&evo); /* do not change content of pop */
// for(int i=0;i<mu;i++)
// 	for(int j=0;j<N;j++)
// 		evo.rgrgx[i][j]=nonDominated[i][j];
// 
// for(int i=0;i<N;i++)
// 	evo.rgxmean[i]=oldMean[i];
// 
// printf("\nccov1 %f", ccov1);
// printf("\nccovmu %f", ccovmu);
// printf("\nsigmasquare %.20f", sigmasquare);
// printf("\nhsig %d", hsig);
// printf("\nccumcov %f", cc);
// printf("\nccov %.20f", cCov);
// printf("\nmucov %.20f", muCov);
// printf("\nmueff: %.20f", muEff);
// printf("\nBDz: ");
// printVector(BDz, N);
// printf("\npSigma: ");
// printVector(sw.pSigma, N);
// printf("\npC: ");
// printVector(sw.pC, N);
// printf("\nD: ");
// for(int i=0;i<N;i++)
// 	printf("%.15f ",sw.D[i]);
// printf("\npsxps: %f", psxps);
// printf("\nchin: %f", chiN);
// printf("\ncSigma: %f", cSigma);
// printf("\ndSigma: %f", dSigma);
// printf("\nsigma: %f", sw.sigma);
// printf("\nmean: ");
// printVector(sw.mean, N);
// printf("\nold: ");
// printVector(oldMean, N);
// printf("\nweights: ");
// printVector(weights, mu);
// 	
// cmaes_UpdateDistribution(&evo, weights);
// 
// if(evo.sp.ccov != cCov){
// 	printf("\ncCov different\n");
// 	exit(1);
// }
// if(evo.sp.mucov != muCov){
// 	printf("\nmuCov different\n");
// 	exit(1);
// }
// if(evo.sp.ccumcov != cc){
// 	printf("\ncc different\n");
// 	exit(1);
// }
// if(evo.sp.mueff != muEff){
// 	printf("\nmuEff different\n");
// 	exit(1);
// }
// if(evo.sp.cs != cSigma){
// 	printf("\ncSigma different\n");
// 	exit(1);
// }
// if(evo.sp.damps != dSigma){
// 	printf("\ndSigma different\n");
// 	exit(1);
// }
// for(int i=0;i<N;i++){
// 	if(evo.rgBDz[i] != BDz[i]){
// 		printf("\nBDz on %d different\n", i);
// 		exit(1);
// 	}
// 	if(evo.rgD[i] != sw.D[i]){
// 		printf("\nD on %d different\n", i);
// 		exit(1);
// 	}
// 	if(evo.rgdTmp[i] != tmp[i]){
// 		printf("\ntmp on %d different\n", i);
// 		exit(1);
// 	}
// 	if(evo.rgps[i] != sw.pSigma[i]){
// 		printf("\npSigma on %d different (%.20f)\n", i, abs(evo.rgps[i] - sw.pSigma[i]));
// 		exit(1);
// 	}
// 	if(evo.rgpc[i] != sw.pC[i]){
// 		printf("\npC on %d different\n", i);
// 		exit(1);
// 	}
// }
// 
// printf("\nevo mean: ");
// printVector(evo.rgxmean, N);
// printf("\nevo old: ");
// printVector(evo.rgxold, N);
// printf("\nevo weights: ");
// printVector(evo.sp.weights, mu);
// 
// printf("\n\norig-cmaes\n");
// for(int i=0; i<cmaes_Get(&evo, "dim");i++){
// 	for(int j=0;j<cmaes_Get(&evo, "dim");j++)
// 		printf("%f ",evo.C[i][j]);
// 	printf("\n");
// }
// 
// printf("\nmodif-cmaes\n");
// for(int i=0; i<N;i++){
// 	for(int j=0;j<N;j++)
// 		printf("%f ",sw.C[i][j]);
// 	printf("\n");
// }
// 
// for(int i=0; i<N;i++){
// 	for(int j=0;j<N;j++){
// 		if(((int)(evo.C[i][j]*1000000))/1000000.0 != ((int)(sw.C[i][j]*1000000))/1000000.0){
// 			printf("Different %f -- %f \t (%e) (%d,%d)\n", evo.C[i][j], sw.C[i][j], abs(evo.C[i][j] - sw.C[i][j]), i, j );
// 			maxDifference=std::max(maxDifference, abs(evo.C[i][j] - sw.C[i][j]));
// // 			exit(1);
// 		}
// 	}
// }
// printf("maxDifference: %.20f, sigmaDifference: %.20f\n", maxDifference, abs(evo.sigma-sw.sigma));
// if(maxDifference > 1 ){
// 	printf("%d\n", mu);
// 	exit(1);
// }

//**************************************************************************************//

}

bool mySample (Swarm &sw){
	int N=decisionNumber;
	/* calculate eigensystem */

	if(logCMAES){
		printf("\n----------------------before eigen calculation ---------------------- ");
		
// 		printf("\nCov mat: \n");
// 		for (int i = 0; i < N; ++i){
// 			printVector(sw.C[i], N);
// 			printf("\n");
// 		}
		
		printf("\nB mat: \n");
		for (int i = 0; i < N; ++i){
			printVector(sw.B[i], N);
			printf("\n");
		}
		
		printf("\nD mat: \n");
		printVector(sw.D, N);
		printf("\n");
		
		printf("---------------------end before ---------------------- \n");
		
		printf("\nCov mat: \n");
		for (int i = 0; i < N; ++i){
			printVector(sw.C[i], N);
			printf("\n");
		}
	}	
	
	double rgdTmp[N];
	Eigen( N, sw.C, sw.D, sw.B, rgdTmp);
	//By using Cholesky decomposition, it turns out that the determinant of a matrix is equal to the product of its eigenvalues
	sw.det=sw.D[0];
	for(int i=1;i<N;i++)
		sw.det*=sw.D[i];
	
	double lg=log(sw.det);
	if(lg >= MAXDOUBLE || lg != lg){
		fprintf(stderr, "WARNING!, log of the covariance matrix determinant: %f. ", lg); //shows warning message
		sw.init=true; //set the init flag, so all the CMA-ES variables are reset
		return false;
	}
	
	for (int i = 0; i < N; ++i)
		sw.D[i] = sqrt(sw.D[i]);
	
	for(int i=0;i<N;i++){
		if(sw.D[i] != sw.D[i] || sw.D[i] >= MAXDOUBLE || sw.D[i] < 0){//tests for 1)null 2)inf 3) positive definiteness of the covariance matrix
			// 		if(sw.det != sw.det || sw.det <= MAXDOUBLE*-1 || sw.D[i] != sw.D[i] || sw.D[i] >= MAXDOUBLE || sw.D[i] < 0){//tests for 1) determinant null or invalid 2)null 3)inf 4) positive definiteness of the covariance matrix
			
			fprintf(stderr, "WARNING!, value: %f in eigenValues vector. ", sw.D[i]); //shows warning message
			sw.init=true; //set the init flag, so all the CMA-ES variables are reset
			return false;
			// 			exit(1);
			//myLearn(sw, neighboringSwarms); //learn again using the default values as base
		}
	}
	
	if(logCMAES){
		printf("---------------------after eigen calculation ---------------------- \n");
		
// 		printf("\nCov mat: \n");
// 		for (int i = 0; i < N; ++i){
// 			printVector(sw.C[i], N);
// 			printf("\n");
// 		}
		
		printf("B mat: \n");
		for (int i = 0; i < N; ++i){
			printVector(sw.B[i], N);
			printf("\n");
		}
		
		printf("\nD mat: \n");
		printVector(sw.D, N);
		printf("\n");
		
		printf("---------------------end after ---------------------- \n");
		
		printf("\nsigma: %f (%e)\n",sw.sigma, sw.sigma);
		
		for (int i = 0; i < N; ++i)
			rgdTmp[i] = sw.D[i] * 3;
		printf("\nD (sample) (* 3) : ");
		for (int i = 0; i < N; ++i){
			double sum = 0.0;
			for (int j = 0; j < N; ++j)
				sum += sw.B[i][j] * rgdTmp[j];
			if((sw.mean[i] + sw.sigma * sum ) >= 0 )
				printf(" ");
			printf("%.3f ", /*fabs*/(sw.mean[i] + sw.sigma * sum ) );
		}
		
		for (int i = 0; i < N; ++i)
			rgdTmp[i] = sw.D[i] * 0;
		printf("\nD (sample) (* 0) : ");
		for (int i = 0; i < N; ++i){
			double sum = 0.0;
			for (int j = 0; j < N; ++j)
				sum += sw.B[i][j] * rgdTmp[j];
			if((sw.mean[i] + sw.sigma * sum ) >= 0 )
				printf(" ");
			printf("%.3f ", /*fabs*/(sw.mean[i] + sw.sigma * sum ) );
		}
		
		for (int i = 0; i < N; ++i)
			rgdTmp[i] = sw.D[i] * -3;
		printf("\nD (sample) (* -3): ");
		for (int i = 0; i < N; ++i){
			double sum = 0.0;
			for (int j = 0; j < N; ++j)
				sum += sw.B[i][j] * rgdTmp[j];
			if((sw.mean[i] + sw.sigma * sum ) >= 0 )
				printf(" ");
			printf("%.3f ", /*fabs*/(sw.mean[i] + sw.sigma * sum ) );
		}
	}
	

// // 	printVector(sw.mean, decisionNumber); //show the average decision vector per iteration
// // 	printf("\n");
	
	
// 	printf("\nD (sample) (difsum): ");
// 	for (int i = 0; i < N; ++i){
// 		for (j = 0, sum = 0.0; j < N; ++j)
// 			sum += sw.B[i][j] * rgdTmp[j];
// 		printf("%f ", fabs(sw.mean[i] + sw.sigma * sw.D[i]*-3 )+fabs(sw.mean[i] + sw.sigma * sqrt(sw.D[i])*3 ) );
// 	}
	
	/*************************/

	double previousSol[N];
	if(sw.getSize() == 1)
		memcpy(previousSol, sw.particles[0].solution.decisionVector, sizeof(double)*N);

// 	original sampling
	for (int iNk = 0; iNk < sw.getSize(); ++iNk){ /* generate scaled cmaes_random vector (D * z)    */
	// 	int rescount=0;
	// 	while(true){//resample when invalid
	// 		bool resample=false;
		
			for (int i = 0; i < N; ++i)
				rgdTmp[i] = sw.D[i] * cmaes_random_Gauss();
			/* add mutation (sigma * B * (D*z)) */
			for (int i = 0; i < N; ++i) {
				double sum = 0.0;
				for (int j = 0; j < N; ++j)
					sum += sw.B[i][j] * rgdTmp[j];

				double value=sw.mean[i] + sw.sigma * sum;//normalized new solution
				
	// 			if(value < 0 || value > 1){
	// 				resample=true;
	// // 				printf("\nresample: %f (%d, %d)\n", value, i, rescount);
	// 			}else
					sw.particles[iNk].solution.evalSolution=true;
					sw.particles[iNk].solution.decisionVector[i] = (value*(superiorPositionLimit[i]-inferiorPositionLimit[i]))+inferiorPositionLimit[i]; //unormalize solution
			}
	// 		if(!resample || rescount > 100)
	// 			break;
	// 		else
	// 			rescount++;
	// 	}
	// 	if(rescount >= 100){
	// 		printf("\n GAVE UP \n");
	// 		exit(1);
	// 	}
	}
	
// 	//ranking for the new surrogate functions
// 	if(repGlobal->getActualSize() >0 ){
// 		repGlobal->organize();
// 		for(int i=0;i<repGlobal->getActualSize();i++)
// 			repGlobal->getSolutions()[i].crowdingDistance=logLikelihood(repGlobal->getSolution(i).decisionVector, sw);//stores temporarily in the crowding distance field
// 		
// 		std::sort(repGlobal->getSolutions(), repGlobal->getSolutions()+repGlobal->getActualSize(), crowdingComparatorSol);
// 		
// 		
// 		double likelihood=logLikelihood(sw.particles[0].solution.decisionVector, sw);
// 		
// 		double rank = likelihoodRanking(sw.particles[0].solution, *repGlobal, sw);
// 		
// 		if(rank >0 && rank < repGlobal->getActualSize())
// 			printf("rank: %.1f lik before: %f lik: %f lik after: %f\n", rank, repGlobal->getSolution((int)rank).crowdingDistance,likelihood,repGlobal->getSolution((int)rank+1).crowdingDistance);
// 		else
// 			printf("rank: %f\n", rank);
// 	}
	
	
	
	
	
// // // 	for (int iNk = 0; iNk < sw.getSize(); ++iNk){//using the likelihood as surrogate
// // // // 		double best[N], bestValue=MAXDOUBLE*-1;//using best as highest likelihood to get convergence
// // // 		double best[N], bestValue=MAXDOUBLE;//using the best as lowest likelihood to get diversity
// // // 		for(int t=0;t<10;t++){//ten trials per particle and choose the best
// // // 			double candidate[N];
// // // 			
// // // 			for (int i = 0; i < N; ++i)
// // // 				rgdTmp[i] = sw.D[i] * cmaes_random_Gauss();
// // // 			/* add mutation (sigma * B * (D*z)) */
// // // 			for (int i = 0; i < N; ++i) {
// // // 				double sum = 0.0;
// // // 				for (int j = 0; j < N; ++j)
// // // 					sum += sw.B[i][j] * rgdTmp[j];
// // // 
// // // 				candidate[i]=sw.mean[i] + sw.sigma * sum;//normalized new solution
// // // 			}
// // // 			double opinion=neighborhoodOpinion(sw, candidate);
// // // // 			printf("opinion: %f of: ",opinion);
// // // // 			printVector(candidate, N);
// // // // 			printf("\n");
// // // // 			if(opinion > bestValue){//using best as highest likelihood to get convergence
// // // 			if(opinion < bestValue){//using the best as lowest likelihood to get diversity
// // // 				for (int i = 0; i < N; ++i)
// // // 					best[i]=candidate[i];
// // // 				bestValue=opinion;
// // // 			}
// // // 		}
// // // 		
// // // 		for (int i = 0; i < N; ++i)
// // // 			sw.particles[iNk].solution.decisionVector[i] = (best[i]*(superiorPositionLimit[i]-inferiorPositionLimit[i]))+inferiorPositionLimit[i]; //unormalize solution
// // // 	}
	
	
	/* Test if function values are identical, escape flat fitness */ //original comment
	//Check if the sampled solutions are equal, or if only one solution is sampled, if it is equal to the previous solution
// 	if (t->rgFuncValue[t->index[0]] == t->rgFuncValue[t->index[(int)t->sp.lambda/2]]) {
// 		t->sigma *= exp(0.2+t->sp.cs/t->sp.damps);
// 		ERRORMESSAGE("Warning: sigma increased due to equal function values\n",
// 				   "   Reconsider the formulation of the objective function",0,0);
// 	}
	if(sw.getSize() == 1){
		if(isEqual(previousSol, sw.particles[0].solution.decisionVector, N)){
			sw.sigma *= exp(0.2+cSigma/dSigma);
			fprintf(stderr, "Warning: sigma increased due to equal decision vectors.\n");
			return false;
		}
// 		else
// 			for(int i=0;i<N;i++)
// 				printf("%.14f ", abs(previousSol[i]-sw.particles[0].solution.decisionVector[i]));
// 			printf("\n");
	}else{
		if(isEqual(sw.particles[0].solution.decisionVector, sw.particles[sw.getSize()-1].solution.decisionVector, N)){
			sw.sigma *= exp(0.2+cSigma/dSigma);
			fprintf(stderr, "Warning: sigma increased due to equal decision vectors.\n");
			return false;
		}
	}

	return true;
} /* SamplePopulation() */

bool mySampleReplacement (Swarm &sw, int trials){
	int N=decisionNumber;
	
	//-------------------------------------------------------------------------------------------
	
	/* calculate eigensystem */
	double rgdTmp[N];
	Eigen( N, sw.C, sw.D, sw.B, rgdTmp);
	//By using Cholesky decomposition, it turns out that the determinant of a matrix is equal to the product of its eigenvalues
	sw.det=sw.D[0];
	for(int i=1;i<N;i++)
		sw.det*=sw.D[i];
	
	double lg=log(sw.det);
	if(lg >= MAXDOUBLE || lg != lg){
		fprintf(stderr, "WARNING!, log of the covariance matrix determinant: %f. ", lg); //shows warning message
		sw.init=true; //set the init flag, so all the CMA-ES variables are reset
		return false;
	}
	
	for (int i = 0; i < N; ++i)
		sw.D[i] = sqrt(sw.D[i]);
	
	for(int i=0;i<N;i++){
		if(sw.D[i] != sw.D[i] || sw.D[i] >= MAXDOUBLE || sw.D[i] < 0){//tests for 1)null 2)inf 3) positive definiteness of the covariance matrix
			// 		if(sw.det != sw.det || sw.det <= MAXDOUBLE*-1 || sw.D[i] != sw.D[i] || sw.D[i] >= MAXDOUBLE || sw.D[i] < 0){//tests for 1) determinant null or invalid 2)null 3)inf 4) positive definiteness of the covariance matrix
			
			fprintf(stderr, "WARNING!, value: %f in eigenValues vector. ", sw.D[i]); //shows warning message
			sw.init=true; //set the init flag, so all the CMA-ES variables are reset
			return false;
			// 			exit(1);
			//myLearn(sw, neighboringSwarms); //learn again using the default values as base
		}
	}
	
	//-------------------------------------------------------------------------------------------

	double previousSol[N];
	if(sw.getSize() == 1)
		memcpy(previousSol, sw.particles[0].solution.decisionVector, sizeof(double)*N);
	
// 	//ranking for the new surrogate functions
// 	if(repGlobal->getActualSize() >0 ){
// 		repGlobal->organize();
// 		for(int i=0;i<repGlobal->getActualSize();i++)
// 			repGlobal->getSolutions()[i].crowdingDistance=logLikelihood(repGlobal->getSolution(i).decisionVector, sw);//stores temporarily in the crowding distance field
// 		
// 		std::sort(repGlobal->getSolutions(), repGlobal->getSolutions()+repGlobal->getActualSize(), crowdingComparatorSol);
// 		
// 		
// 		double likelihood=logLikelihood(sw.particles[0].solution.decisionVector, sw);
// 		
// 		double rank = likelihoodRanking(sw.particles[0].solution, *repGlobal, sw);
// 		
// 		if(rank >0 && rank < repGlobal->getActualSize())
// 			printf("rank: %.1f lik before: %f lik: %f lik after: %f\n", rank, repGlobal->getSolution((int)rank).crowdingDistance,likelihood,repGlobal->getSolution((int)rank+1).crowdingDistance);
// 		else
// 			printf("rank: %f\n", rank);
// 	}
	
	for (int iNk = 0; iNk < sw.getSize(); ++iNk){//using the likelihood as surrogate
		double best[N], bestValue=MAXDOUBLE*-1;//using best as highest likelihood to get convergence
// 		double best[N], bestValue=MAXDOUBLE;//using the best as lowest likelihood to get diversity
		for(int t=0;t<trials;t++){//trials per particle and choose the best
			double candidate[N];
			
			for (int i = 0; i < N; ++i)
				rgdTmp[i] = sw.D[i] * cmaes_random_Gauss();
			/* add mutation (sigma * B * (D*z)) */
			for (int i = 0; i < N; ++i) {
				double sum = 0.0;
				for (int j = 0; j < N; ++j)
					sum += sw.B[i][j] * rgdTmp[j];

				candidate[i]=sw.mean[i] + sw.sigma * sum;//normalized new solution
			}
			double opinion=neighborhoodOpinion(sw, candidate);
// 			printf("opinion: %f of: ",opinion);
// 			printVector(candidate, N);
// 			printf("\n");
			if(opinion > bestValue){//using best as highest likelihood to get convergence
// 			if(opinion < bestValue){//using the best as lowest likelihood to get diversity
				for (int i = 0; i < N; ++i)
					best[i]=candidate[i];
				bestValue=opinion;
			}
		}
		
		sw.particles[iNk].solution.evalSolution=true;
		for (int i = 0; i < N; ++i)
			sw.particles[iNk].solution.decisionVector[i] = (best[i]*(superiorPositionLimit[i]-inferiorPositionLimit[i]))+inferiorPositionLimit[i]; //unormalize solution
	}

	if(sw.getSize() == 1){
		if(isEqual(previousSol, sw.particles[0].solution.decisionVector, N)){
			sw.sigma *= exp(0.2+cSigma/dSigma);
			fprintf(stderr, "Warning: sigma increased due to equal decision vectors.\n");
			return false;
		}
// 		else
// 			for(int i=0;i<N;i++)
// 				printf("%.14f ", abs(previousSol[i]-sw.particles[0].solution.decisionVector[i]));
// 			printf("\n");
	}else{
		if(isEqual(sw.particles[0].solution.decisionVector, sw.particles[sw.getSize()-1].solution.decisionVector, N)){
			sw.sigma *= exp(0.2+cSigma/dSigma);
			fprintf(stderr, "Warning: sigma increased due to equal decision vectors.\n");
			return false;
		}
	}
	return true;
}

void cma_es(Swarm &sw){
	
	if(logCMAES){
		printf("\n--------------------------------------------------------------\n");
		printf("Repository (%d)\n", sw.repository.getActualSize());
		for(int i=0;i<sw.repository.getActualSize();i++){
			printVector(sw.repository.getSolution(i).decisionVector, decisionNumber);
			printf("\n");
		}
	}
	
	double sigmaPrev=sw.sigma;
	myLearn(sw);
	
	// 		if(fabs(sigmaPrev-sw.sigma) > 10000){
	// 			fprintf(stderr, "\nWARNING!, sigma changed too much: %f -> %f. resetting...\n",sigmaPrev, sw.sigma); //shows warning message
	// 			sw.init=true; //set the init flag, so all the CMA-ES variables are reset
	// 			myLearn(sw); //learn again using the default values as base
	// 		}
	
	if(sw.sigma != sw.sigma || sw.sigma >= MAXDOUBLE){ //check for invalid numbers NaN or inf
		fprintf(stderr, "WARNING!, sigma: %f. resetting...\n", sw.sigma); //shows warning message
		// 		exit(1);
		sw.init=true; //set the init flag, so all the CMA-ES variables are reset
		myLearn(sw); //learn again using the default values as base
	}
	
	//four reset criteria as in the paper "injecting cma-es into moea/d"
	//NoEffectCoord
	for (int iKoo = 0; iKoo < decisionNumber; ++iKoo){
		if (sw.mean[iKoo] == sw.mean[iKoo] + 0.2*sw.sigma*sqrt(sw.C[iKoo][iKoo])){
			fprintf(stderr, "NoEffectCoordinate: standard deviation 0.2*%7.2e in coordinate %d without effect\n", sw.sigma*sqrt(sw.C[iKoo][iKoo]), iKoo); //shows warning message
			sw.init=true; //set the init flag, so all the CMA-ES variables are reset
			myLearn(sw); //learn again using the default values as base
			break;
		}
	}
	//NoEffectAxis
	int iKoo;
	for (int iAchse = 0; iAchse < decisionNumber; ++iAchse){
		double fac = 0.1 * sw.sigma * sw.D[iAchse];
		for (iKoo = 0; iKoo < decisionNumber; ++iKoo){ 
			if (sw.mean[iKoo] != sw.mean[iKoo] + fac * sw.B[iKoo][iAchse])
				break;
		}
		if (iKoo == decisionNumber){
			/* t->sigma *= exp(0.2+t->sp.cs/t->sp.damps); */
			fprintf(stderr, "NoEffectAxis: standard deviation 0.1*%7.2e in principal axis %d without effect\n", fac/0.1, iAchse);
			sw.init=true; //set the init flag, so all the CMA-ES variables are reset
			myLearn(sw); //learn again using the default values as base
			break;
		}
	}
	//TolXUp
	double stopTolUpXFactor=1e3;
	double initialStds=0.3;
	int i;
	for(i=0; i<decisionNumber; ++i) {
		if (sw.sigma * sqrt(sw.C[i][i]) > stopTolUpXFactor * initialStds)
			break;
	}
	if (i < decisionNumber) {
		fprintf(stderr, "TolUpX: standard deviation increased by more than %7.2e, larger initial standard deviation recommended \n", stopTolUpXFactor);
		sw.init=true; //set the init flag, so all the CMA-ES variables are reset
		myLearn(sw); //learn again using the default values as base
	}
	//ConditionCov
	double dMaxSignifKond=1e13;
	if (rgdouMax(sw.D, decisionNumber) >= rgdouMin(sw.D, decisionNumber) * dMaxSignifKond) {
		fprintf(stderr, "ConditionNumber: maximal condition number %7.2e reached. maxEW=%7.2e,minEW=%7.2e\n", dMaxSignifKond, rgdouMax(sw.D, decisionNumber), rgdouMin(sw.D, decisionNumber));
		sw.init=true; //set the init flag, so all the CMA-ES variables are reset
		myLearn(sw); //learn again using the default values as base
	}
	
	
	//************************************* RESAMPLE SOLUTIONS FROM THE GOOD FRONTS ******************************//
	//re-learn and re-sample when errors in the covariance matrix are detected in the sampling phase
	bool success=mySample(sw);
	while(!success){
		myLearn(sw);
		success=mySample(sw);
		fprintf(stderr,"Resample\n");
	}
	
	
// 	Neighbor orderedSwarms[swarmNumber];
// 	for(int i=0;i<swarmNumber;i++){
// 		orderedSwarms[i].index=i;
// 		orderedSwarms[i].distance=swarms[i].modelQuality*-1;//inverted model quality because we will use the preexisting sort, that sorts based on the smallest value (distance)
// 	}
// 	std::sort(orderedSwarms, orderedSwarms+swarmNumber, neighborsComparator);
// 	
// // 	for(int i=0;i<swarmNumber;i++){
// // 		printf("swarm: %d has quality %f sz: %d\n", orderedSwarms[i].index, orderedSwarms[i].distance*-1, repGlobal->getActualSize());
// // 	}
// // 	printf("\n\n\n");
// 
// 	bool good=false;
// 	for(int i=0;i<100;i++){
// 		if(sw.neighborhood[0].index == orderedSwarms[i].index){
// 			good=true;
// 			break;
// 		}
// 	}
// 	
// 	if(good){
// 		//re-learn and re-sample when errors in the covariance matrix are detected in the sampling phase
// 		bool success=mySampleReplacement(sw, 10);
// 		while(!success){
// 			myLearn(sw);
// 			success=mySampleReplacement(sw, 10);
// 			fprintf(stderr,"Resample\n");
// 		}
// 	}else{
// 		for(int i=0;i<sw.getSize();i++)
// 			sw.particles[i].solution.evalSolution=false;
// // 		//re-learn and re-sample when errors in the covariance matrix are detected in the sampling phase
// // 		bool success=mySample(sw);
// // 		while(!success){
// // 			myLearn(sw);
// // 			success=mySample(sw);
// // 			fprintf(stderr,"Resample\n");
// // 		}
// 	}
	
		//************************************* END OF RESAMPLING SOLUTIONS FROM THE GOOD FRONTS ******************************//
	
	if(logCMAES){
		printf("\n\npop after \n\n");
		for(int i=0; i<sw.getSize();i++){
			//printVector(evo.rgrgx[i],N+2);
			printVector(sw.particles[i].solution.decisionVector,decisionNumber);
			printf("\n");
		}
	}
}

double cmaes_random_Gauss(){
	double x1, x2, rquad, fac;
	
	do {
		x1 = 2.0 * rand()/(double)RAND_MAX - 1.0;
		x2 = 2.0 * rand()/(double)RAND_MAX - 1.0;
		rquad = x1*x1 + x2*x2;
	} while(rquad >= 1 || rquad <= 0);
	fac = sqrt(-2.0*log(rquad)/rquad);
	return fac * x2;
}
static double myhypot(double a, double b) {
/* sqrt(a^2 + b^2) numerically stable. */
	double r = 0;
	if (fabs(a) > fabs(b)) {
		r = b/a;
		r = fabs(a)*sqrt(1+r*r);
	} else if (b != 0) {
		r = a/b;
		r = fabs(b)*sqrt(1+r*r);
	}
	return r;
}

static void QLalgo2 (int n, double *d, double *e, double** V) {
	/*
	 * - > n     : Dimen*sion. 
	 * -> d     : Diagonale of tridiagonal matrix. 
	 * -> e[1..n-1] : off-diagonal, output from Householder
	 * -> V     : matrix output von Householder
	 * <- d     : eigenvalues
	 * <- e     : garbage?
	 * <- V     : basis of eigenvectors, according to d
	 * 
	 * Symmetric tridiagonal QL algorithm, iterative 
	 * Computes the eigensystem from a tridiagonal matrix in roughtly 3N^3 operations
	 * 
	 * code adapted from Java JAMA package, function tql2. 
	 */
	
	int i, k, l, m;
	double f = 0.0;
	double tst1 = 0.0;
	double eps = 2.22e-16; /* Math.pow(2.0,-52.0);  == 2.22e-16 */
	
	/* shift input e */
	for (i = 1; i < n; i++) {
		e[i-1] = e[i];
	}
	e[n-1] = 0.0; /* never changed again */
	
	for (l = 0; l < n; l++) { 
		
		/* Find small subdiagonal element */
		
		if (tst1 < fabs(d[l]) + fabs(e[l]))
			tst1 = fabs(d[l]) + fabs(e[l]);
		m = l;
		while (m < n) {
			if (fabs(e[m]) <= eps*tst1) {
				/* if (fabs(e[m]) + fabs(d[m]+d[m+1]) == fabs(d[m]+d[m+1])) { */
				break;
			}
			m++;
		}
		
		/* If m == l, d[l] is an eigenvalue, */
		/* otherwise, iterate. */
		
		if (m > l) {  /* TODO: check the case m == n, should be rejected here!? */
			int iter = 0;
			do { /* while (fabs(e[l]) > eps*tst1); */
				double dl1, h;
				double g = d[l];
				double p = (d[l+1] - g) / (2.0 * e[l]); 
				double r = myhypot(p, 1.); 
				
				iter = iter + 1;  /* Could check iteration count here */
				
				/* Compute implicit shift */
				
				if (p < 0) {
					r = -r;
				}
				d[l] = e[l] / (p + r);
				d[l+1] = e[l] * (p + r);
				dl1 = d[l+1];
				h = g - d[l];
				for (i = l+2; i < n; i++) {
					d[i] -= h;
				}
				f = f + h;
				
				/* Implicit QL transformation. */
				
				p = d[m];
				{
					double c = 1.0;
					double c2 = c;
					double c3 = c;
					double el1 = e[l+1];
					double s = 0.0;
					double s2 = 0.0;
					for (i = m-1; i >= l; i--) {
						c3 = c2;
						c2 = c;
						s2 = s;
						g = c * e[i];
						h = c * p;
						r = myhypot(p, e[i]);
						e[i+1] = s * r;
						s = e[i] / r;
						c = p / r;
						p = c * d[i] - s * g;
						d[i+1] = h + s * (c * g + s * d[i]);
						
						/* Accumulate transformation. */
						
						for (k = 0; k < n; k++) {
							h = V[k][i+1];
							V[k][i+1] = s * V[k][i] + c * h;
							V[k][i] = c * V[k][i] - s * h;
						}
					}
					p = -s * s2 * c3 * el1 * e[l] / dl1;
					e[l] = s * p;
					d[l] = c * p;
				}
				
				/* Check for convergence. */
				
			} while (fabs(e[l]) > eps*tst1);
		}
		d[l] = d[l] + f;
		e[l] = 0.0;
	}
	
	/* Sort eigenvalues and corresponding vectors. */
	#if 1
	/* TODO: really needed here? So far not, but practical and only O(n^2) */
	{
		int j; 
		double p;
		for (i = 0; i < n-1; i++) {
			k = i;
			p = d[i];
			for (j = i+1; j < n; j++) {
				if (d[j] < p) {
					k = j;
					p = d[j];
				}
			}
			if (k != i) {
				d[k] = d[i];
				d[i] = p;
				for (j = 0; j < n; j++) {
					p = V[j][i];
					V[j][i] = V[j][k];
					V[j][k] = p;
				}
			}
		}
	}
	#endif 
} /* QLalgo2 */ 


/* ========================================================= */
static void Householder2(int n, double** V, double *d, double *e) {
	/* 
	 * H ouseholder tra*nsformation of a symmetric matrix V into tridiagonal form. 
	 * -> n             : dimension
	 * -> V             : symmetric nxn-matrix
	 * <- V             : orthogonal transformation matrix:
	 * tridiag matrix == V * V_in * V^t
	 * <- d             : diagonal
	 * <- e[0..n-1]     : off diagonal (elements 1..n-1) 
	 * 
	 * code slightly adapted from the Java JAMA package, function private tred2()  
	 * 
	 */
	
	int i,j,k; 
	
	for (j = 0; j < n; j++) {
		d[j] = V[n-1][j];
	}
	
	/* Householder reduction to tridiagonal form */
	
	for (i = n-1; i > 0; i--) {
		
		/* Scale to avoid under/overflow */
		
		double scale = 0.0;
		double h = 0.0;
		for (k = 0; k < i; k++) {
			scale = scale + fabs(d[k]);
		}
		if (scale == 0.0) {
			e[i] = d[i-1];
			for (j = 0; j < i; j++) {
				d[j] = V[i-1][j];
				V[i][j] = 0.0;
				V[j][i] = 0.0;
			}
		} else {
			
			/* Generate Householder vector */
			
			double f, g, hh;
			
			for (k = 0; k < i; k++) {
				d[k] /= scale;
				h += d[k] * d[k];
			}
			f = d[i-1];
			g = sqrt(h);
			if (f > 0) {
				g = -g;
			}
			e[i] = scale * g;
			h = h - f * g;
			d[i-1] = f - g;
			for (j = 0; j < i; j++) {
				e[j] = 0.0;
			}
			
			/* Apply similarity transformation to remaining columns */
			
			for (j = 0; j < i; j++) {
				f = d[j];
				V[j][i] = f;
				g = e[j] + V[j][j] * f;
				for (k = j+1; k <= i-1; k++) {
					g += V[k][j] * d[k];
					e[k] += V[k][j] * f;
				}
				e[j] = g;
			}
			f = 0.0;
			for (j = 0; j < i; j++) {
				e[j] /= h;
				f += e[j] * d[j];
			}
			hh = f / (h + h);
			for (j = 0; j < i; j++) {
				e[j] -= hh * d[j];
			}
			for (j = 0; j < i; j++) {
				f = d[j];
				g = e[j];
				for (k = j; k <= i-1; k++) {
					V[k][j] -= (f * e[k] + g * d[k]);
				}
				d[j] = V[i-1][j];
				V[i][j] = 0.0;
			}
		}
		d[i] = h;
	}
	
	/* Accumulate transformations */
	
	for (i = 0; i < n-1; i++) {
		double h; 
		V[n-1][i] = V[i][i];
		V[i][i] = 1.0;
		h = d[i+1];
		if (h != 0.0) {
			for (k = 0; k <= i; k++) {
				d[k] = V[k][i+1] / h;
			}
			for (j = 0; j <= i; j++) {
				double g = 0.0;
				for (k = 0; k <= i; k++) {
					g += V[k][i+1] * V[k][j];
				}
				for (k = 0; k <= i; k++) {
					V[k][j] -= g * d[k];
				}
			}
		}
		for (k = 0; k <= i; k++) {
			V[k][i+1] = 0.0;
		}
	}
	for (j = 0; j < n; j++) {
		d[j] = V[n-1][j];
		V[n-1][j] = 0.0;
	}
	V[n-1][n-1] = 1.0;
	e[0] = 0.0;
	
} /* Housholder() */

/* ========================================================= */
static void Eigen( int N,  double **C, double *diag, double **Q, double *rgtmp)
/* 
 * Calculating eigenvalues and vectors. 
 * Input: 
 *  N: dimension.
 *  C: symmetric (1:N)xN-matrix, solely used to copy data to Q
 *  niter: number of maximal iterations for QL-Algorithm. 
 *  rgtmp: N+1-dimensional vector for temporal use. 
 * Output: 
 *  diag: N eigenvalues. 
 *  Q: Columns are normalized eigenvectors.
 */
{
	int i, j;
	
	if (rgtmp == NULL){ /* was OK in former versions */
		//FATAL("cmaes_t:Eigen(): input parameter double *rgtmp must be non-NULL", 0,0,0);
		fprintf(stderr,"cmaes_t:Eigen(): input parameter double *rgtmp must be non-NULL");
		exit(1);
	}
	
	/* copy C to Q */
	if (C != Q) {
		for (i=0; i < N; ++i)
			for (j = 0; j <= i; ++j)
				Q[i][j] = Q[j][i] = C[i][j];
	}
	
	#if 0
	Householder( N, Q, diag, rgtmp);
	QLalgo( N, diag, Q, 30*N, rgtmp+1);
	#else
	Householder2( N, Q, diag, rgtmp);
	QLalgo2( N, diag, rgtmp, Q);
	#endif
}

static double rgdouMax( const double *rgd, int len){
	int i;
	double max = rgd[0];
	for (i = 1; i < len; ++i)
		max = (max < rgd[i]) ? rgd[i] : max;
	return max;
}

static double rgdouMin( const double *rgd, int len){
	int i;
	double min = rgd[0];
	for (i = 1; i < len; ++i)
		min = (min > rgd[i]) ? rgd[i] : min;
	return min;
}