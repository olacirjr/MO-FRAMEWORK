int numMixComponents=3; //parametro do aprendizado
int selsize; // numero de elementos selecionados para treinar o modelo
int stringlength; //numero de variaveis de decisao
int	kappa=5; //number of parents
int	saveModel=0; //salvar ou nao o modelo em disco
int	maximumAmountOfClusters=200;  // Maximum number of possible clusters for leader algorithm.
double	leaderThreshold=0.3; //randomized leader algorithm (RLA) with a threshold value of 0.3
int useClusters=1;
int	generation; //contador pra geracao atual
int	popsize; //tamanho da populacao
//double 	**population;
Repository *rep_eda;
Swarm* swarm_eda;
bool init=true;

//#define GETREAL(INDEX,MEMBER) (population[(MEMBER)][(INDEX)])  // obtem um membro da populacao
//#define GETREALSELECTED(INDEX,MEMBER) (selected[(MEMBER)][(INDEX)]) // obtem um membro selecionado da populacao (fronteira?)
//#define SETREAL(REAL,INDEX,MEMBER) (population[(MEMBER)][(INDEX)] = REAL) //insere um membro na populacao
//#define SETREALSELECTED(REAL,INDEX,MEMBER) (selected[(MEMBER)][(INDEX)] = REAL) //insere um membro selecionado na populacao (fronteira?)
#define SMALLER_VALUE_BETTER(X1,X2) (X1 < X2)

double GETREALSELECTED(int j, int i){
	if(j>= stringlength)
		printf("ERRO!!! stringlength");
	if(i>= rep_eda->getActualSize())
		printf("ERRO!!! selsize (%d)\n", i);
	//return selected[i][j];
	//return 0;
// 	return swarm_eda->repository.getSolution(i).decisionVector[j];
	return rep_eda->getSolution(i).decisionVector[j];
}

double GETREAL(int j, int i){
	if(j>= stringlength)
		printf("ERRO!!! stringlength");
	if(i>= popsize)
		printf("ERRO!!! popsize");
	//return selected[i][j];
	//return 0;
	return swarm_eda->particles[i].solution.decisionVector[j];
}
void SETREAL(double real, int j, int i){
	if(j>= stringlength)
		printf("ERRO!!! stringlength");
	if(i>= popsize)
		printf("ERRO!!! popsize");
	//return selected[i][j];
	//return 0;
// 	printf("\na: %f", swarm_eda.particles[i].solution.decisionVector[j]);
// 	printf("\nd: %f", real);
	
	
	swarm_eda->particles[i].solution.decisionVector[j]=real;
}

#include "Miscellaneous.h"
#include "Matrix.h"
#include "RandomNumber.h"
#include "rBOABody.h"

void r_BOA(Swarm &sw, int swarm_edaSize, bool init){
	swarm_eda=&sw;
	Repository repTemp;
	
	if(decomposition){
		mergeNeighboringRepositories(repTemp, sw);
		repTemp.organize();
		
		selsize=repTemp.getActualSize();
		rep_eda=&repTemp;
		
	}else{
		selsize=swarm_eda->repository.getActualSize();
		rep_eda=&sw.repository;
	}
	
	stringlength=decisionNumber;
	popsize=swarm_edaSize;
    
	if(init)
		initRandom();

	clusteringWholeProblem();
	
	modelSelection();
	
	decomposeProblem();
	
	clusteringSubproblems( useClusters, leaderThreshold, maximumAmountOfClusters);
	
	modelFitting();
		
	selsize=0;
	
	generateNewSolutions();
	
	//memcpy(&sw, &swarm_eda, sizeof(Swarm));
	
	//free
	for(int i = 0; i < maxNumClustersWholeProblem; i++ )
		free(wholeClusters[i]);
	free(wholeClusters);
	wholeClusters=0;
	free(wholeClustersSize);

	//free memory of clustering subproblems
	free(clustersUsedOverall);
	free(clustersUsed);
	for(int i = 0; i < stringlength; i++ ){
		for(int j = 0; j < maximumAmountOfClusters; j++ )
			free(clusters[i][j]);

		free(clusters[i]);
		free(clustersSize[i]);
		free(alpha[i]);
	}
	free(clusters);
	clusters=0;
	free(clustersSize);
	free(alpha);
	
	//new free additions to comply with changing the number of decision variables there is a huge memory leak
	for(int  i = 0; i < maxNumClustersWholeProblem; i++ ){
		for(int  j = 0; j < stringlength; j++ )
			free(logarithmicProbNew[i][j]);

		free(logarithmicProbCurrent[i]);
		free(logarithmicProbNew[i]);
	}
	free(logarithmicProbCurrent);
	free(logarithmicProbNew);
	logarithmicProbCurrent=0;
		
	for(int i = 0; i < maxNumClustersWholeProblem; i++ ) {
		free(mu[i]);
		for(int x = 0; x < stringlength; x++ )//matrix free
			free(sigma_check[i][x]);
		free(sigma_check[i]);
	}
	free(mu);
	free(sigma_check);
	mu=0;
		
	for(int i = 0; i < maxNumClustersWholeProblem; i++ ){
		
		for(int j = 0; j < stringlength; j++ ){
			free(d_new_family[i][j]);
			free(d_new_parents[i][j]);
		}
		free(d_current_family[i]);
		free(d_current_parents[i]);
		free(d_new_family[i]);
		free(d_new_parents[i]);
	}
	free(d_current_family);
	free(d_current_parents);
	free(d_new_family);
	free(d_new_parents);
	d_new_family=0;
	
	for(int i = 0; i < stringlength; i++ )
		free(pi[i]);
	free(pi);
	free(piLen);
	free(omega);
	pi=0;
	
	for(int i = 0; i < stringlength; i++ )
		free(SubProblems[i]);
	free(SubProblemsLen);
	free(SubProblems);
	SubProblemsLen=0;

	for(int i = 0; i < stringlength; i++ ) 
		free(omega_sub[i]);
		
	free(omega_sub);
	omega_sub=0;
	
	for(int i = 0; i < stringlength; i++ ) {
		for(int j = 0; j < maximumAmountOfClusters; j++ ) {
			if( (mu_sub[i][j]) ){
				free(mu_sub[i][j]);
				for(int x = 0; x < stringlength; x++ )//matrix free
					free(sigma_sub_check[i][j][x]);
				free(sigma_sub_check[i][j]);
			}
		}
	}
	
	for(int i = 0; i < stringlength; i++ ) {
		for(int j = 0; j < maximumAmountOfClusters; j++ ) {
			if( (sigma_sub_i[i][j]) ){
				for(int x = 0; x < stringlength; x++ )//matrix free
					free(sigma_sub_i[i][j][x]);
				free(sigma_sub_i[i][j]);
				free(sigma_sub_tilde[i][j]);
				free(mu_sub_tilde[i][j]);
			}
		}
	}
	
	for(int i = 0; i < stringlength; i++ ) {
		free(mu_sub[i]);
		free(mu_sub_tilde[i]);
		free(sigma_sub_check[i]);
		free(sigma_sub_i[i]);
		free(sigma_sub_tilde[i]);
	}
	free(mu_sub);
	free(mu_sub_tilde);
	free(sigma_sub_check);
	free(sigma_sub_i);
	free(sigma_sub_tilde);
	mu_sub=0;
	
	free(sci);
	scilen = 0;
	sci=0;
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
// 	for( int i = 0; i < maxNumClustersWholeProblem; i++ ) {
// 		free(mu[i]);
// 		free(sigma_check[i]);
// 	}
// 	free(mu);
// 	mu=0;
// 	free(sigma_check);
	
// 	for(int i = 0; i < maxNumClustersWholeProblem; i++ ) {
// 		for(int j = 0; j < stringlength; j++ ){
// 			free(d_new_family[i][j]);
// 			free(d_new_parents[i][j]);
// 		}
// 		free(d_current_family[i]);
// 		free(d_current_parents[i]);
// 		free(d_new_family[i]);
// 		free(d_new_parents[i]);
// 	}
// 	free(d_current_family);
// 	free(d_current_parents);
// 	free(d_new_family);
// 	d_new_family=0;
// 	free(d_new_parents);
// 	
// 	for(int i = 0; i < maxNumClustersWholeProblem; i++ ){
// 		for(int j = 0; j < stringlength; j++ )
// 			free(logarithmicProbNew[i][j]);
// 		
// 		free(logarithmicProbCurrent[i]);
// 		free(logarithmicProbNew[i]);
// 	}
// 	free(logarithmicProbCurrent);
// 	free(logarithmicProbNew);
// 	logarithmicProbCurrent=0;
// 	
// 
// 	for(int i = 0; i < stringlength; i++ )
// 		free(SubProblems[i]);
// 	
// 	free(SubProblemsLen);
// 	free(SubProblems);
// 	SubProblemsLen=0;
// 	
// 	for(int i = 0; i < stringlength; i++ ) 
// 		free(omega_sub[i]);
// 		
// 	free(omega_sub);
// 	omega_sub=0;
// 
// 		
// 	for(int i = 0; i < stringlength; i++ ) {
// 		for(int j = 0; j < maximumAmountOfClusters; j++ ){
// 			if(mu_sub[i][j]){
// 				free(mu_sub[i][j]);
// 				free(sigma_sub_check[i][j]);
// 			}
// 		}
// 		free(mu_sub[i]);
// 		free(mu_sub_tilde[i]);
// 		free(sigma_sub_check[i]);
// 		free(sigma_sub_i[i]);
// 		free(sigma_sub_tilde[i]);
// 	}
// 	free(mu_sub);
// 	free(mu_sub_tilde);
// 	free(sigma_sub_check);
// 	free(sigma_sub_i);
// 	free(sigma_sub_tilde);
// 	mu_sub=0;
// 	
// 	for(int i = 0; i < stringlength; i++ )
// 		free(pi[i]);
// 	free(pi);
// 	free(piLen);
// 	free(omega);
// 	pi=0;
	
	
// 	printf("\n particle after: ");
// 	printVector(swarm_eda.particles[10].solution.objectiveVector, objectiveNumber);
// 	printf(" -> ");
// 	printVector(swarm_eda.particles[10].solution.decisionVector, decisionNumber);
}