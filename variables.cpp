//prototipation of some classes to be declared in this file
class Repository;
class Swarm;
class Solution;

//prototipation of functions to be used throughout the algorithm, since this is the first file included in main
int dominance(const double *sol1, const double *sol2, const int &objectiveNumber);
int compute_b_mga(Solution* solutions, int repActualSize);
bool isEqual(double* vec1, double* vec2, int size);
bool vectorZero(double* x, int size);
void IMultiOperations(int iteration);
void hmopsoBefore();
void hmopsoAfter();
void mostrador(int iteration);
void algorithmInitializationOperations(const int argc, const char* argv[]);
void readParameters(const int argc, const char *inputFile);
void finalizeRun(timeval startTime);
void calculateSigmaVector(double* out, double* in);
void normalizeObjectives(double* out, const double* in);
void printVector(double* vec, int vecSize);
void obtainIdealSolution(Solution* solutions, Solution ideal, int objectiveNumber, int repActualSize);
void box_mga(Solution &solution, double* output, int b, int objectiveNumber);
void updateCrowdingDistances(Solution* solutions, int solutionNumber);
void updateContributingR2(Repository &rep);
void updateContributingHypervolume(Repository &rep);
void updateLargestSmallestDecVectors(Swarm &swarm);
void updateMaxMinObjs();
double round_5(const double &in);
double getScalarValue(double* sol, double* weight);
double PBI(double* sol, double* weight);
double TCH(double* sol, double* weight);
double combination(const int &m, const int &n);
double calculateNWSum(Solution &solution, Solution &GBest);
double calculateEuclideanDistance(double* vector1, double* vector2, int &size);
int localOptimizerObjFunc(int n, int m, double *x, double *f, double *con, void *state_);
bool vectorNan(double* x, int &size);
double averageWeightedLogLikelihood();
double averageScalarObjFromPopulation();
double averageScalarObjFromRepository();
void printVectorToFile(double* vec, int vecSize, char* fileName);
void mergeNeighboringRepositories(Repository &repOut, Swarm &sw);
void clearBetweenRuns();
double cmaes_random_Gauss();

//defining PSO general parameters
const double MAX_PHI=1;
// const double MAX_OMEGA=0.8;
const double INERTIA=0.1;
double *deltaMax=NULL, *deltaMin=NULL; //JMetal - used in the velocity constriction
double *inferiorPositionLimit=NULL, *superiorPositionLimit=NULL; //position limits of the decision variables
// const double VELOCITY_REDUCTION=0.001;
// const double MAX_MUT=0.5;
// const double superiorVelocityLimit=5; //Andre
// const double inferiorVelocityLimit=-5; //Andre
// double mutationIndice=0.15; //percentage of the particles to be influenced by the turbulence operator
//end of defining PSO general parameters


//defining general global variables to the entire search
int objectiveNumber=-1, decisionNumber=-1; //number of objectives and decision variables
int originalPopSize=-1;
Swarm *swarms=NULL; //Swarms
char problem[10]; //DTLZ(1-7), WFG(1-9) or tandem(1-24)
char outputFileName[300]; //name of the output file
int repositorySize=-1; //maximum number of solutions that the repository can hold
int swarmNumber=-1; //the number of swarms used in the search
int maxIterations=-1; //number of iterations of of the search
int maxEvals=-1;
int numExec=-1; //number of independent runs
char archiver[50]; //random (rnd), crowding distance (cd), ideal, mga
char leader[50]; //random (rnd), crowding distance (cd), NWSum, Sigma
char algorithm[100]; //which algorithm is to be executed
// char line[10000000];
double *minObjectives=NULL, *maxObjectives=NULL; //maximum and minimum objective values found so far in the search (used for normalization)
struct timeval startTime, endTime, tmpTime;// structures to store the running time
int funcEvals=0;

bool showTimes=false; //if enabled, show the time elapsed to do each operation
bool printHV=false; //if enabled, print the hypervolume per iteration
bool printIGD=false; //if enabled, print the IGD per iteration
bool printMOEAD_CMAES=false;//if enabled, print the average weighted log likelihood and the average scalarized value for all subproblems (only for decomposition)
//end of defining global variables


//defining names of files that can be used on output
char _s[300]; //name of file to store the optimal solutions
char _f[300]; //name of file to store the optimal front
char _p[300]; //name of file to store the probabilities at each iteration
char _hv[300]; //name of file to store the hypervolume
char _igd[300]; //name of file to store the IGD
char _MOEAD_CMAES[300]; //name of the file to store the MOEA/D-CMAES performance data
// char _partEnter[300]; //name of file to store the optimal solutions
// char _trunc[300]; //name of file to store the optimal solutions
// char _matCov[300]; //name of the file to store the covariance matrices
// char _solObj[300]; //name of the file to store the solutions/objectives as Roberto requested
//end of defining names of files


//multi-swarm parameters
int subSwarmTotalPopulation=-1; //total population in the multi-swarm phase, this will be divided by the swarmNumber - Imulti original = 750
int diversityIterations=-1; //number of diversity iterations before the multi-swarm phase begins - Imulti original = 100
int numParticoes=-1; //number of partitioning iterations throughout the multi-swarm phase - IMulti original = 5
char archSubSwarms[50]; //archiver to be used in the multi-swarm phase. Normal archivers: random (rnd), crowding distance (cd), ideal, mga. Decomposition archivers: pbi, tch, w-ideal, r-ideal
double rangeInicial=0.5; // initial value of the range to truncate the solutions around the centroid - IMulti original = 0.5
double rangeFinal=0.1; // final value of the the range to truncate the solutions around the centroid - IMulti original = 0.1
double increment=0; //increment value to go from the rangeInicial to the rangeFinal, set automatically during the run
double range=rangeInicial;// current range used to to truncate the solutions around the centroid, this value is automatically updated during the run
int subSwarmPopulation=-1;//number of solutions per sub swarm. This will be automatically updated during the run;
int originalSwarmNumber=-1;//stores the original swarm number, the number of swarms can change during the search if there is not enough non-dominated solutions
char clusteringType[50], clusteringDistance[50]; // space used for applying the clustering method (objectives, decision, both) and clustering distance used (euclidean, thebycheff or minkowski(0.5,4) )
char truncType[50];//type of truncation applied when a solution is out of allowed bounds (random, rdmInSwarm, imulti and extremes)
char weightDistribution[50]; //equal, log, linear (cma-es)
char solSet[50]; //repository, population, both (cma-es)
char cmaRanking[50]; //cd, hv, r2 (cma-es)
bool diversityPhase=false; //flag indicating if this is one of the diversisy iterations on multi-swarm algorithms, set automatically during the run
bool decomposition=false; //flag indicating if decomposition approach is being used, controls several steps in the algorithm where the use of decomposition makes sense. set automatically during the run
// end of multi-swarm parameters


//H-MOPSO, parameters
// Repository repPrevious;
const int combinations=9;//number of combinations between leader selection and archiving methods
double probabilities[combinations];//stores the current probabilities of the roulette
double **reference=NULL;//store the reference points(weights) to be used in the r2 calculation
char file[50];//name of the file containing the reference points to calculate r2
int refSize=-1;//number of reference points used to calculate r2
int used=-1;//index of the low-level heuristic used in the current iteration
//end of H-MOPSO, parameters

// Defining parameters for the trajectory problems
int seq[5];// sequence of planets for the tandem problem, each sequence defines one problem = {3, 2, 3, 3, 6}; // Earth, Venus, Earth,Earth, Staturn (EVEES)
const double pi_const = acos(-1.0);//used in tandem
double bestTandem=MAXDOUBLE;//best value found so far for the first objective
// int b=0; w=0; e=0;
// end of parameters for the trajectory problems

//decomposition parameters
int maxReplacements=2; //maximum number of neighbors a solution can update - original 2
int neighborhoodSize=-1; // size of the neighborhood - 0 means auto
Repository *repGlobal=NULL, *repTemp=NULL; //repGlobal, a global external repository to store all the nondominated solutions found so far. RepTemp, a temporary repository to be used in different situations
double* globalTmpWeight=NULL;// used to transfer the weight vector for the local optimizer for tandem calculation, set automatically
int updateNeighborhoodMetric=0;
double delta=0.9;//parameter that sets the probability of learning and updating just from the neighborhood, otherwise learn from all the subproblems
//end of decomposition parameters