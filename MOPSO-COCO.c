/**
 USING THE MOPSO ALGORITHM IN THE COCO FRAMEWORK
 *
 * Set the global parameter BUDGET_MULTIPLIER to suit your needs.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "MOPSO.cpp"

#include "additionalCode/coco/code-experiments/build/c/coco.h"
#include "additionalCode/coco/code-experiments/build/c/coco.c"

/**
 * The maximal budget for evaluations done by an optimization algorithm equals dimension * BUDGET_MULTIPLIER.
 * Increase the budget multiplier value gradually to see how it affects the runtime.
 */
static const size_t BUDGET_MULTIPLIER = 2;

/**
 * The maximal number of independent restarts allowed for an algorithm that restarts itself.
 */
static const size_t INDEPENDENT_RESTARTS = 1e5;

/**
 * The random seed. Change if needed.
 */
static const uint32_t RANDOM_SEED = 0xdeadbeef;

/**
 * A function type for evaluation functions, where the first argument is the vector to be evaluated and the
 * second argument the vector to which the evaluation result is stored.
 */
typedef void (*evaluate_function_t)(const double *x, double *y);

/**
 * A pointer to the problem to be optimized (needed in order to simplify the interface between the optimization
 * algorithm and the COCO platform).
 */
static coco_problem_t *PROBLEM;

/**
 * The function that calls the evaluation of the first vector on the problem to be optimized and stores the
 * evaluation result in the second vector.
 */
static void evaluate_function(const double *x, double *y) {
  coco_evaluate_function(PROBLEM, x, y);
}

/* Declarations of all functions implemented in this file (so that their order is not important): */
void example_experiment(const char *suite_name, const char *observer_name);

/**
 * The main method initializes the random number generator and calls the example experiment on the
 * bi-objective suite.
 */
int main(const int argc, const char* argv[]) {
	coco_random_state_t *random_generator = coco_random_new(RANDOM_SEED);

	/* Change the log level to "warning" to get less output */
	coco_set_log_level("warning");

	printf("Running the example experiment... (might take time, be patient)\n");
	fflush(stdout);
	
	//INITIALIZE THE MOPSO FRAMEWORK
	
	readParameters(argc, argv[1]);
	algorithmInitializationOperations(argc, argv);

	example_experiment("bbob-biobj", "bbob-biobj");

	/* Uncomment the line below to run the same example experiment on the bbob suite
	example_experiment("bbob", "bbob", random_generator); */

	printf("Done!\n");
	fflush(stdout);

	coco_random_free(random_generator);

	return 0;
}

void Problem::calculateCoco(double* decisionVector, double* objectiveVector){
	coco_evaluate_function(PROBLEM, decisionVector, objectiveVector);
}

/**
 * A simple example of benchmarking random search on a suite with instances from 2016.
 *
 * @param suite_name Name of the suite (use "bbob" for the single-objective and "bbob-biobj" for the
 * bi-objective suite).
 * @param observer_name Name of the observer (use "bbob" for the single-objective and "bbob-biobj" for the
 * bi-objective observer).
 * @param random_generator The random number generator.
 */
void example_experiment(const char *suite_name, const char *observer_name) {
	size_t r;
	coco_suite_t *suite;
	coco_observer_t *observer;
	
	//********************** DETERMINE THE PROBLEM*************************//
	int problemId=-1;
	if(!strncmp(problem, "coco", 4)){//if the problem is from the coco framework
		if(strlen(problem) > 4 && strlen(problem) < 7){
			char tmp[3]="";
			strncpy(tmp, problem+4, strlen(problem)-4);
			problemId=atoi(tmp)-1;
			if(problemId < 0 || problemId > 54){
				printf("\nERROR! invalid problem ID for coco, must be in 1 - 55, you set %d\n", problemId+1);
				exit(1);
			}
		}else{
			printf("\nERROR! on problem ID of coco (%s) (%d)\n", problem, problemId);
			exit(1);
		}
	}else{
		printf("\nERROR! non coco problem in the coco framework\n");
		exit(1);
	}
	//************************END OF DETERMINING THE PROBLEM********************//

	/* Set some options for the observer. See documentation for other options. */
	char *observer_options =
		coco_strdupf("result_folder: %s-%s "
				"algorithm_name: %s "
				"algorithm_info: \"%s Framework\"", outputFileName, problem, algorithm, algorithm);

	/* Initialize the suite and observer */
	suite = coco_suite(suite_name, "year: 2016", "dimensions: 2,3,5,10,20,40"); //original
// 	suite = coco_suite(suite_name, "year: 2016", "dimensions: 2,3,5,10,20");
// 	suite = coco_suite(suite_name, "year: 2016", "dimensions: 40");
// 	suite = coco_suite(suite_name, "year: 2016", "dimensions: 20");
// 	suite = coco_suite(suite_name, "year: 2016", "dimensions: 5");
// 	suite = coco_suite(suite_name, "year: 2016", "dimensions: 2");
// 	suite = coco_suite(suite_name, "year: 2016", "dimensions: 2,3,5");

	observer = coco_observer(observer_name, observer_options);
	coco_free_memory(observer_options);
	
	for(int i=0;i<6;i++){
// 	for(int i=0;i<5;i++){
// 	for(int i=0;i<1;i++){
// 	int i=5;{
		for(int j=0;j<10;j++){
			if (suite->current_problem) {
				coco_problem_free(suite->current_problem);
			}
			PROBLEM = coco_suite_get_problem_from_indices(suite, problemId, i, j);//coco_suite_get_problem_from_indices(suite, function_idx, dimension_idx, instance_idx);
			PROBLEM = coco_problem_add_observer(PROBLEM, observer);
			suite->current_problem = PROBLEM;
		
		/* Iterate over all problems in the suite */
	// 	while ((PROBLEM = coco_suite_get_next_problem(suite, observer)) != NULL) {
			size_t dimension = coco_problem_get_dimension(PROBLEM);
			size_t evaluations_done = coco_problem_get_evaluations(PROBLEM);
			long evaluations_remaining = (long) (dimension * BUDGET_MULTIPLIER) - (long) evaluations_done;

			//SET THE PROBLEM VARIABLES AND PARAMETERS
			sprintf(problem, "coco%d",(int)PROBLEM->suite_dep_function);//force set the problem
			decisionNumber=(int)dimension;//set the decision number
			objectiveNumber=(int)coco_problem_get_number_of_objectives(PROBLEM);//set the objective number - for now only bi-objective problems
		// 	maxEvals=(int)evaluations_remaining;//if set stop by max evals, otherwise stop by iteration number
			if(inferiorPositionLimit != NULL)//unset previous inferior position limits if any
				delete[] inferiorPositionLimit;
			if(superiorPositionLimit != NULL)//unset previous superior position limits if any
				delete[] superiorPositionLimit;
			
			inferiorPositionLimit = new double[decisionNumber];
			superiorPositionLimit = new double[decisionNumber];
			for (int var = 0; var < decisionNumber; var++) {
				inferiorPositionLimit[var]=coco_problem_get_smallest_values_of_interest(PROBLEM)[var];
				superiorPositionLimit[var]=coco_problem_get_largest_values_of_interest(PROBLEM)[var];
			}
					
			//FINALLY CALL THE MOPSO FRAMEWORK
			run();
			printf("\n");
			
		// 	printf("%d == %d -- %d\n", funcEvals, (int)coco_problem_get_evaluations(PROBLEM), (int)evaluations_remaining);


			/* Break the loop if the algorithm performed no evaluations or an unexpected thing happened */
			if (coco_problem_get_evaluations(PROBLEM) == evaluations_done) {
			printf("WARNING: Budget has not been exhausted (%lu/%lu evaluations done)!\n", evaluations_done,
				dimension * BUDGET_MULTIPLIER);
			break;
			}
			else if (coco_problem_get_evaluations(PROBLEM) < evaluations_done)
			coco_error("Something unexpected happened - function evaluations were decreased!");
		//     }
		}
	}
	coco_observer_free(observer);
	coco_suite_free(suite);
}

/**
 * A random search algorithm that can be used for single- as well as multi-objective optimization.
 *
 * @param evaluate The evaluation function used to evaluate the solutions.
 * @param dimension The number of variables.
 * @param number_of_objectives The number of objectives.
 * @param lower_bounds The lower bounds of the region of interested (a vector containing dimension values).
 * @param upper_bounds The upper bounds of the region of interested (a vector containing dimension values).
 * @param max_budget The maximal number of evaluations.
 * @param random_generator Pointer to a random number generator able to produce uniformly and normally
 * distributed random numbers.
 */
void my_random_search(evaluate_function_t evaluate, const size_t dimension, const size_t number_of_objectives, const double *lower_bounds, const double *upper_bounds, const size_t max_budget, coco_random_state_t *random_generator) {
	double *x = coco_allocate_vector(dimension);
	double *y = coco_allocate_vector(number_of_objectives);
	double range;
	size_t i, j;

	for (i = 0; i < max_budget; ++i) {
		/* Construct x as a random point between the lower and upper bounds */
		for (j = 0; j < dimension; ++j) {
			range = upper_bounds[j] - lower_bounds[j];
			x[j] = lower_bounds[j] + coco_random_uniform(random_generator) * range;
		}

		/* Call the evaluate function to evaluate x on the current problem (this is where all the COCO logging
			* is performed) */
		evaluate(x, y);
	}
	coco_free_memory(x);
	coco_free_memory(y);
}