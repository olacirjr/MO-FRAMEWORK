/*#############################################################
  ##                                                         ##
  ##                      RandomNumber.h                     ##
  ##           -------------------------------------         ##
  ##             Copyright (c) 2004 Chang Wook Ahn           ##  
  ##     Gwangju Institute of Science & Technology (GIST)    ##
  ##                                                         ##
  ##             Main Work: Random number generator          ##
  ##                                                         ##
  ##  It generates a random number with a seed which is set  ##
  ##    by the system time, so no need to manually input it. ##
  ##                                                         ##
  #############################################################*/


#include <assert.h>

int RandomSeed;

/****************************************************************
**       random computes a pseudo-random                       **
**       double value between 0 and 1, excluding 1.            **
*****************************************************************/

#define RANDOM01() ((double)uniform01())

#define NEXTRANDOMNUMBER() ((int) (RANDOM01()*(RAND_MAX+1)))

#define RANDOM_MAX_VALUE   (RAND_MAX+1)

#define RANDOMNUMBER(MAX)  ((long) (RANDOM01()*((double) MAX)))

#define SETSEED(SEED) {RandomSeed = (SEED); srand( (unsigned int) (SEED) );}

void initRandom( void );
double uniform01( void );

long M = 2147483647;
long A = 16807;
long Q = ((long) M/A);
long R = ((long) M%A);


/****** Initialize the random number seed by means of the current CPU clock *******/
void initRandom( void ) {   
    int seed;
    struct timeb t;
    ftime(&t);

    double longSeed = double(t.millitm)/1000.0; //For Random Seed among (0 ~ 1)
    seed = int( 1 + longSeed * (M-2) );
    SETSEED(seed);

    assert( seed>0 && seed<M );
}

double uniform01() {
    long lo, hi, test;

    hi = RandomSeed/Q;
    lo = RandomSeed%Q;
    test = A*lo - R*hi;
    if( test>0 )
        RandomSeed = test;
    else
        RandomSeed = test+M;
    return double(RandomSeed)/double(M);
}
