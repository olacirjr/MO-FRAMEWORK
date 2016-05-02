#include "additionalCode/problem/trajobjfuns.h"
#include "additionalCode/problem/trajobjfuns.cpp"
#include "additionalCode/problem/Astro_Functions.h"
#include "additionalCode/problem/Astro_Functions.cpp"
#include "additionalCode/problem/mga.h"
#include "additionalCode/problem/mga.cpp"
#include "additionalCode/problem/time2distance.h"
#include "additionalCode/problem/time2distance.cpp"
#include "additionalCode/problem/Lambert.h"
#include "additionalCode/problem/Lambert.cpp"
#include "additionalCode/problem/mga_dsm.h"
#include "additionalCode/problem/mga_dsm.cpp"
#include "additionalCode/problem/Pl_Eph_An.h"
#include "additionalCode/problem/Pl_Eph_An.cpp"
#include "additionalCode/problem/propagateKEP.h"
#include "additionalCode/problem/propagateKEP.cpp"
#include "additionalCode/problem/misc4Tandem.h"
#include "additionalCode/problem/misc4Tandem.cpp"
#include "additionalCode/problem/PowSwingByInv.h"
#include "additionalCode/problem/PowSwingByInv.cpp"
#include "additionalCode/problem/ZeroFinder.h"
#include "additionalCode/problem/ZeroFinder.cpp"

class Problem{
	public:
		//swich between the problems implemented and evaluate the solution with the chosen
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void evaluate(double* decisionVector, double* objectiveVector);
		
	private:
	// 	Problem(const Problem &source){}//copy
		Problem& operator= (const Problem &source){}//assignment
		
		//function used in the WFG problems
		//params - the decision vector to be evaluated
		void calculate_x(double* t);
		//evaluate a solution using the DTLZ1 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateDTLZ1(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the DTLZ2 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateDTLZ2(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the DTLZ3 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateDTLZ3(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the DTLZ4 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateDTLZ4(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the DTLZ5 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateDTLZ5(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the DTLZ6 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateDTLZ6(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the DTLZ7 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateDTLZ7(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the WFG1 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateWFG1(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the WFG2 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateWFG2(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the WFG3 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateWFG3(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the WFG4 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateWFG4(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the WFG5 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateWFG5(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the WFG6 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateWFG6(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the WFG7 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateWFG7(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the WFG8 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateWFG8(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the WFG9 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateWFG9(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the UF1 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateUF1(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the UF2 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateUF2(double* decisionVector, double* objectiveVector);
				//evaluate a solution using the UF3 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateUF3(double* decisionVector, double* objectiveVector);
				//evaluate a solution using the UF4 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateUF4(double* decisionVector, double* objectiveVector);
				//evaluate a solution using the UF5 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateUF5(double* decisionVector, double* objectiveVector);
				//evaluate a solution using the UF6 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateUF6(double* decisionVector, double* objectiveVector);
				//evaluate a solution using the UF7 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateUF7(double* decisionVector, double* objectiveVector);
				//evaluate a solution using the UF8 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateUF8(double* decisionVector, double* objectiveVector);
				//evaluate a solution using the UF9 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateUF9(double* decisionVector, double* objectiveVector);
				//evaluate a solution using the UF10 problem
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateUF10(double* decisionVector, double* objectiveVector);
		//evaluate a solution using the Tandem function as implemented in  trajobjfuns.cpp
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateTandem(double* decisionVector, double* objectiveVector);
		//evaluate a solution using one of the functions of the coco framework
		//params - the decision vector to be evaluated
		//		 - the objective vector to hold the evaluation values
		void calculateCoco(double* decisionVector, double* objectiveVector);
};


void printVector(double* vec, int vecSize);

//swich between the problems implemented and evaluate the solution with the chosen
//params - the decision vector to be evaluated
//		 - the objective vector to hold the evaluation values
void Problem::evaluate(double* decisionVector, double* objectiveVector){
	if(maxEvals > 0 && maxEvals <= funcEvals){
		for(int i=0;i<objectiveNumber;i++)
			objectiveVector[i]=MAXDOUBLE;
		return;
	}
	
	if(!strcmp(problem, "dtlz1"))
		calculateDTLZ1(decisionVector, objectiveVector);
	if(!strcmp(problem, "dtlz2"))
		calculateDTLZ2(decisionVector, objectiveVector);
	if(!strcmp(problem, "dtlz3"))
		calculateDTLZ3(decisionVector, objectiveVector);
	if(!strcmp(problem, "dtlz4"))
		calculateDTLZ4(decisionVector, objectiveVector);
	if(!strcmp(problem, "dtlz5"))
		calculateDTLZ5(decisionVector, objectiveVector);
	if(!strcmp(problem, "dtlz6"))
		calculateDTLZ6(decisionVector, objectiveVector);
	if(!strcmp(problem, "dtlz7"))
		calculateDTLZ7(decisionVector, objectiveVector);
	if(!strcmp(problem, "wfg1"))
		calculateWFG1(decisionVector, objectiveVector);
	if(!strcmp(problem, "wfg2"))
		calculateWFG2(decisionVector, objectiveVector);
	if(!strcmp(problem, "wfg3"))
		calculateWFG3(decisionVector, objectiveVector);
	if(!strcmp(problem, "wfg4"))
		calculateWFG4(decisionVector, objectiveVector);
	if(!strcmp(problem, "wfg5"))
		calculateWFG5(decisionVector, objectiveVector);
	if(!strcmp(problem, "wfg6"))
		calculateWFG6(decisionVector, objectiveVector);
	if(!strcmp(problem, "wfg7"))
		calculateWFG7(decisionVector, objectiveVector);
	if(!strcmp(problem, "wfg8"))
		calculateWFG8(decisionVector, objectiveVector);
	if(!strcmp(problem, "wfg9"))
		calculateWFG9(decisionVector, objectiveVector);
	if(!strcmp(problem, "uf1"))
		calculateUF1(decisionVector, objectiveVector);
	if(!strcmp(problem, "uf2"))
		calculateUF2(decisionVector, objectiveVector);
	if(!strcmp(problem, "uf3"))
		calculateUF3(decisionVector, objectiveVector);
	if(!strcmp(problem, "uf4"))
		calculateUF4(decisionVector, objectiveVector);
	if(!strcmp(problem, "uf5"))
		calculateUF5(decisionVector, objectiveVector);
	if(!strcmp(problem, "uf6"))
		calculateUF6(decisionVector, objectiveVector);
	if(!strcmp(problem, "uf7"))
		calculateUF7(decisionVector, objectiveVector);
	if(!strcmp(problem, "uf8"))
		calculateUF8(decisionVector, objectiveVector);
	if(!strcmp(problem, "uf9"))
		calculateUF9(decisionVector, objectiveVector);
	if(!strcmp(problem, "uf10"))
		calculateUF10(decisionVector, objectiveVector);
	if(!strncmp(problem, "tandem",6))
		calculateTandem(decisionVector, objectiveVector);
	if(!strncmp(problem, "coco",4))
		calculateCoco(decisionVector, objectiveVector);
	
	if(strcmp(problem, "dtlz1") && strcmp(problem, "dtlz2") && strcmp(problem, "dtlz3") && strcmp(problem, "dtlz4") && 
	strcmp(problem, "dtlz5") && strcmp(problem, "dtlz6") && strcmp(problem, "dtlz7") && strcmp(problem, "wfg1") && 
	strcmp(problem, "wfg2") && strcmp(problem, "wfg3") && strcmp(problem, "wfg4") && strcmp(problem, "wfg5") &&
	strcmp(problem, "wfg6") && strcmp(problem, "wfg7") && strcmp(problem, "wfg8") && strcmp(problem, "wfg9") && 
	strcmp(problem, "uf1") && strcmp(problem, "uf2") && strcmp(problem, "uf3") && strcmp(problem, "uf4") && strcmp(problem, "uf5") &&
	strcmp(problem, "uf6") && strcmp(problem, "uf7") && strcmp(problem, "uf8") && strcmp(problem, "uf9") && strcmp(problem, "uf10") && strncmp(problem, "tandem",6) && strncmp(problem, "coco",4))
		printf("INVALID OPTIMIZATION PROBLEM! (%s)\n", problem);
	funcEvals++;
}
//evaluate a solution using the DTLZ1 problem
//params - the decision vector to be evaluated
//		 - the objective vector to hold the evaluation values
void Problem::calculateDTLZ1(double* decisionVector, double* objectiveVector){

	int k= decisionNumber - objectiveNumber + 1;
	double g=0.0;

	for(int i=decisionNumber-k;i<decisionNumber;i++)
		g+=(decisionVector[i]-0.5)*(decisionVector[i]-0.5)-cos(20.0*M_PI*(decisionVector[i]-0.5));

	g=100*(k+g);
	for(int i=0;i<objectiveNumber;i++)
		objectiveVector[i] = (1.0+g)*0.5;

	for(int i=0;i<objectiveNumber;i++){
		for(int j=0;j<objectiveNumber-(i+1);j++)
			objectiveVector[i] *= decisionVector[j];
			if(i != 0){
				int aux = objectiveNumber - (i+1);
				objectiveVector[i] *= 1-decisionVector[aux] ;
			}
	}
}
//evaluate a solution using the DTLZ2 problem
//params - the decision vector to be evaluated
//		 - the objective vector to hold the evaluation values
void Problem::calculateDTLZ2(double* decisionVector, double* objectiveVector){

	int k= decisionNumber - objectiveNumber + 1;
	double g=0.0;

	for(int i=decisionNumber-k;i<decisionNumber;i++)
		g+=(decisionVector[i]-0.5)*(decisionVector[i]-0.5);

	for(int i=0;i<objectiveNumber;i++)
		objectiveVector[i] = (1.0+g);

	for(int i=0;i<objectiveNumber;i++){
		for(int j=0;j<objectiveNumber-(i+1);j++)
			objectiveVector[i] *= cos(decisionVector[j]*0.5*M_PI);
			if(i != 0){
				int aux = objectiveNumber - (i+1);
				objectiveVector[i] *= sin(decisionVector[aux]*0.5*M_PI);
			}
	}
}
//evaluate a solution using the DTLZ3 problem
//params - the decision vector to be evaluated
//		 - the objective vector to hold the evaluation values
void Problem::calculateDTLZ3(double* decisionVector, double* objectiveVector){

	int k= decisionNumber - objectiveNumber + 1;
	double g=0.0;

	for(int i=decisionNumber-k;i<decisionNumber;i++)
		g+=(decisionVector[i]-0.5)*(decisionVector[i]-0.5)-cos(20.0*M_PI*(decisionVector[i]-0.5));

	g=100.0*(k+g);
	for(int i=0;i<objectiveNumber;i++)
		objectiveVector[i] = 1.0+g;

	for(int i=0;i<objectiveNumber;i++){
		for(int j=0;j<objectiveNumber-(i+1);j++)
			objectiveVector[i] *= cos(decisionVector[j]*0.5*M_PI);
			if(i != 0){
				int aux = objectiveNumber - (i+1);
				objectiveVector[i] *= sin(decisionVector[aux]*0.5*M_PI);
			}
	}
}
//evaluate a solution using the DTLZ4 problem
//params - the decision vector to be evaluated
//		 - the objective vector to hold the evaluation values
void Problem::calculateDTLZ4(double* decisionVector, double* objectiveVector){

	int k= decisionNumber - objectiveNumber + 1;
	double g=0.0;
	double alpha=100.0;

	for(int i=decisionNumber-k;i<decisionNumber;i++)
		g+=(decisionVector[i]-0.5)*(decisionVector[i]-0.5);

	for(int i=0;i<objectiveNumber;i++)
		objectiveVector[i] = (1.0+g);

	for(int i=0;i<objectiveNumber;i++){
		for(int j=0;j<objectiveNumber-(i+1);j++)
			objectiveVector[i] *= cos(pow(decisionVector[j],alpha)*(M_PI/2));
			if(i != 0){
				int aux = objectiveNumber - (i+1);
				objectiveVector[i] *= sin(pow(decisionVector[aux],alpha)*(M_PI/2));
			}
	}
}
//evaluate a solution using the DTLZ5 problem
//params - the decision vector to be evaluated
//		 - the objective vector to hold the evaluation values
void Problem::calculateDTLZ5(double* decisionVector, double* objectiveVector){

	int k= decisionNumber - objectiveNumber + 1;
	double g=0.0;
	double theta[objectiveNumber-1];

	for(int i=decisionNumber-k;i<decisionNumber;i++)
		g+=(decisionVector[i]-0.5)*(decisionVector[i]-0.5);

	double t = M_PI/(4.0*(1.0+g));

	theta[0]=decisionVector[0]*(M_PI/2.0);
	for(int i=1;i<(objectiveNumber-1);i++)
		theta[i] = t*(1.0+2.0*g*decisionVector[i]);

	for(int i=0;i<objectiveNumber;i++)
		objectiveVector[i] = (1.0+g);

	for(int i=0;i<objectiveNumber;i++){
		for(int j=0;j<objectiveNumber-(i+1);j++)
			objectiveVector[i] *= cos(theta[j]);
			if(i != 0){
				int aux = objectiveNumber - (i+1);
				objectiveVector[i] *= sin(theta[aux]) ;
			}
	}
}
//evaluate a solution using the DTLZ6 problem
//params - the decision vector to be evaluated
//		 - the objective vector to hold the evaluation values
void Problem::calculateDTLZ6(double* decisionVector, double* objectiveVector){

	int k= decisionNumber - objectiveNumber + 1;
	double g=0.0;
	double theta[objectiveNumber-1];

	for(int i=decisionNumber-k;i<decisionNumber;i++)
		g+=pow(decisionVector[i],0.1);

	double t = M_PI/(4.0*(1.0+g));
	theta[0]=decisionVector[0]*(M_PI/2.0);
	for(int i=1;i<(objectiveNumber-1);i++)
		theta[i] = t*(1.0+2.0*g*decisionVector[i]);

	for(int i=0;i<objectiveNumber;i++)
		objectiveVector[i] = (1.0+g);

	for(int i=0;i<objectiveNumber;i++){
		for(int j=0;j<objectiveNumber-(i+1);j++)
			objectiveVector[i] *= cos(theta[j]);
			if(i != 0){
				int aux = objectiveNumber - (i+1);
				objectiveVector[i] *= sin(theta[aux]) ;
			}
	}
}
//evaluate a solution using the DTLZ7 problem
//params - the decision vector to be evaluated
//		 - the objective vector to hold the evaluation values
void Problem::calculateDTLZ7(double* decisionVector, double* objectiveVector){
	int k= decisionNumber - objectiveNumber + 1;
	double g=0.0;

	for(int i=decisionNumber-k;i<decisionNumber;i++)
		g+=decisionVector[i];

	g=1+(9.0*g)/k;

	for(int i=0;i<objectiveNumber-1;i++)
		objectiveVector[i] = decisionVector[i];

	double h=0.0;
	for(int i=0;i<objectiveNumber-1;i++)
		h+=(objectiveVector[i]/(1.0+g))*(1+sin(3.0*M_PI*objectiveVector[i]));

	h=objectiveNumber-h;

	objectiveVector[objectiveNumber-1]=(1+g)*h;
}

/**************************************************************************************  WFG FUNCTIONS ************************************************************************/

double correct_to_01(double a){
	double min = (double)0.0;
	double max = (double)1.0;
	double epsilon = (double)1e-7;

	double min_epsilon = min - epsilon;
	double max_epsilon = max + epsilon;

	if (( a <= min && a >= min_epsilon ) || (a >= min && a <= min_epsilon)) {
		return min;        
	} else if (( a >= max && a <= max_epsilon ) || (a <= max && a >= max_epsilon)) {
		return max;        
	} else {
		return a;        
	}
}
void normalize(double* decisionVector, double* result){
	for (int i = 0; i < decisionNumber; i++){
		double bound = (double)2.0 * (i + 1);
		result[i] = decisionVector[i] / bound;
		result[i]=correct_to_01(result[i]);
	}
}    
double s_linear(double y, double A){
	return correct_to_01(fabs(y - A) /(double)fabs(floor(A - y) + A));
}

double s_multi(double y, int A, int B, double C){                
	double tmp1, tmp2;
		
	tmp1 = ((double)4.0 * A + (double)2.0) * (double)M_PI * ((double)0.5 - fabs(y - C) / ((double)2.0 * ((double)floor(C - y) + C)));
	tmp2 = (double)4.0 * B * (double)pow(fabs(y - C) / ((double)2.0 * ((double)floor(C - y) + C)), (double)2.0);
		
	return correct_to_01(((double)1.0 + (double)cos(tmp1) + tmp2) / (B + (double)2.0));
}

double s_decept(double y, double A, double B, double C){        
	double tmp, tmp1, tmp2;
		
	tmp1 = (double)floor(y - A + B) * ((double)1.0 - C + (A - B)/B) / (A - B);
	tmp2 = (double)floor(A + B - y) * ((double)1.0 - C + ((double)1.0 - A - B) / B) / ((double)1.0 - A - B);
		
	tmp = fabs(y - A) - B;
		
	return correct_to_01((double)1 + tmp * (tmp1 + tmp2 + (double)1.0/B));
}

double b_flat(double y, double A, double B, double C){    
	double tmp1 = std::min((double)0, (double)floor(y - B))* A*(B-y)/B;
	double tmp2 = std::min((double)0, (double)floor(C - y))* (1 - A)*(y - C)/(1 - C);
		
	return correct_to_01(A + tmp1 - tmp2);
}
double b_poly(double y, double alpha){
	if (! ( y>=0 || y<=1 || alpha>0 || alpha != 1 ) )
		printf("ERROR ON WFG FUNCTION! (b_poly)");
	
	return correct_to_01((double)pow(y,alpha));
}
double b_param(double y, double u, double A, double B, double C){
	double result, v, exp;
		
	v = A - ((double)1.0 - (double)2.0 * u) * fabs((double)floor((double)0.5 - u) + A);
	exp = B + (C - B)*v;
	result = (double)pow(y,exp);
		
	return correct_to_01(result);                  
}
void subVector(double* z, int head, int tail, double* result){	
	for( int i = head; i <= tail; i++ ){
		result[i-head]=z[i];
	}
}

double r_sum(double* y, double* w, int length){
	double tmp1 = (double)0.0, tmp2 =(double) 0.0;
	for (int i = 0; i < length; i++){
		tmp1 += y[i]*w[i];
		tmp2 += w[i];
	}
		
	return correct_to_01(tmp1 / tmp2);
}
double r_nonsep(double* y, int A, int size){
	double tmp, denominator, numerator;
		
	tmp = (double)ceil(A/(double)2.0);        
	denominator = size * tmp * ((double)1.0 + (double)2.0*A - (double)2.0*tmp)/A;        
	numerator = (double)0.0;
	for (int j = 0; j < size; j++){
		numerator += y[j];
		for (int k = 0; k <= A-2; k++){
		numerator += fabs( y[j] - y[( j+k+1 ) % size]);
		}
	}
		
	return correct_to_01(numerator/denominator);
}

double convex(double* x, int m){
	double result = (double)1.0;
		
	for (int i = 1; i <= objectiveNumber - m; i++)
		result *= (1 - cos(x[i-1] * M_PI * 0.5));
					
	if (m != 1)
		result *= (1 - sin(x[objectiveNumber - m] * M_PI * 0.5));
		
	return result;
}
double concave(double* x, int m){
	double result = (double)1.0;
		
	for (int i = 1; i <= objectiveNumber - m; i++)
		result *= sin(x[i-1] * M_PI * 0.5);
		
	if (m != 1)
		result *= cos(x[objectiveNumber - m] * M_PI * 0.5);
			
	return result;
}
double linear(double* x, int m){
	double  result = (double)1.0;        

	for (int i = 1; i <= objectiveNumber - m; i++)
		result *= x[i-1];
		
	if (m != 1)        
		result *= (1 - x[objectiveNumber - m]);
				
	return result;
}

double mixed(double* x, int A, double alpha){
	double tmp;        
	tmp =(double) cos((double)2.0 * A * (double)M_PI * x[0] + (double)M_PI * (double)0.5);
	tmp /= (2.0 * (double) A * M_PI);
	return (double)pow(((double)1.0 - x[0] - tmp),alpha);
}

double disc(double* x, int A, double alpha, double beta){
	double tmp;
	tmp = (double)cos((double)A * pow(x[0], beta) * M_PI);
		
	return (double)1.0 - (double)pow(x[0],alpha) * (double)pow(tmp,2.0);
}

void Problem::calculate_x(double* t){
	double tmp[objectiveNumber];
	double A[objectiveNumber-1];
	
	if(!strcmp(problem, "wfg3")){
		A[0] = 1;
		for (int i = 1; i < objectiveNumber-1; i++) 
			A[i] = 0;
	}else
		for (int i = 0; i < objectiveNumber-1; i++)
			A[i] = 1;
	
	
		
	for (int i = 0; i < objectiveNumber-1; i++){
		tmp[i] = std::max(t[objectiveNumber-1],A[i]) * (t[i]  - (double)0.5) + (double)0.5;
	}
		
	tmp[objectiveNumber-1] = t[objectiveNumber-1];
		
	memcpy(t, &tmp, sizeof(double)*objectiveNumber);
}

void WFG1_t1(double* z, int k){
	double tmp[decisionNumber];
	
	memcpy(&tmp, z, sizeof(double)*decisionNumber);
		
	for (int i = k; i < decisionNumber; i++) {
		tmp[i]=s_linear(z[i],0.35);
	}
	
	memcpy(z, &tmp, sizeof(double)*decisionNumber);
}

void WFG4_t1(double* z){
	double tmp[decisionNumber];
	
	for (int i = 0; i < decisionNumber; i++) {
		tmp[i] = s_multi(z[i],30,10,(double)0.35);
	}
	memcpy(z, &tmp, sizeof(double)*decisionNumber);
}

void WFG5_t1(double* z){
	double tmp[decisionNumber];
	
	for (int i = 0; i < decisionNumber; i++) {
		tmp[i] = s_decept(z[i],(double)0.35,(double)0.001,(double)0.05);
	}
	memcpy(z, &tmp, sizeof(double)*decisionNumber);
}

void WFG7_t1(double* z, int k){
	double tmp[decisionNumber];
	double w[decisionNumber];
	double subZ[decisionNumber];
	double subW[decisionNumber];
	
	memcpy(&tmp, z, sizeof(double)*decisionNumber);

	for (int i = 0; i < decisionNumber; i++) {
		w[i] = (double)1.0;
	}

	for (int i = 0; i < k; i++){
		int head = i+1;
		int tail = decisionNumber-1;
		subVector(z,head,tail,subZ);
		subVector(w,head,tail,subW);
		double aux = r_sum(subZ,subW,(tail-head+1));

		tmp[i] = b_param(z[i],aux,(double)0.98/(double)49.98,(double)0.02,(double)50);
	}
	memcpy(z, &tmp, sizeof(double)*decisionNumber);
}

void WFG8_t1(double* z, int k){
	double tmp[decisionNumber];
	double w[decisionNumber];
	double subZ[decisionNumber];
	double subW[decisionNumber];
		
	for (int i = 0; i < decisionNumber; i++) {
		w[i] = (double)1.0;
	}

	memcpy(&tmp, z, sizeof(double)*decisionNumber);
		
	for (int i = k; i < decisionNumber; i++){
		int head = 0;
		int tail = i - 1;
		subVector(z,head,tail,subZ);
		subVector(w,head,tail,subW);            
		double aux = r_sum(subZ,subW,(tail-head+1));
			
		tmp[i] = b_param(z[i],aux,(double)0.98/(double)49.98,(double)0.02,50);
	}
		
	memcpy(z, &tmp, sizeof(double)*decisionNumber);
}

void WFG9_t1(double* z, int k){
	double tmp[decisionNumber];
	double w[decisionNumber];
	double subZ[decisionNumber];
	double subW[decisionNumber];
		
	for (int i = 0; i < decisionNumber; i++) {
		w[i] = (double)1.0;
	}

	memcpy(&tmp, z, sizeof(double)*decisionNumber);
		
	for (int i = 0; i < decisionNumber-1; i++){
      int head = i+1;
      int tail = decisionNumber-1;
      subVector(z,head,tail, subZ);
      subVector(w,head,tail,subW);
      double aux = r_sum(subZ,subW,(tail-head+1));
      tmp[i] = b_param(z[i],aux,(double)0.98/(double)49.98,(double)0.02,(double)50);
    }
        
    tmp[decisionNumber-1] = z[decisionNumber-1];
		
	memcpy(z, &tmp, sizeof(double)*decisionNumber);
}

void WFG1_t2(double* z, int k){
	double tmp[decisionNumber];

	memcpy(&tmp, z, sizeof(double)*decisionNumber);
		
	for (int i = k; i < decisionNumber; i++) {
		tmp[i]=b_flat(z[i],(double)0.8,(double)0.75,(double)0.85);
	}
		
	memcpy(z, &tmp, sizeof(double)*decisionNumber);
}

void WFG2_t2(double* z, int k){
	double tmp[decisionNumber];
	double subZ[decisionNumber];
	memcpy(&tmp, z, sizeof(double)*decisionNumber);
	
	int l = decisionNumber - k;
		
	for (int i = k+1; i <= k + l/2; i++){
		int head = k + 2*(i - k) - 1;
		int tail = k + 2*(i - k);              
		subVector(z,head-1,tail-1,subZ);
		
		tmp[i-1] = r_nonsep(subZ,2,(tail-head+1) );
	}
	memcpy(z, &tmp, sizeof(double)*decisionNumber);
}

void WFG4_t2(double* z, int k){
	double tmp[objectiveNumber];
	double w[decisionNumber];
	double subZ[decisionNumber];
	double subW[decisionNumber];

	for (int i = 0; i < decisionNumber; i++) {
		w[i] = (double)1.0;
	}
			
	for (int i = 1; i <= objectiveNumber-1; i++) {
		int head = (i - 1)*k/(objectiveNumber-1) + 1;
		int tail = i * k / (objectiveNumber - 1);
		subVector(z,head-1,tail-1,subZ);
		subVector(w,head-1,tail-1,subW);
		
		tmp[i-1] = r_sum(subZ,subW,(tail-head+1));
	}
		
	int head = k + 1;
	int tail = decisionNumber;
		
	subVector(z,head-1,tail-1,subZ);
	subVector(w,head-1,tail-1,subW);
	tmp[objectiveNumber-1] = r_sum(subZ,subW,(tail-head+1));
		
	memcpy(z, &tmp, sizeof(double)*objectiveNumber);
}

void WFG6_t2(double* z, int k){
	double tmp[objectiveNumber];
	double subZ[decisionNumber];
		
	for (int i = 1; i <= objectiveNumber-1; i++){
		int head = (i - 1)*k/(objectiveNumber-1) + 1;
		int tail = i * k / (objectiveNumber - 1);
		subVector(z,head-1,tail-1, subZ);            
			
		tmp[i-1] = r_nonsep(subZ,k/(objectiveNumber-1),(tail-head+1));            
	}
		
	int head = k + 1;
	int tail = decisionNumber;
	int l = decisionNumber - k;
			
	subVector(z,head-1,tail-1, subZ);              
	tmp[objectiveNumber-1] = r_nonsep(subZ,l,(tail-head+1));
				
	memcpy(z, &tmp, sizeof(double)*objectiveNumber);
}

void WFG9_t2(double* z, int k){
	double tmp[decisionNumber];
		
	for (int i = 0; i < k; i++) {
		tmp[i] = s_decept(z[i],(double)0.35,(double)0.001,(double)0.05);
	}
		
	for (int i = k; i < decisionNumber; i++) {
		tmp[i] = s_multi(z[i],30,95,(double)0.35);
	}        

	memcpy(z, &tmp, sizeof(double)*decisionNumber);
}

void WFG1_t3(double* z){
	double tmp[decisionNumber];
		
	for (int i = 0; i < decisionNumber; i++) {
		tmp[i]=b_poly(z[i],(double)0.02);
	}
		
	memcpy(z, &tmp, sizeof(double)*decisionNumber);
}

void WFG2_t3(double* z, int k){
	double tmp[objectiveNumber];
	double w[decisionNumber];
	
	double subZ[decisionNumber];
	double subW[decisionNumber];

	for (int i = 0; i < decisionNumber; i++) {
		w[i] = (double)1.0;
	}
		
	for (int i = 1; i <= objectiveNumber-1; i++){
		int head = (i - 1)*k/(objectiveNumber-1) + 1;
		int tail = i * k / (objectiveNumber - 1);
		subVector(z,head-1,tail-1,subZ);
		subVector(w,head-1,tail-1,subW);
		
		tmp[i-1] = r_sum(subZ,subW,(tail-head+1));
	}
		
	int l = decisionNumber - k;
	int head = k + 1;
	int tail = k + l / 2;
	subVector(z,head-1,tail-1, subZ);
	subVector(w,head-1,tail-1, subW);

	tmp[objectiveNumber-1] = r_sum(subZ,subW,(tail-head+1));
				
	memcpy(z, &tmp, sizeof(double)*objectiveNumber);
}

void WFG7_t3(double* z, int k){
	double tmp[objectiveNumber];
	double w[decisionNumber];
	double subZ[decisionNumber];
	double subW[decisionNumber];
		
	for (int i = 0; i < decisionNumber; i++) {
		w[i] = (double)1.0;
	}
		
	for (int i = 1; i <= objectiveNumber-1; i++){
		int head = (i - 1)*k/(objectiveNumber-1) + 1;
		int tail = i * k / (objectiveNumber - 1);
		subVector(z,head-1,tail-1, subZ);
		subVector(w,head-1,tail-1, subW);
			
		tmp[i-1] = r_sum(subZ,subW,(tail-head+1));
	}
		
	int head = k + 1;
	int tail = decisionNumber;
	subVector(z,head-1,tail-1, subZ);
	subVector(w,head-1,tail-1, subW);
	tmp[objectiveNumber-1] = r_sum(subZ,subW,(tail-head+1));

	memcpy(z, &tmp, sizeof(double)*objectiveNumber);
}

void WFG1_t4(double* z, int k){
	double tmp[objectiveNumber];
	double w[decisionNumber];
	double subZ[decisionNumber];
	double subW[decisionNumber];
				
	for (int i = 0; i < decisionNumber; i++) {
		w[i] = (double)2.0 * (i + 1);
	}
		
	for (int i = 1; i <= objectiveNumber-1; i++){
		int head = (i - 1)*k/(objectiveNumber-1) + 1;
		int tail = i * k / (objectiveNumber - 1);
		
		subVector(z,head-1,tail-1,subZ);
		subVector(w,head-1,tail-1,subW);		
		
		tmp[i-1]=r_sum(subZ,subW,(tail-head+1));
	}
		
	int head = k + 1 - 1;
	int tail = decisionNumber - 1;    
	
	subVector(z,head,tail,subZ);
	subVector(w,head,tail,subW);
	
	tmp[objectiveNumber-1]=r_sum(subZ,subW,(tail-head+1));
	
	memcpy(z, &tmp, sizeof(double)*objectiveNumber);
}
// bool vectorIn01(double* x){
// 	for( int i = 0; i < decisionNumber; i++ ){
// 		if( x[i] < 0.0 || x[i] > 1.0 ){
// 			return false;
// 		}
// 	}
// 
//   return true;
// }

void Problem::calculateWFG1(double* decisionVector, double* objectiveVector){
	double y[decisionNumber];
	int k=4;
	if(objectiveNumber > 2)
		k=2*(objectiveNumber-1);

	normalize(decisionVector, y);	
	WFG1_t1(y,k);
	WFG1_t2(y,k);	
	WFG1_t3(y);
	WFG1_t4(y,k);
	
	calculate_x(y);
	
	for (int m = 1; m <= objectiveNumber - 1 ; m++) {
		objectiveVector[m-1] = y[objectiveNumber-1] + superiorPositionLimit[m-1] * convex(y,m);
	}
	objectiveVector[objectiveNumber-1] = y[objectiveNumber-1] + superiorPositionLimit[objectiveNumber-1] * mixed(y,5,(double)1.0);
}

void Problem::calculateWFG2(double* decisionVector, double* objectiveVector){
	double y[decisionNumber];
	int k=4;
	if(objectiveNumber > 2)
		k=2*(objectiveNumber-1);
	
	normalize(decisionVector, y);
	WFG1_t1(y,k); // = wfg2_t1
	WFG2_t2(y,k);
	WFG2_t3(y,k);
		
	for (int m = 1; m <= objectiveNumber - 1 ; m++) {
		objectiveVector[m-1] = y[objectiveNumber-1] + superiorPositionLimit[m-1] * convex(y,m);
	}        
	objectiveVector[objectiveNumber-1] = y[objectiveNumber-1] + superiorPositionLimit[objectiveNumber-1] * disc(y,5,(double)1.0,(double)1.0);
}

void Problem::calculateWFG3(double* decisionVector, double* objectiveVector){
	double y[decisionNumber];
	int k=4;
	if(objectiveNumber > 2)
		k=2*(objectiveNumber-1);	
	
	normalize(decisionVector, y);
	WFG1_t1(y,k); //=wfg3_t1
	WFG2_t2(y,k); //=wfg3_t2
	WFG2_t3(y,k); //=wfg3_t3

	calculate_x(y);        
	for (int m = 1; m <= objectiveNumber ; m++) {
		objectiveVector [m-1] = y[objectiveNumber-1] + superiorPositionLimit[m-1] * linear(y,m);
	}
}
void Problem::calculateWFG4(double* decisionVector, double* objectiveVector){
	double y[decisionNumber];
	int k=4;
	if(objectiveNumber > 2)
		k=2*(objectiveNumber-1);
	
	normalize(decisionVector, y);
	WFG4_t1(y);
	WFG4_t2(y,k);

	calculate_x(y);
	for (int m = 1; m <= objectiveNumber ; m++) {
		objectiveVector [m-1] = y[objectiveNumber-1] + superiorPositionLimit[m-1] * concave(y,m);
	}
}

void Problem::calculateWFG5(double* decisionVector, double* objectiveVector){
	double y[decisionNumber];
	int k=4;
	if(objectiveNumber > 2)
		k=2*(objectiveNumber-1);
	
	normalize(decisionVector, y);
	WFG5_t1(y);
	WFG4_t2(y,k); //=wfg5_t2
		
	calculate_x(y);
	for (int m = 1; m <= objectiveNumber ; m++) {
		objectiveVector [m-1] = y[objectiveNumber-1] + superiorPositionLimit[m-1] * concave(y,m);
	}
}
void Problem::calculateWFG6(double* decisionVector, double* objectiveVector){
	double y[decisionNumber];
	int k=4;
	if(objectiveNumber > 2)
		k=2*(objectiveNumber-1);
	
	normalize(decisionVector, y);
	WFG1_t1(y,k); //=wfg6_t1
	WFG6_t2(y,k);

	calculate_x(y);
	for (int m = 1; m <= objectiveNumber ; m++) {
		objectiveVector[m-1] = y[objectiveNumber-1] + superiorPositionLimit[m-1] * concave(y,m);
	}
}

void Problem::calculateWFG7(double* decisionVector, double* objectiveVector){
	double y[decisionNumber];
	int k=4;
	if(objectiveNumber > 2)
		k=2*(objectiveNumber-1);
	
	normalize(decisionVector, y);
	
	WFG7_t1(y,k);
	WFG1_t1(y,k); //=wfg7_t2
	WFG7_t3(y,k);
		
	calculate_x(y);
	
	for (int m = 1; m <= objectiveNumber ; m++) {
		objectiveVector[m-1] = y[objectiveNumber-1] + superiorPositionLimit[m-1] * concave(y,m);
	}
}
void Problem::calculateWFG8(double* decisionVector, double* objectiveVector){
	double y[decisionNumber];
	int k=4;
	if(objectiveNumber > 2)
		k=2*(objectiveNumber-1);
	
	normalize(decisionVector, y);
	WFG8_t1(y,k);
	WFG1_t1(y,k); //=wfg8_t2
	WFG4_t2(y,k); //=wfg8_t3
		
	calculate_x(y);        
	for (int m = 1; m <= objectiveNumber ; m++) {
		objectiveVector [m-1] = y[objectiveNumber-1] + superiorPositionLimit[m-1] * concave(y,m);                
	}        
}

void Problem::calculateWFG9(double* decisionVector, double* objectiveVector){
	double y[decisionNumber];
	int k=4;
	if(objectiveNumber > 2)
		k=2*(objectiveNumber-1);
	
	normalize(decisionVector, y);
	WFG9_t1(y,k);
	WFG9_t2(y,k);
	WFG6_t2(y,k); //=wfg9_t3
			
	calculate_x(y);        
	for (int m = 1; m <= objectiveNumber ; m++) {
		objectiveVector [m-1] = y[objectiveNumber-1] + superiorPositionLimit[m-1] * concave(y,m);                
	}        
}

/**************************************************************************************  END OF WFG FUNCTIONS **********************************************************************/
/**************************************************************************************  CEC09 UF FUNCTIONS ************************************************************************/

void Problem::calculateUF1(double* decisionVector, double* objectiveVector) {
	int count1, count2;
	double sum1, sum2, yj;
	sum1   = sum2   = 0.0;
	count1 = count2 = 0;

	for (int j = 2 ; j <= decisionNumber; j++) {
		yj = decisionVector[j-1] - sin(6.0*M_PI*decisionVector[0] + j*M_PI/decisionNumber);
		yj = yj * yj;
		if(j % 2 == 0) {
			sum2 += yj;
			count2++;
		} else {
			sum1 += yj;
			count1++;
		}      
	}

	objectiveVector[0]=decisionVector[0] + 2.0 * sum1 / (double)count1;
	objectiveVector[1]=1.0 - sqrt(decisionVector[0]) + 2.0 * sum2 / (double)count2;
	
// 	solution.setObjective(0, decisionVector[0] + 2.0 * sum1 / (double)count1);
// 	solution.setObjective(1, 1.0 - sqrt(decisionVector[0]) + 2.0 * sum2 / (double)count2);
}

void Problem::calculateUF2(double* decisionVector, double* objectiveVector) {
	int count1, count2;
	double sum1, sum2, yj;
	sum1   = sum2   = 0.0;
	count1 = count2 = 0;

	for (int j = 2 ; j <= decisionNumber; j++) {
		if(j % 2 == 0) {
			yj = decisionVector[j-1] -
				(0.3 * decisionVector[0] * decisionVector[0] * cos(24 * M_PI * decisionVector[0] + 4 * j * M_PI / decisionNumber) + 0.6 * decisionVector[0])*
					sin(6.0 * M_PI* decisionVector[0] + j * M_PI / decisionNumber);
			sum2 += yj*yj;
			count2++;
		} else {

			yj = decisionVector[j-1] -
				(0.3 * decisionVector[0] * decisionVector[0] * cos(24 * M_PI * decisionVector[0] + 4 * j * M_PI / decisionNumber) + 0.6 * decisionVector[0])*
					cos(6.0 * M_PI* decisionVector[0] + j * M_PI / decisionNumber);

			sum1 += yj*yj;
			count1++;
		}
	}
	
	objectiveVector[0]=decisionVector[0] + 2.0 * sum1 / (double)count1;
	objectiveVector[1]=1.0 - sqrt(decisionVector[0]) + 2.0 * sum2 / (double)count2;

// 	solution.setObjective(0, decisionVector[0] + 2.0 * sum1 / (double)count1);
// 	solution.setObjective(1, 1.0 - sqrt(decisionVector[0]) + 2.0 * sum2 / (double)count2);
}

void Problem::calculateUF3(double* decisionVector, double* objectiveVector) {
	int count1, count2;
	double sum1, sum2, prod1, prod2, yj, pj;
	sum1   = sum2   = 0.0;
	count1 = count2 = 0;
	prod1  = prod2  = 1.0;


	for (int j = 2 ; j <= decisionNumber; j++) {
		yj = decisionVector[j-1]-pow(decisionVector[0],0.5*(1.0+3.0*(j-2.0)/(decisionNumber-2.0)));
		pj = cos(20.0*yj*M_PI/sqrt(j));
		if (j % 2 == 0) {
			sum2  += yj*yj;
			prod2 *= pj;
			count2++;
		} else {
			sum1  += yj*yj;
			prod1 *= pj;
			count1++;
		}
	}

	objectiveVector[0]=decisionVector[0] + 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
	objectiveVector[1]=1.0 - sqrt(decisionVector[0]) + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
// 	solution.setObjective(0,  decisionVector[0] + 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1);
// 	solution.setObjective(1, 1.0 - sqrt(decisionVector[0]) + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2);
}

void Problem::calculateUF4(double* decisionVector, double* objectiveVector) {
	int count1, count2;
	double sum1, sum2, yj, hj ;
	sum1   = sum2   = 0.0;
	count1 = count2 = 0;

	for (int j = 2 ; j <= decisionNumber; j++) {
		yj = decisionVector[j-1]-sin(6.0*M_PI*decisionVector[0]+j*M_PI/decisionNumber);
		hj = abs(yj)/(1.0+exp(2.0*abs(yj)));
		if (j % 2 == 0) {
			sum2  += hj;
			count2++;
		} else {
			sum1  += hj;
			count1++;
		}
	}

	objectiveVector[0]=decisionVector[0]	+ 2.0*sum1 / (double)count1;
	objectiveVector[1]=1.0 - decisionVector[0]*decisionVector[0]	+ 2.0*sum2 / (double)count2;
// 	solution.setObjective(0, decisionVector[0]	+ 2.0*sum1 / (double)count1);
// 	solution.setObjective(1, 1.0 - decisionVector[0]*decisionVector[0]	+ 2.0*sum2 / (double)count2);
}

void Problem::calculateUF5(double* decisionVector, double* objectiveVector) {
	int n=10;//from JMetal
	double epsilon=0.1;//from JMetal
	int count1, count2;
	double sum1, sum2, yj, hj ;
	sum1   = sum2   = 0.0;
	count1 = count2 = 0;

	for (int j = 2 ; j <= decisionNumber; j++) {
		yj = decisionVector[j-1]-sin(6.0*M_PI*decisionVector[0]+j*M_PI/decisionNumber);
		hj = 2.0*yj*yj - cos(4.0*M_PI*yj) + 1.0;
		if (j % 2 == 0) {
			sum2  += hj;
			count2++;
		} else {
			sum1  += hj;
			count1++;
		}
	}
	hj = (0.5/n + epsilon)*abs(sin(2.0*n*M_PI*decisionVector[0]));

	objectiveVector[0]=decisionVector[0] + hj + 2.0*sum1 / (double)count1;
	objectiveVector[1]=1.0 - decisionVector[0] + hj + 2.0*sum2 / (double)count2;
// 	solution.setObjective(0, decisionVector[0] + hj + 2.0*sum1 / (double)count1);
// 	solution.setObjective(1, 1.0 - decisionVector[0] + hj + 2.0*sum2 / (double)count2);
}

void Problem::calculateUF6(double* decisionVector, double* objectiveVector) {
	int n=2;//from JMetal
	double epsilon=0.1;//from JMetal
	int count1, count2 ;
	double prod1, prod2 ;
	double sum1, sum2, yj, hj, pj ;
	sum1   = sum2   = 0.0;
	count1 = count2 = 0;
	prod1  = prod2  = 1.0;

	for (int j = 2 ; j <= decisionNumber; j++) {
		yj = decisionVector[j-1]-sin(6.0*M_PI*decisionVector[0]+j*M_PI/decisionNumber);
		pj = cos(20.0*yj*M_PI/sqrt(j));
		if (j % 2 == 0) {
			sum2  += yj*yj;
			prod2 *= pj;
			count2++;
		} else {
			sum1  += yj*yj;
			prod1 *= pj;
			count1++;
		}
	}
	hj = 2.0*(0.5/n + epsilon)*sin(2.0*n*M_PI*decisionVector[0]);
	if (hj < 0.0) 
		hj = 0.0;

	objectiveVector[0]=decisionVector[0] + hj + 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
	objectiveVector[1]=1.0 - decisionVector[0] + hj + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
// 	solution.setObjective(0, decisionVector[0] + hj + 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1);
// 	solution.setObjective(1, 1.0 - decisionVector[0] + hj + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2);
}

void Problem::calculateUF7(double* decisionVector, double* objectiveVector) {
	int count1, count2;
	double sum1, sum2, yj;
	sum1   = sum2   = 0.0;
	count1 = count2 = 0;

	for (int j = 2 ; j <= decisionNumber; j++) {
		yj = decisionVector[j-1] - sin(6.0*M_PI*decisionVector[0]+j*M_PI/decisionNumber);
		if (j % 2 == 0) {
			sum2  += yj*yj;
			count2++;
		} else {
			sum1  += yj*yj;
			count1++;
		}
	}
	yj = pow(decisionVector[0],0.2);

	objectiveVector[0]=yj + 2.0*sum1 / (double)count1;
	objectiveVector[1]=1.0 - yj + 2.0*sum2 / (double)count2;
// 	solution.setObjective(0, yj + 2.0*sum1 / (double)count1);
// 	solution.setObjective(1, 1.0 - yj + 2.0*sum2 / (double)count2);
}

void Problem::calculateUF8(double* decisionVector, double* objectiveVector) {
	int count1, count2, count3;
	double sum1, sum2, sum3, yj;
	sum1   = sum2 = sum3 = 0.0;
	count1 = count2 = count3 = 0;

	for (int j = 3 ; j <= decisionNumber; j++) {
		yj = decisionVector[j-1] - 2.0*decisionVector[1]*sin(2.0*M_PI*decisionVector[0]+j*M_PI/decisionNumber);
		if(j % 3 == 1) {
			sum1  += yj*yj;
			count1++;
		} else if(j % 3 == 2) {
			sum2  += yj*yj;
			count2++;
		} else {
			sum3  += yj*yj;
			count3++;
		}
	}

	objectiveVector[0]=cos(0.5*M_PI*decisionVector[0])*cos(0.5*M_PI*decisionVector[1]) + 2.0*sum1 / (double)count1;
	objectiveVector[1]=cos(0.5*M_PI*decisionVector[0])*sin(0.5*M_PI*decisionVector[1]) + 2.0*sum2 / (double)count2;
	objectiveVector[2]=sin(0.5*M_PI*decisionVector[0])                       + 2.0*sum3 / (double)count3;
// 	solution.setObjective(0, cos(0.5*M_PI*decisionVector[0])*cos(0.5*M_PI*decisionVector[1]) + 2.0*sum1 / (double)count1);
// 	solution.setObjective(1, cos(0.5*M_PI*decisionVector[0])*sin(0.5*M_PI*decisionVector[1]) + 2.0*sum2 / (double)count2);
// 	solution.setObjective(2, sin(0.5*M_PI*decisionVector[0])                       + 2.0*sum3 / (double)count3) ;
}

void Problem::calculateUF9(double* decisionVector, double* objectiveVector) {
	double epsilon=0.1;//from JMetal
	int count1, count2, count3;
	double sum1, sum2, sum3, yj;
	sum1   = sum2 = sum3 = 0.0;
	count1 = count2 = count3 = 0;

	for (int j = 3 ; j <= decisionNumber; j++) {
		yj = decisionVector[j-1] - 2.0*decisionVector[1]*sin(2.0*M_PI*decisionVector[0]+j*M_PI/decisionNumber);
		if(j % 3 == 1) {
			sum1  += yj*yj;
			count1++;
		} else if(j % 3 == 2) {
			sum2  += yj*yj;
			count2++;
		} else {
			sum3  += yj*yj;
			count3++;
		}
	}

	yj = (1.0+epsilon)*(1.0-4.0*(2.0*decisionVector[0]-1.0)*(2.0*decisionVector[0]-1.0));
		if (yj < 0.0) 
		yj = 0.0;
	
	objectiveVector[0]=0.5*(yj + 2*decisionVector[0])*decisionVector[1]		+ 2.0*sum1 / (double)count1;
	objectiveVector[1]=0.5*(yj - 2*decisionVector[0] + 2.0)*decisionVector[1] + 2.0*sum2 / (double)count2;
	objectiveVector[2]=1.0 - decisionVector[1]                   + 2.0*sum3 / (double)count3;
// 	solution.setObjective(0, 0.5*(yj + 2*decisionVector[0])*decisionVector[1]		+ 2.0*sum1 / (double)count1);
// 	solution.setObjective(1, 0.5*(yj - 2*decisionVector[0] + 2.0)*decisionVector[1] + 2.0*sum2 / (double)count2);
// 	solution.setObjective(2, 1.0 - decisionVector[1]                   + 2.0*sum3 / (double)count3) ;
}

void Problem::calculateUF10(double* decisionVector, double* objectiveVector) {
	int count1, count2, count3;
	double sum1, sum2, sum3, yj, hj;
	sum1   = sum2 = sum3 = 0.0;
	count1 = count2 = count3 = 0;

	for (int j = 3 ; j <= decisionNumber; j++) {
		yj = decisionVector[j-1] - 2.0*decisionVector[1]*sin(2.0*M_PI*decisionVector[0]+j*M_PI/decisionNumber);
		hj = 4.0*yj*yj - cos(8.0*M_PI*yj) + 1.0;
		if(j % 3 == 1) {
			sum1  += hj;
			count1++;
		} else if(j % 3 == 2) {
			sum2  += hj;
			count2++;
		} else {
			sum3  += hj;
			count3++;
		}
	}

	objectiveVector[0]=cos(0.5*M_PI*decisionVector[0])*cos(0.5*M_PI*decisionVector[1]) + 2.0*sum1 / (double)count1;
	objectiveVector[1]=cos(0.5*M_PI*decisionVector[0])*sin(0.5*M_PI*decisionVector[1]) + 2.0*sum2 / (double)count2;
	objectiveVector[2]=sin(0.5*M_PI*decisionVector[0])                       + 2.0*sum3 / (double)count3;
// 	solution.setObjective(0, cos(0.5*M_PI*decisionVector[0])*cos(0.5*M_PI*decisionVector[1]) + 2.0*sum1 / (double)count1);
// 	solution.setObjective(1, cos(0.5*M_PI*decisionVector[0])*sin(0.5*M_PI*decisionVector[1]) + 2.0*sum2 / (double)count2);
// 	solution.setObjective(2, sin(0.5*M_PI*decisionVector[0])                       + 2.0*sum3 / (double)count3) ;
}

/**************************************************************************************  END OF CEC09 UF FUNCTIONS *****************************************************************/



// double thefunct(double* x){
// 	double f,obj_f ;
// 	ConvertFromCont(x);
// 	MGA_DSM(t, dim, sequence, c_o, rp, e, Isp, mass, AUdist, DVtotal, DVonboard, obj_f, Delta_V);	
// 	f=obj_f ;
// 	
// 	return (f) ;
// }
// double myvalue(asa_objective *asa){/* evaluate the objective function */
// 	double f, t,objf ;
// 	double* x;
// 	
// 	int i, n ;
// 	x = asa->x ;
// 	n = asa->n ;
// 	
// 	// 	f = thefunct(x);
// 	
// // 	std::vector<double> X;
// 	double tof;
// // 	for (int i = 0; i < n; i++)
// // 		X.push_back(x[i]);
// // 	f=tandem(X,tof,seq); 
// // 	f=2000-tandem(X,tof,seq); 
// 	f = 2000-tandem_wbounds(x,tof,seq); 
// 	
// 	return (f) ;
// }
// void mygrad(asa_objective *asa){/* evaluate the gradient of the objective function */
// 	double tt, *g, *x ;
// 	int i, n ;
// 	x = asa->x ;
// 	g = asa->g ;
// 	n = asa->n ;
// 	for (i = 0; i < n; i++){
// 		tt = i + 1 ;
// 		tt = sqrt (tt) ;
// 		//g[i] = exp(x[i]) -  tt ;
// // 		g[i] = 0; //myrand();
// // 		g[i] = (rand()/(double)RAND_MAX); //myrand();
// 		g[i] = 1; //myrand();
// 
// 	}
// // 	g[randomint(n)]=1;
// // 	g[rand()%n]=1;
// 	return ;
// }

// int calcfc(int n, int m, double *x, double *f, double *con, void *state_){
// // 	double obj_f ;
// 	int i;
// 	
// // 	gettimeofday(&endTime, NULL);
// // // 	if(( (double)endTime.tv_sec - (double)tandemTime.tv_sec) > 300){ //if the local optimization is running for over 5 minutes
// // 	if(( (double)endTime.tv_sec - (double)tandemTime.tv_sec) >= 10){ //if the local optimization is running for over 1s
// // 		fprintf(stderr, "Local optimization stop: over 10s.\n");
// // 		return 1; //ask the cobyla to stop
// // 	}
// // 	printf("%f\n", ((double)endTime.tv_sec - (double)tandemTime.tv_sec) );
// 	
// 	
// 	/* Parameter adjustments */
// 	--con;
// 	--x;
// 	bool violate=false;
// 	for(i=0;i<n;i++){
// 		con[i+1] = x[i+1]- inferiorPositionLimit[i];
// 		con[n+i+1] =  superiorPositionLimit[i]- x[i+1];
// 		
// 		if(x[i] > superiorPositionLimit[i]){
// 			violate=true;
// 			break;
// 		}
// 		if(x[i] < inferiorPositionLimit[i]){
// 			violate=true;
// 			break;
// 		}
// 	}
// 	x++;
// 
// // 	std::vector<double> X;
// 	
//         
// //         i=0;
// // 	while(violate==false && i < n){
// // 		if(x[i] > superiorPositionLimit[i]){
// // // 			printf("\t\t\tVariable %d (%f) is larger than %f\n", i, x[i], superiorPositionLimit[i]);
// // 			// 			x[i]=superiorPositionLimit[i];
// // 			violate=true;
// // 		}
// // 		else if(x[i] < inferiorPositionLimit[i]){
// // // 			printf("\t\t\tVariable %d (%f) is smaller than %f\n", i, x[i], inferiorPositionLimit[i]);
// // 			// 			x[i]=inferiorPositionLimit[i];
// // 			violate=true;
// // 		}
// // 		
// // // 		X.push_back(x[i]);
// // 		i++;
// // 	}
// 	
// 	double tof, obj_f;
// 	if(!violate){
// // 		double obj=tandem_wbounds(X,tof,seq);
// 		double obj = tandem_wbounds(x,tof,seq); 
// 		obj_f=2000-obj;
// // 		obj_f=obj;
// 	}else
// 		obj_f = 4000;
// 	
// // 	printf("\ndec_during: (%f) ", obj_f);
// // 	printVector(x, decisionNumber);
// // 	
// // 	printf("\nmin: (%f) ", obj_f);
// // 	printVector(inferiorPositionLimit, decisionNumber);
// // 	
// // 	printf("\nmax: (%f) ", obj_f);
// // 	printVector(superiorPositionLimit, decisionNumber);
// 
// 	if (isnan(obj_f)) obj_f = 4000;
// 	if (obj_f == 0) obj_f = 4000;
// 	*f=obj_f ;
// 	// cout<<setprecision(10)<<endl<<" eval  = " << obj_f<<endl;
// 	return 0;
// } /* calcfc */
typedef struct {
	int nprob;
} example_state;

void Problem::calculateTandem(double* decisionVector, double* objectiveVector){
	int i;
	double tof,obj;         
	//double sol[] = {8614.67190542681782972068, 3.23881168881835801443, 0.49828711089951693847, 0.55247385652360181396, 2323.11617248155198467430, 1911.28129676196658692788, 2104.75935736787505447865, 2499.99977990785100701032, 0.90606524729349469105, 0.83398203560248784783, 0.77167564570127900048, 0.18281414705776333207, 1.41932189054689850138, 1.46704208497941945843, 1.05000000000010906831, -1.24053864588253404122, -1.62734953864368758758, -1.25424110003836464244};
	
	//**********************************************************************************************************************//
	
// 	asacg_parm cgParm ;
// 	asa_parm asaParm ;
// 	asa_cg_default (&cgParm) ;
// 	asa_default (&asaParm) ;
// 	cgParm.PrintParms = FALSE ;
// 	cgParm.PrintLevel = 0 ;
// 	
// 	asaParm.PrintParms = FALSE ;
// 	asaParm.PrintLevel = 0 ;
// 	asaParm.PrintFinal = FALSE;

// 	std::vector<double> aa;
// 	for (int i = 0; i < decisionNumber; i++){aa.push_back(decisionVector[i]);}
// 	double objB = tandem(aa,tof,seq); 
// 	printf("\n\ndec_bef: (%5.8f)", objB);
// 	printVector(decisionVector, decisionNumber);
// 	
// 	asa_cg(decisionVector, inferiorPositionLimit, superiorPositionLimit, decisionNumber, NULL, &cgParm, &asaParm, 1.e-8, myvalue, mygrad, NULL, NULL) ;
// 	
// 	std::vector<double> ab;
// 	for (int i = 0; i < decisionNumber; i++){ab.push_back(decisionVector[i]);}
// 	double objA = tandem(ab,tof,seq); 
// 	printf("\ndec_aft: (%5.8f)", objA);
// 	printVector(decisionVector, decisionNumber);
// 	if(objB < objA)
// 		printf("\nImproved (%5.8f)", objB - objA);
// 	if(objB > objA)
// 		printf("\nWorsen (%5.8f)",objB - objA);
// 	if(objB == objA)
// 		printf("\nEqual (%5.8f)",objB - objA);
// 	
// 	printf("\n\n");
	
	//**********************************************************************************************************************//
	//rc = cobyla(vars, 2*vars, TheBestsol, rhobeg, rhoend, iprint, &maxfun, calcfc, &state);
// 	double  rhobeg = 0.5, rhoend = 0.0001;
// 	int iprint=0, maxfun = 1000;
// 	example_state state;
// 	state.nprob = 0;
// 		
// 	gettimeofday(&tandemTime, NULL);
// 	obj = tandem_wbounds(decisionVector,tof,seq);
// 	objectiveVector[0]=2000-obj;
// 	objectiveVector[1] = tof;
// // 	if(objectiveVector[0] < best+(best/100.0)){ //only if at least 1% worse
// // 	if(objectiveVector[0] < 1999){ //only if at least 1% worse
// // 		printf("%f smaller than %f, best: %f\n", objectiveVector[0],1999.0, best);
// // 		printf("%f --> ", objectiveVector[0]);
// 		int rc = cobyla(decisionNumber, 2*decisionNumber, decisionVector, rhobeg, rhoend, iprint, &maxfun, calcfc, &state);

	//**********************************************************************************************************************//
// 	double tof4=(decisionVector[4]*(2500-20))+20;
// 	double alpha=(decisionVector[8]*(0.99-0.01))+0.01;
// 	double time = tof4 * (1 - alpha) * 86400;
	
	bool violate=false;
	for(int i=0;i<decisionNumber;i++){
		if(decisionVector[i] > superiorPositionLimit[i]){
			violate=true;
			break;
		}
		if(decisionVector[i] < inferiorPositionLimit[i]){
			violate=true;
			break;
		}
	}
	if(!violate){
		obj = tandem_wbounds(decisionVector,tof,seq);
// 		if(time <= 0)
// 			printf("lambert error on time -- %f \t %f\n", tof4, alpha);
		objectiveVector[0]=2000-obj;
		objectiveVector[1] = tof;
	}else
		objectiveVector[0]=objectiveVector[1]=4000;
	// 		printf("%f\n", objectiveVector[0]);

	if(bestTandem > objectiveVector[0]){
		bestTandem = objectiveVector[0];
// 		printf("New Best: %f --> ", bestTandem);
// 		printVector(decisionVector, decisionNumber);
// 		printf("\n");
	}
	
	
	if(objectiveVector[1] != objectiveVector[1]){
		printf("\nError on obj1: ");
		printVector(objectiveVector, objectiveNumber);
		printf(" --> ");
		printVector(decisionVector, decisionNumber);
		printf("\n");
// 		exit(1);
		objectiveVector[1]=4000;
	}
	
	if( objectiveVector[0] == 0 ){
// 	if(objectiveVector[0] == 0 &&  objectiveVector[0]/objectiveVector[0] == 1){ //can be both at the same time? only if it is nan. Works better than std::isnan() that is not compliant with --fast-math
		printf("Error on obj0: ");
		printVector(objectiveVector, objectiveNumber);
		printf(" --> ");
		printVector(decisionVector, decisionNumber);
		printf("\n");
// 		exit(1);
		objectiveVector[0]=4000;
	}
	
	if( objectiveVector[0] != objectiveVector[0] ){
		// 	if(objectiveVector[0] == 0 &&  objectiveVector[0]/objectiveVector[0] == 1){ //can be both at the same time? only if it is nan. Works better than std::isnan() that is not compliant with --fast-math
		printf("Error on obj0: ");
		printVector(objectiveVector, objectiveNumber);
		printf(" --> ");
		printVector(decisionVector, decisionNumber);
		printf("\n");
		// 		exit(1);
		objectiveVector[0]=4000;
	}
}
