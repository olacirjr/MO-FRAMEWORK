#include"trajobjfuns.h"
#include<vector>
#include<iostream>

int main(){
	
	//Put here your solution
	
	double sol[] = {8614.67190542681782972068, 3.23881168881835801443, 0.49828711089951693847, 0.55247385652360181396, 2323.11617248155198467430, 1911.28129676196658692788, 2104.75935736787505447865, 2499.99977990785100701032, 0.90606524729349469105, 0.83398203560248784783, 0.77167564570127900048, 0.18281414705776333207, 1.41932189054689850138, 1.46704208497941945843, 1.05000000000010906831, -1.24053864588253404122, -1.62734953864368758758, -1.25424110003836464244}; 
	const int seq[] = {3, 2, 3, 3, 6}; // Earth, Venus, Earth,Earth, Staturn
	int dimension = sizeof(sol) / sizeof(sol[0]);

	//Casting the array to a vector
	std::vector<double> X;
	for (int i = 0; i < dimension; i++){
		X.push_back(sol[i]);
	}

	double tof;
	//Evaluating the point:
	double obj = tandem(X,tof,seq);

	//Printing the point
	for (int i = 0; i < dimension; i++){
		std::cout << "No.: " << i << " = =" << X[i] << std::endl;
	}
	
	//Printing the result:
	std::cout << "Function Value: " << obj << std::endl;




}
