/*
 *  GOProblems.h
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/17/08.
 *  Copyright 2008 ¿dvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#ifndef TRAJOBJFUNS_H
#define TRAJOBJFUNS_H
#include<vector>

const double MinTandemBounds[] = {5475.0,2.5,0,0,20,20,20,20,0.01,0.01,0.01,0.01,1.05,1.05,1.05,-pi_const, -pi_const,-pi_const};
const double MaxTandemBounds[] = {9131,4.9,1,1,2500,2500,2500,2500,0.99,0.99,0.99,0.99,10,10,10,pi_const,pi_const,pi_const};

//NOTE: the functions here have passing by reference + const as they are called a lot of time during execution and thus
//it is worth trying to save time by avoiding to make a copy of the variable passed

double gtoc1        (const std::vector<double>& x, std::vector<double>& rp);
double cassini1     (const std::vector<double>& x, std::vector<double>& rp);
double sagas        (const std::vector<double>& x, double& DVtot, double& DVonboard);
double rosetta      (const std::vector<double>& x);
double cassini2     (const std::vector<double>& x);
double messenger    (const std::vector<double>& x);
double messengerfull(const std::vector<double>& x);
double tandem       (const std::vector<double>& x, double& tof, const int sequence_[]);
double tandem_wbounds      (const std::vector<double>& x, double& tof, const int sequence_[]);


#endif
