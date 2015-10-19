/*
 * ewald.h
 *
 *  Created on: Aug 24, 2015
 *      Author: Arun
 */

#ifndef EWHEAD_H_
#define EWHEAD_H_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

const double J = 3.72; 
const double D = 1.41;
const double rnn = 0.25*sqrt(2);

inline double structure2(double x, double y, double z, double i, double j, double k)
{
	return cos(M_PI*2.0*(x*((double) i) + y*((double) j) + z*((double) k)));
}

const double is3 = 1/sqrt(3);

inline double B(double r, double alpha){
	return ( (1 - erf(alpha*r))/(r*r*r) + (2*alpha/sqrt(M_PI))*exp(-alpha*alpha*r*r)/(r*r) );
}

inline double C(double r, double alpha){
	return ( 3*(1 - erf(alpha*r))/(r*r*r*r*r) + (2*alpha/sqrt(M_PI))*(2*alpha*alpha + 3/(r*r))*exp(-alpha*alpha*r*r)/(r*r) );
}

double realsum(double x, double y, double z, int m, double alpha, int real_cut, double cellsize, int bsize, double* NNenergy, double* intmat);
double recsum(double x, double y, double z, int m, double alpha, int recip_cut, double cellsize, int bsize, double* intmat);

double selfint(double alpha, double cellsize, int bsize);

#endif /* EWHEAD_H_ */
