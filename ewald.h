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

const double J = 3.72; //before had the negative sign
const double D = 1.41;
const double rnn = 0.25*sqrt(2);


double B(double r, double alpha);
double C(double r, double alpha);

double structure2(double x, double y, double z, double i, double j, double k);

double realsum(double x, double y, double z, double m, double alpha, int real_cut, double cellsize, int bsize, bool charray[], double* NNenergy, double* intmat);
double recsum(double x, double y, double z, double m, double alpha, int recip_cut, double cellsize, int bsize, bool charray[], double* intmat);

double selfint(double alpha, double cellsize, int bsize);

void boffsetassign(double* boff, int m);
void dipassign(double* dip, int m, bool ori);

#endif /* EWHEAD_H_ */
