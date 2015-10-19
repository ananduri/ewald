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

double B(double r, double alpha);
double C(double r, double alpha);

double structure2(double x, double y, double z, double i, double j, double k);

double realsum(double x, double y, double z, double alpha, int recip_cut, int real_cut, double cellsize, bool charray[]);
double recsum(double x, double y, double z, double alpha, int recip_cut, int real_cut, double cellsize, bool charray[]);

double selfint(double alpha, double cellsize);



#endif /* EEHEAD_H_ */
