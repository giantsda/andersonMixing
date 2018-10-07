/*
 * NR_chen.h
 *
 *  Created on: Oct 6, 2018
 *      Author: chen
 */

#ifndef NR_CHEN_H_
#define NR_CHEN_H_

#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include <stdlib.h>

struct parameterDumper
{
  int i1, i2, i3, i4, i5;
  double d1, d2, d3, d4, d5;
  int *pi1, *pi2, *pi3, *pi4, *pi5;
  double* pd1, *pd2, *pd3, *pd4, *pd5;
  char* s1, *s2, *s3, *s4, *s5;
};

int
adm_chen (void
(*f) (double*, double*, int, struct parameterDumper*),
	  double* x_old, double tol, int maxIteration, int n,
	  struct parameterDumper* p);

#endif /* NR_CHEN_H_ */
