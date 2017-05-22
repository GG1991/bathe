
#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int cmp_int(void *a, void *b);
int cmp_dou(void *a, void *b);
int strBin2Dec(char *str, int *dec);
int voigt2mat(double strvoi[6],double strmat[3][3]);
double normvec(double * vec, int n);

#endif
