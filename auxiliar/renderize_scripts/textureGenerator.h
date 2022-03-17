#include <iostream>
#include <string.h>
#include <math.h>

void info(char **argv);
void readfile(FILE *input, double *n, int *pt, int *N);
void molp(double* lightInt,double* n,int* pt,int* N);
void printData(double *lightInt, int *N, char fname[]);
void gera_plot (char arquivo[]);
