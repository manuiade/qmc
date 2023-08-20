#ifndef QMC_H
#define QMC_H
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define ING 6 //Max dimension imposed from stack
#define DIM 1 << ING
#define TRUE 1
#define FALSE 0

int min_terms[DIM];
int tab[DIM*DIM][ING];
int flag[DIM];
int reduction[DIM*DIM][ING];
int not_ones_1[DIM];
int not_ones_2[DIM];
int range[DIM*DIM][DIM*DIM];
int intmp[DIM*DIM][DIM*DIM];
int matrbin[DIM*DIM][DIM];
int mincolumn[DIM*DIM];
int check[DIM*DIM];
int prime_ess[DIM*DIM];
int solution[DIM*DIM];
int indiff[DIM];
int tabind[DIM*DIM][ING];

extern void upperTerm(int x, int pos, int num);
extern void lowerTerm(int pos, int num);
extern void outputTerm(int x, int pos, int num);
extern void swap(int i, int j, int k,int num);
extern void sort(int nmin, int nind,int num);
extern void fillRange(int a, int b, int step, int *dimint);
extern int hemming(int*tab1, int*tab2, int x, int a, int b, int step, int num, int *dimint);
extern int copy_minterms(int *x, int nmin,int num);
extern int fill_matrix(int x, int nmin, int step);
extern void bab(int s, int nmin, int x, int best, int*curr);
extern void printsp(int x,int num);
extern void printps(int x, int num);
extern void printrid(int x, int passo,int num);

#endif
