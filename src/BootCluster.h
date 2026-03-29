#ifndef BOOTCLUSTER_H
#define BOOTCLUSTER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "array.h"
#include <assert.h>
#include "mat_vec-VM.h"
#include "statutils.h"

void hc(int n, int m, int iopt, double **data, int *ia, int *ib, double *crit);
void kmeans(double **a, int m, int n, double **c, int k, int *ic1,int *nc,
	    int iter,double *wss,int *ifault);
void EigValDec(int size, double *W, double **A, double (*determinant));
void hclassify(int n,int m, double **x,int hcrit,int nclass,int *class);
void kmeans_hclustgroup(double **x, int n, int p, int nclass, double **Mu, 
			int *ia, int *ib, int iter, double *wss, int *nc, 
			int *class);
double kmnsgroup(double **x,int n,int p,int nclass,double **Mu, int *nc, 
		 int *class);

double HierClust(double **X, int N, int p, double **c, int K, int *ic1, int *nc,
		 double *wss, int hcrit);

double RunKMeansNew(double **X, int N, int p, double **c, int K, int *ic1, int *nc, int iter, double *wss, int *ifault, int hcrit, int flag);

double RunKMeans(double **X, int N, int p, double **c, int K, int *ic1, 
		 int *nc, int iter, double *wss, int *ifault, int Nrepl, 
		 int flag, double **means, short usemeans);

void StandardErr(double **E, int *Xic1, int N, int p, double **VarCov, int K, double **DeSt);
void DeStandardErr(double **W, int *Xic1, int N, int p, double **DeSt, int K);
double TestKvsKstar_HIER(char *outdir, int K, int Kstar, double **X, int N, int p, 
			 int RefDistN, int Nrepl, int iter, int hcrit);
double TestKvsKstar_SPH(char *outdir, int K, int Kstar, double **X, int N, int p, 
			int RefDistN, int Nrepl, int iter, int hcrit);

int TestForward(char *outdir, int minK, int maxK, double **X, int N, int p, int RefDistN, 
		int Nrepl, int iter, double alpha, double *PVforw, int hcrit, 
		short meth);

#endif /* BOOTCLUSTER_H */
