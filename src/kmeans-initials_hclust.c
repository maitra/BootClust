/*

  Bootstrapping for clustering significance
  Copyright (C) 2012  Ranjan Maitra and Volodymyr Melnykov

    DISCLAIMER: THIS SOURCE CODE IS SUPPLIED "AS IS" WITHOUT WARRANTY OF
    ANY KIND, AND ITS AUTHOR AND THE JOURNAL OF MACHINE LEARNING RESEARCH
    (JMLR) AND JMLR'S PUBLISHERS AND DISTRIBUTORS, DISCLAIM ANY AND ALL
    WARRANTIES, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES OF
    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, AND ANY
    WARRANTIES OR NON INFRINGEMENT. THE USER ASSUMES ALL LIABILITY AND
    RESPONSIBILITY FOR USE OF THIS SOURCE CODE, AND NEITHER THE AUTHOR NOR
    JMLR, NOR JMLR'S PUBLISHERS AND DISTRIBUTORS, WILL BE LIABLE FOR
    DAMAGES OF ANY KIND RESULTING FROM ITS USE. Without limiting the
    generality of the foregoing, neither the author, nor JMLR, nor JMLR's
    publishers and distributors, warrant that the Source Code will be
    error-free, will operate without interruption, or will meet the needs
    of the user.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Authors' contact information:
    Ranjan Maitra						Volodymyr Melnykov
    maitra@iastate.edu					vmelnykov@ua.edu
    Department of Statistics			Department of Information Systems,
    Iowas State Unviersity				Statistics, and Management Science
	Ames, IA 50011						The University of Alabama
										Tuscaloosa, AL 35487
*/


#include<stdlib.h>
#include<stdio.h>
#include "array.h"
#include<math.h>

#define PI 3.141593
#define Inf 1e+140

void hclass(int n, int *ia, int *ib, int lev, int *iclass);

void kmeans(double **a, int m, int n, double **c, int k, int *ic1,int *nc,
	    int iter,double *wss,int *ifault);

void meansonly(double **x, int n, int p, double *mu)
{
  /* This routine calculates (only) the means of a homogeneous sample */
	int i, j;
	
	if ((n==0) || (n==1)) {
		if (n==1) {
			for(i=0;i<p;i++) {
				mu[i]=x[0][i];
			}
		}
		else {
			for(i=0;i<p;i++) {
				mu[i]=0.0;
			}
		}
	}
	else {
		for (i=0;i<p;i++) {
			mu[i]=0.;
		}
		for (i=0;i<n;i++) {
			for (j=0;j<p;j++) {
				mu[j]+=x[i][j];
			}
		}
		for (i=0;i<p;i++) {
			mu[i]/=n;
		}
	}
	return;
}

void initialmeans(double **x,int n,int p,int nclass,int *nc,
		  double **Mu,int *class)
{
	double **y;
	int i, j, l;
	MAKE_MATRIX(y, n, p);
	for (i = 0; i < nclass ; i++) {
		nc[i] = 0;
		for(l = 0; l < n; l++) {
			if (class[l]==i) {
				for (j = 0; j < p; j++) {
					y[nc[i]][j]=x[l][j];
				}
				nc[i]++;
			}
		}
		meansonly(y, nc[i], p, Mu[i]);
	}
	FREE_MATRIX(y);
	return;
}

void kmeanshcstarters(double **x, int n, int p, int nclass, double **Mu, 
		      int *ia, int *ib)
/* This routine takes the output from hc and cuts the tree into the first 
   nclass groups giving nclass groups. Initial estimates of Mu are obtained
   from these groups. */
{
	int *nc, *class;
  
	MAKE_VECTOR(nc, nclass);
	MAKE_VECTOR(class, n);

	hclass(n, ia, ib, nclass, class);

	initialmeans(x, n, p, nclass, nc, Mu, class);
	
	FREE_VECTOR(nc);
	FREE_VECTOR(class);

	return;
}

void kmeans_hclustgroup(double **x, int n, int p, int nclass, double **Mu,
		       int *ia, int *ib, int iter, double *wss, int *nc, 
		       int *class)

/* This function returns the value of a call to ifault. It initializes given 
   the output of hc.c and then uses kmeans. For kmeans, the best choice for 
   a hierarchically clustered tree is given by hclust with Ward's criterion. */
{
	int ifault;
	
	kmeanshcstarters(x, n, p, nclass, Mu, ia, ib);
	
	kmeans(x, n, p, Mu, nclass, class, nc, iter, wss, &ifault);
	
	return;

}
