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




#include "statutils.h"

#define SQ(x) ((x) * (x))

int srswor(int n, int k, int *y)
{
	/* Provide k out of n indices sampled at random without replacement */
	
	if (k > n) {
		printf("Error: k greater than n in srswor()");
		return 1;
	}
	else {
		
		int i, j, *x;
		
		MAKE_VECTOR(x, n);
		for (i = 0; i < n; i++)	x[i] = i;
		
		for (i = 0; i < k; i++) {
			j = n * runif(0, 1);
			y[i] = x[j];
			x[j] = x[--n];
		}
		FREE_VECTOR(x);
	}
	return 0;
}


void SampleVarCov(double **X, int N, int p, int K, double **c, int *ic1, int *nc, double **VarCov){


	int i,j,k,h,sch;
	int m;


	m = p * (p + 1) / 2;

	for (k = 0; k < K; k++) {

		for (sch = 0; sch < m; sch++) {

			VarCov[k][sch] = 0.0;

		}

	}


	sch = 0;

	for (h = 0; h < p; h++) {

		for (j = h; j < p; j++) {	

			if (h == j){

				for (i = 0; i < N; i++) {

					VarCov[ic1[i]][sch] = VarCov[ic1[i]][sch] + (X[i][j] - c[ic1[i]][j]) * (X[i][j] - c[ic1[i]][j]);

				}				

			} else {

				for (i = 0; i < N; i++) {

					VarCov[ic1[i]][sch] = VarCov[ic1[i]][sch] + (X[i][h] - c[ic1[i]][h]) * (X[i][j] - c[ic1[i]][j]);

				}				

			}

	
			for (k = 0; k < K; k++) {

				VarCov[k][sch] = VarCov[k][sch] / (nc[k] - 1);				

			}

			sch++;

		}
	}


	return;
}


double trW(double *wss, int K){

	double trW;
	int k;

	trW = 0;

	for (k=0;k<K;k++) {
	
		trW = trW + wss[k];

	}

	return trW;

}

double range(double *x, size_t n) {
	double minmax[2] = {INFINITY, -INFINITY};
	for (size_t i = 0; i < n; i++) {
		if (x[i] < minmax[0])
			minmax[0] = x[i];
		else 
			if (x[i] > minmax[1])
				minmax[1] = x[i];
	}
	return minmax[1] - minmax[0];
}

double sd(double *x, size_t n) {
	double mean = 0, sd;
	for (size_t i = 0; i < n; i++) 
		mean += x[i];
	mean /= n;

	sd = 0;
	for (size_t i = 0; i < n; i++) 
		sd += SQ(x[i] - mean);
	sd /= n - 1;
	return sqrt(sd);
}
