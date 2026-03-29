

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

#include <stdio.h> 
#include <stdlib.h> 
#include "array.h"

#define MIN(i,j) ((i)<(j) ? (i) : (j))
#define MAX(i,j) ((i)>(j) ? (i) : (j))

void dgesdd_(char *jobz,int *m,int *n,double *a,int *lda,double *s,double *u,
	     int *ldu,double *vt,int *ldvt,double *work,int *lwork,int *iwork,
	     int *info);
int svdd(double **a, int m, int n, double *d, double **u, double **v)
{
	double *A, *U, *VT;
	int lwork = -1;
	int liwork = 8*MIN(m,n);
	char jobz = 'S';
	double dw,*work=&dw; /*points to a temporary cell*/
	int *iwork;
	int i, j, k, info,minmn=MIN(m,n);

	MAKE_VECTOR(A, m*n);
	for (j=0, k=0; j<n; j++) {
		for (i=0; i<m; i++) A[k++] = a[i][j];
	}

	MAKE_VECTOR(U, m*minmn);
	MAKE_VECTOR(VT,minmn*n);
	MAKE_VECTOR(iwork, liwork);

	lwork=-1;
	
	dgesdd_(&jobz, &m, &n, A, &m, d, U, &m, VT, &n,
		work, &lwork, iwork, &info); /*call to get optimal lwork*/
	
	if (info!=0) {
	  printf("error: allocating LWORK in svdd\n");
	  exit(1);
	}

	lwork=(int)*work;
	MAKE_VECTOR(work, lwork);

	dgesdd_(&jobz, &m, &n, A, &m, d, U, &m, VT, &minmn,
			work, &lwork, iwork, &info);
	FREE_VECTOR(A);
	FREE_VECTOR(work);
	FREE_VECTOR(iwork);

	for (j=0, k=0; j<minmn; j++) {
	  for (i=0; i<m; i++) u[i][j]=U[k++];
	}

	/* VT, as calculated by dgesdd_(), is the transpose of the right
	 * multiplier.  Here we undo the transpose so that the matrix
	 * v[][] returned by this function is not transposed anymore.
	 */

	for (i=0,k=0; i<n; i++) {
	  for (j=0; j<minmn; j++)   v[i][j]=VT[k++];
	}

	FREE_VECTOR(U);
	FREE_VECTOR(VT);

	return info;
}


