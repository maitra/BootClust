
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
#include <string.h>
#include "array.h"

void dsyev_(char *JOBZp, char *UPLOp,int *Np, double *A, int *LDAp, double *Wp, double *WORK, int *LWORK, int *INFOp);



void EigValDec(int size, double *W, double **A, double (*determinant))

/*
	returns:
	W - vector of eigenvalues
	A - matrix of eigenvectors 
*/

{
  int i, j, INFO, N, LDA;
  char uplo='L';
  double *AT;
  char JOBZ='V';
  double *WORK;
  int LWORK;

  MAKE_VECTOR(AT,size*size);
  for (i=0; i<size; i++){
    for(j=0; j<size; j++) AT[j+size*i]=A[j][i];
  }

  N=size;
  LDA=size;

  LWORK=3*size-1;
  MAKE_VECTOR(WORK,LWORK);

  dsyev_ (&JOBZ, &uplo, &N, AT, &LDA, W, WORK, &LWORK, &INFO);

  if (INFO==0){
    int i;
    (*determinant)=1.0;
    for (i=0;i<N;i++){
      (*determinant)*=W[i];
    }
/* 
   printf("Eigenvalues:\n ");
    for (i=0; i<size; i++){
      printf("%f \n",W[i]);
    }
*/

  }

  for (i=0; i<size; i++){
    for(j=0; j<size; j++) A[j][i]=AT[j+size*i];
  }

  if (INFO!=0){
      printf("Problem in EigValDec:  error %d\n",INFO);
  }

  FREE_VECTOR(AT);
  FREE_VECTOR(WORK);

  return;
}

