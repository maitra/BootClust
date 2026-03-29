

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


/*
  This routine clusters a dataset into groups using the E-M algorithm with 
  start points provided by the starts_via_svd routine. The function returns 
  the classification ids as well as the parameter estimates.
*/

#include<stdlib.h>
#include<stdio.h>
#include "array.h"
#include<math.h>
#define PI 3.141593
#define Inf 1e+140

double determinant(double *LTSigma,int n);
void emcluster(int n,int p,int k,double *pi,double **X,double **Mu, 
	       double **LTSigma,int maxiter,double eps,double *llhdval);
int starts_via_svd(int n,int m,double **Mu,double **x,int nclus,int *ningrp,
                    double *pi,int *grpids,double **LTSigma,double alpha,
                    int llhdnotW);
void assign(int n, int p,int k,double **X,double *pi,double **Mu,
	    double **LTSigma,int *class,int *nc);
void meandispersion(double **x, int n, int p, double *mu, double *ltsigma);
int emgroup(double **x,int n,int p,int nclass,double *pi,double **Mu,
	     double **LTSigma,double *llhdval,int *nc,int *class)
{
  int j,flag=0;
  double like;
  
  if (nclass==1) {
    nc[nclass-1]=n;
    pi[nclass-1]=1.0;
    for (j=0;j<n;j++) class[j]=0;
    meandispersion(x,n,p,Mu[0],LTSigma[0]);
    like=-0.5*n*p-0.5*n*log(determinant(LTSigma[0],p))-0.5*n*p*log(2*PI); 
  }
  else {
    if(!starts_via_svd(n,p,Mu,x,nclass,nc,pi,class,LTSigma,0.99,1))      {
      for(j=0;j<nclass;j++) pi[j]=nc[j]/(double)n;
      emcluster(n,p,nclass,pi,x,Mu,LTSigma,1000,0.0001,&like);
      assign(n,p,nclass,x,pi,Mu,LTSigma,class,nc);
    }
    else flag=1;
  }
  (*llhdval)=like;
  printf("like =  %f\n",like);
  return flag;
} 

