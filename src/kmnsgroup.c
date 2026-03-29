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
  This routine clusters a dataset into groups using the kmeans algorithm with 
  start points provided by the starts_via_svd routine. The function returns 
  the classification ids as well as the estimated centers.
*/

#include<stdlib.h>
#include<stdio.h>
#include "array.h"
#include<math.h>
#define PI 3.141593
#define Inf 1e+140

void kmeans(double **a, int m, int n, double **c, int k, int *ic1,int *nc,
	    int iter,double *wss,int *ifault);
int starts_via_svd(int n,int m,double **Mu,double **x,int nclus,int *ningrp,
		    double *pi,int *grpids,double **LTSigma,double alpha,
		    int llhdnotW);
void initials(double **x,int n,int p,int nclass,int *nc,
              double **Mu,double **LTSigma,int *class);
double determinant(double *LTSigma,int n);

double trW(double *wss, int K);

double kmnsgroup(double **x,int n,int p,int nclass,double **Mu,int *nc,
		 int *class)
{
  int j,flag=0,ifault;
  double *pi,**LTSigma,*WSS,detw,*temp;
  double trOfW;
  
  MAKE_VECTOR(pi,nclass);
  MAKE_MATRIX(LTSigma,nclass,p*(p+1)/2);
  if (nclass==1) {
    double *temp;
    nc[nclass-1]=n;
    pi[nclass-1]=1.0;
    for (j=0;j<n;j++) class[j]=0;
    initials(x,n,p,nclass,nc,Mu,LTSigma,class);
    MAKE_VECTOR(temp,p*(p+1)/2);
    for (j=0;j<(p*(p+1)/2);j++) temp[j]=(n-1)*LTSigma[0][j];
    detw=log(determinant(temp,p)); 
    FREE_VECTOR(temp);
  }
  else {
    double *pi;
    MAKE_VECTOR(pi,nclass);
    j=starts_via_svd(n,p,Mu,x,nclass,nc,pi,class,LTSigma,0.99,0);
    FREE_VECTOR(pi);
    for(j=0;((j<nclass) && (flag==0));j++) {
      if(nc[j]==0) flag=1;
    }   
    if (flag==1) {
      printf("Note: singular starting point");
      detw=Inf;
    }
    else {
      int i;
      MAKE_VECTOR(WSS,nclass);
      kmeans(x,n,p,Mu,nclass,class,nc,100000,WSS,&ifault);
      trOfW = trW(WSS,nclass);
      initials(x,n,p,nclass,nc,Mu,LTSigma,class);
      MAKE_VECTOR(temp,p*(p+1)/2);
      for(j=0;j<p*(p+1)/2;j++)  temp[j]=0;
      for(i=0;i<nclass;i++) {
	for (j=0;j<(p*(p+1)/2);j++) temp[j]+=(nc[i]-1)*LTSigma[i][j];
      }
      detw=log(determinant(temp,p)); 
      FREE_VECTOR(temp);
      FREE_VECTOR(WSS);
    }
  }
  return trOfW;
} 

