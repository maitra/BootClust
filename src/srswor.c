
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
#include "order.h"
#include "array.h"
#define MATHLIB_STANDALONE 1 /*It is essential to have this before the call 
                               to the Rmath's header file because this decides
                               the definitions to be set. */
#include <Rmath.h>

/* Note that use of this function involves a prior call to the Rmath library to
   get the seeds in place. It is assumed that this is done in the calling 
   function. */

int srswor(int n,int *ranordr)
{
  /*Given a sample of size n, provide a random permutation ordering in the 
    integer array ranordr.*/
  double *x;
  size_t *ord,i;
  MAKE_VECTOR(x,n);
  for(i=0;i<n;i++) x[i]=runif(0.,1.);
  ord=orderDouble(x,n);
  FREE_VECTOR(x);
  for(i=0;i<n;i++) ranordr[i]=ord[i];
  FREE_VECTOR(ord);
  return 0;
}

int WRSampleUnequalProb(int n, int k, double *prob, int *smpl)
{
  /*
    Given n and k, provide a random sample with replacement of size n from 
    k classes, each occurring with unequal probabilities, or frequencies.
    Input parameters are as follows:

    n = sample size
    k = numer of class ids
    prob = probability of occurrence of class ids (also can handle frequencies)
    smpl = n-variate vector containing ids (0 through k-1) of the sample

    written by Ranjan Maitra, January 21, 2007.  
*/

  double *cdf, *sam, *ordsam;
  int i, j;
  size_t *ord;

  MAKE_VECTOR(cdf, k);
  cdf[0] = prob[0];
  for (i = 1; i < k; i++) {
    cdf[i] = cdf[i-1];
    cdf[i] += prob[i];
  }
  if (cdf[k-1] != 1) { 
    for (i = 1; i < k; i++) cdf[i] /= cdf[k-1];
  }

  MAKE_VECTOR(sam, n);
  for (i=0; i<n; i++) sam[i] = runif(0., 1.);
  ord = orderDouble(sam, n);

  MAKE_VECTOR(ordsam, n);
  for (i = 0; i < n; i++) ordsam[i] = sam[ord[i]];
  j=0;
  for (i = 0; i < k ; i++) {
    for (; (j < n) && (ordsam[j] < cdf[i]) ; j++) {
      sam[j] = i;
    }
  }
  for (i = 0; i < n; i++)  smpl[i] = sam[ord[i]];
  FREE_VECTOR(sam);
  FREE_VECTOR(ord);
  FREE_VECTOR(ordsam);
  FREE_VECTOR(cdf);
  return 0;
}
