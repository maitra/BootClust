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
    Ranjan Maitra			Volodymyr Melnykov
    maitra@iastate.edu			vmelnykov@ua.edu
    Department of Statistics		Department of Information Systems,
    Iowas State Unviersity	 	Statistics, and Management Science
    Ames, IA 50011			The University of Alabama
					Tuscaloosa, AL 35487
*/

#include "BootCluster.h"
#include "readopt.h"
#include <getopt.h>
#include <ctype.h>

int main(int argc, char **argv)
{
	int N=0, p=0, K, RefDistN = 1000, Nrepl = 1000, minK=1, 
		maxK=0, i, j, hcrit = 1, iter = 1000;
	double **X, *PVforw, alpha = 0.05;
	FILE *fdata, *fout;
	char fname[256];
	char *filename = NULL;	 /* input filename path */
	char *output_dir = NULL; /* directory for output */
	short verbose, scaled, meth;
	int seq;
	
	verbose = read_options(argc, argv, &filename, &output_dir, &N, &p, &minK, &maxK, &RefDistN, &Nrepl, &seq, &alpha, &meth, &hcrit, &scaled);

/*	iter :  # of iterations for kmeans 
	Nrepl : # of initial trials of kmeans 
	N :     Sample size 
	p :     Dimension 
	RefDistN:# of points of Reference distribution 
	minK:	Min # of clusters (under H0) 
	maxK:   Max # of clusters (under Ha) 
	alpha:  FDR (or p-value for rejection)
	hcrit:  Hierarchical clustering criterion 
	1 - Ward, 2 - single, 3 - complete, 4 - average, 5 - McQuitty, 6 - median, 7 - centroid 
	meth  meth = 0 is spherical (kmeans), 1 is hierarchical */

	sprintf(fname, "%s", filename); /* file name */
	
	MAKE_MATRIX(X, N, p);
	MAKE_VECTOR(PVforw, maxK-1);

	if (verbose) 
		printf("Dataset file name: %s:\n", fname);

	fdata = fopen(fname,"r");
	if (fdata == NULL) {
		fprintf(stderr, "Missing input file: %s. Does it exist?\n", fname);
		exit(1);
	}
	if (verbose)
		fprintf(stderr, "Reading in data from file: %s:\n", fname);
	for (i=0;i<N;i++) 		
		for (j=0;j<p;j++) 
			fscanf(fdata,"%lf ", &X[i][j]);
	fclose(fdata);	
	if (verbose)
		fprintf(stderr, "Done reading in datam from file: %s:\n", fname);
	if (scaled) {
		double ra, *y;
		MAKE_VECTOR(y, N);
		for (i =0; i < p; i++) {
			for (j = 0; j < N; j++)
				y[j] = X[j][i];
			ra = sd(y, N);
			if (ra != 0)
				for (j = 0; j < N; j++)
					X[j][i] /= ra;
		}
		FREE_VECTOR(y);
	}
	
	K = TestForward(output_dir, minK, maxK, X, N, p, RefDistN, Nrepl, iter, alpha, PVforw, hcrit, meth);

	sprintf(fname, "%s/pvals.dat", output_dir); /* file name */
	fout = fopen(fname, "w");
	for (i=0;i<maxK-1;i++) 
		fprintf(fout,"%f ", PVforw[i]);
	fprintf(fout,"\n");
	fclose(fout);

	printf("   Result: %i clusters have been detected\n", K);

	FREE_VECTOR(PVforw);
	FREE_MATRIX(X);
	return EXIT_SUCCESS;

}
