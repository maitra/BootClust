#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <ctype.h>

void usage(char *s)
{
	fprintf(stderr, "NAME\n %s - Run the Boostrapped Clustering algorithm of Maitra et al, JASA 2012. \n  SYNOPSIS \n \t%s  -i <file> -O <dir> -n <int> -p <int> -# <int> -N -k -K -s -q -m -# -l -S \n", s, s);
	fprintf(stderr, "\nOPTIONS\n \t-i <file> datafile (with n (row) observations of p (column) attributes\n \t-O <dir>  output directory\n \t-n <int>  number of observations (rows) \n \t-p <int>  number of dimensions (columns) \n \t-N <int>  number of bootstrap replications (1000 by default)\n \t-k <int> (minimum) number of clusters under the null hypothesis\n \t-K <int> (maximum) number of clusters under the alternative hypothesis\n \t-s        if series of tests should be performed (none by default)\n \t-q <real> FDR level (for the sequence of tests)\n \t-m <0|1>  if kmeans or hierarchical clustering is the clustering algorithm\n \t-# <int>  number of initializations for k-means (if kmeans chosen, default 0 which equates to n * p * k)\n \t-l <int>  linkage for hierarchical clustering (if chosen) with:");
	fprintf(stderr, " \n \t \t           1 = Ward (default) \n \t \t           2 = single \n \t \t           3 = complete \n \t \t           4 = average \n \t \t           5 = McQuitty \n \t \t           6 = median \n \t \t           7 = centroid \n \t -S         should observations in each coordinate be scaled by their ranges (maximum - minimum values) (default: no)\n ");
} /* print_usage */

short read_options(int agc, char **agv, char **infil, char **oudir, int *n, int *p, int *k, int *K, int *RefDistN, int *N, int *s, double *q, short *meth, int *linkage, short *scaled) {
	char c;
	short verbose = 0;

	*s = 0;
	*meth = 0;
	*scaled = 0;
	*linkage = 1;
	*RefDistN = 1000;
	*N = 0;

	while ((c = getopt(agc, agv, "vi:O:n:p:B:k:K:s:q:m:#:l:S")) != EOF) {
		switch (c) {
		case 'v':
			verbose = 1;
			break;
		case 'i':
			*infil = optarg;
			break;
		case 'O':
 			*oudir = optarg;
			break;
		case 'n':
			*n = atoi(optarg);
			break;
		case 'p':
			*p = atoi(optarg);
			break;
		case 'B':
			*RefDistN = atoi(optarg);
			break;
		case 'k':
			*k = atoi(optarg);
			break;
		case 'K':
			*K = atoi(optarg);
			break;
		case 's':
			*s = 1;
			break;
		case 'q':
			*q = atof(optarg);
			break;
		case 'm':
			*meth = atoi(optarg);
			break;
		case '#':
			*N = atoi(optarg);
			break;
 		case 'l':
			*linkage = atoi(optarg);
			break;
		case 'S':
			*scaled = 1;
			break;
		case '?':
			if (optopt == 'i' || optopt == 'O' || optopt == 'n' || optopt == 'p' || optopt == 'B' || optopt == 'k' || optopt == 'K' || optopt == 'N' || optopt == '#' || optopt == 'q' || optopt == 'm')
				fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else 
				if ((*meth == 1) && (optopt == 'l'))
					fprintf (stderr, "Option -%c requires an argument for hierarchical clustering.\n", optopt);
				else
					if (isprint (optopt))
						fprintf (stderr, "Unknown option `-%c'.\n", optopt);
					else 
						fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
			usage(agv[0]);
			exit(1);
		default:
			fprintf(stderr, "-i -O -n -p -K needed, -v -B -k -s -q -l -m -S optional\n");
			usage(agv[0]);
			exit(1);
		}
	}
	if (!*infil || !*oudir || !*K || !*n || !*p) {
		printf("Missing mandatory arguments\n");
		usage(agv[0]);
		exit(1);
	}
	if (verbose) {
		fprintf (stderr, "input file = %s, output directory = %s\n", *infil, *oudir);
		fprintf(stderr, "n = %d, p = %d, k = %d, K = %d, B = %d, # = %d, s = %d, q = %f, m = %d, l = %d\n", *n, *p, *k, *K, *RefDistN, *N, *s, *q, *meth, *linkage);
	}
	return verbose;
}
