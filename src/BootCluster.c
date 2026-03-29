/*
  NOTES: by Ranjan Maitra
  changed k-means done using hierarchical clustering (with Ward's criterion --
  which I also changed) to k-means done using random initializations.

  This code now should correspond to the experiments done using Maitra et al,
  JASA 2012.) 
*/

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
    Iowas State Unviersity		Statistics, and Management Science
    Ames, IA 50011			The University of Alabama
					Tuscaloosa, AL 35487
*/

#include "BootCluster.h"

double HierClust(double **X, int N, int p, double **c, int K, int *ic1, int *nc, double *wss, int hcrit)
{
	double TrW;
	int i,j,k;
	
	hclassify(N, p, X, hcrit, K, ic1);

	for (k=0; k<K; k++){
		nc[k] = 0;
		wss[k] = 0.0;
		for (j=0; j<p; j++){
			c[k][j] = 0;
		}
	}

	for (i=0; i<N; i++){
		nc[ic1[i]]++;
		for (j=0; j<p; j++){
			c[ic1[i]][j] = c[ic1[i]][j] + X[i][j];
		}
	}

	for (k=0; k<K; k++){
		for (j=0; j<p; j++){
			c[k][j] = c[k][j] / nc[k];
		}
	}

	for (i=0; i<N; i++){
		for (j=0; j<p; j++){
			wss[ic1[i]] = wss[ic1[i]] + (X[i][j] - c[ic1[i]][j]) * (X[i][j] - c[ic1[i]][j]);
		}
	}

	TrW = trW(wss, K);
	return TrW;
}


double RunKMeansNew(double **X, int N, int p, double **c, int K, int *ic1, int *nc, int iter, double *wss, int *ifault, int hcrit, int flag)
{
	int *ia, *ib;
  	double *crit;
  	double minTrW, trWHier;
	double **minc;
	int *minic1, *minnc;
	int i, k, j;
  	
  	MAKE_VECTOR(ia, N);
  	MAKE_VECTOR(ib, N);

	MAKE_MATRIX(minc,K,p);
	MAKE_VECTOR(minic1,N);
	MAKE_VECTOR(minnc,K);

	
	if (K == 1){
			
		for (j=0;j<p;j++){

			c[0][j] = 0.0;

			for (i=0;i<N;i++){

				c[0][j] = c[0][j] + X[i][j];
				ic1[i] = 0;

			}

			c[0][j] = c[0][j] / N;
	
		}

		nc[0] = N;
		
		minTrW = 0.0;
		for (i=0; i<N; i++){
			for (j=0; j<p; j++){
				minTrW = minTrW + (X[i][j] - c[0][j]) * (X[i][j] - c[0][j]);
			}
		}
		
	} else {

		MAKE_VECTOR(crit, N);

		hc(N, p, hcrit, X, ia, ib, crit);
	
		FREE_VECTOR(crit);

		kmeans_hclustgroup(X, N, p, K, minc, ia, ib, iter, wss, minnc, 
				   minic1);

		trWHier = trW(wss,K);

		minTrW = 1e+20;
	
	
		if (trWHier < minTrW){
		
			minTrW = trWHier;
	
			if (flag == 1){
	
				for (k=0;k<K;k++){
	
					for (j=0;j<p;j++){

						c[k][j] = minc[k][j];

					}

					nc[k] = minnc[k];

				}

				for (i=0;i<N;i++){

					ic1[i] = minic1[i];

				}

			}	
	
		}

	}


  	FREE_VECTOR(ia);
  	FREE_VECTOR(ib);
  	
	FREE_MATRIX(minc);
	FREE_VECTOR(minic1);
	FREE_VECTOR(minnc);
			 
	return minTrW;

}



double RunKMeans(double **X, int N, int p, double **c, int K, int *ic1, 
		 int *nc, int iter, double *wss, int *ifault, int Nrepl, 
		 int flag, double **means, short usemeans)
{
/*
	Returns the sum of within-cluster sums of squares
	
	Parameters:
		X - data
		N - # of observations
		p - dimension of data
		K - # of clusters
		
		parameters of kmeans

		Nrepl - # if initial trials for kmeans
		flag - 0 if we just need the sum of within-cluster sums of squares
		     - 1 if we also need centers and classifications 
*/

	double **minc, *minwss, minTrW, curr;
	int *minic1, *minnc, *ind, code, i, j, k, t;

	if (!Nrepl)
		Nrepl = K * N * p;
	
	MAKE_MATRIX(minc,K,p);
	MAKE_VECTOR(minic1,N);
	MAKE_VECTOR(minnc,K);
	MAKE_VECTOR(minwss,K);
	MAKE_VECTOR(ind,K);
	
	minTrW = 1e20;
	
	if (K == 1) {
		for (j=0;j<p;j++){
			c[0][j] = 0.0;
			for (i=0;i<N;i++){
				c[0][j] = c[0][j] + X[i][j];
				ic1[i] = 0;
			}
			c[0][j] = c[0][j] / N;
		}
		nc[0] = N;
		minTrW = 0.0;
		for (i=0; i<N; i++)
			for (j=0; j<p; j++)
				minTrW = minTrW + (X[i][j] - c[0][j]) * (X[i][j] - c[0][j]);			
	}
	else {
/*		printf("K = %d, usemeans = %i\n", K, usemeans);*/
		if (usemeans) {
		        for (i = 0; i < K; i++)
				for (j = 0; j < p; j++)
					c[i][j] = means[i][j];
			kmeans(X, N, p, c, K, ic1, nc, iter, wss, ifault);
			minTrW = trW(wss, K);
		}
		else {
			for (t=0;t<Nrepl;t++){	
				code = srswor(N, K, ind);
				if (code != 0)
					printf("Problem with function srswor\n");
				for (k=0;k<K;k++)
					for (j=0;j<p;j++)
						c[k][j] = X[ind[k]][j];
				
				kmeans(X, N, p, c, K, ic1, nc, iter, wss, ifault);
				curr = trW(wss, K);
				if (curr < minTrW){
					minTrW = curr;					
					if (flag == 1){
						for (k=0; k<K; k++){
						for (j=0; j<p; j++)
							minc[k][j] = c[k][j];
						minnc[k] = nc[k];
						minwss[k] = wss[k];
						}	
						for (i=0; i<N; i++)
							minic1[i] = ic1[i];
					}
				}
			}
			if (flag == 1){
				for (k=0; k<K; k++){
					for (j=0; j<p; j++)
						c[k][j] = minc[k][j];
					nc[k] = minnc[k];
					wss[k] = minwss[k];
				}
				for (i=0;i<N;i++)
					ic1[i] = minic1[i];
			}
		}
	}
	FREE_MATRIX(minc);
	FREE_VECTOR(minic1);
	FREE_VECTOR(minnc);
	FREE_VECTOR(minwss);
	FREE_VECTOR(ind);
	
	return minTrW;
}

void StandardErr(double **E, int *Xic1, int N, int p, double **VarCov, int K, double **DeSt)
{

	double *Eig, *resid, *new,  **VC, **D, **De, **H, dtmt; 
	int i,j,h,k, m;

	m = p * (p+1) / 2;

	MAKE_VECTOR(Eig,p);
	MAKE_VECTOR(resid,p);
	MAKE_VECTOR(new,p);

	MAKE_MATRIX(VC,p,p);
	MAKE_MATRIX(D,p,p);
	MAKE_MATRIX(De,p,p);
	MAKE_MATRIX(H,p,p);


	for (k=0;k<K;k++){
	
		i = 0;
		j = -1;

		for (h=0;h<m;h++){

			j++;

			if (j == p){
			
				i++;
				j = i;

			}

			VC[i][j] = VarCov[k][h];

			if (i != j){
				
				VC[j][i] = VarCov[k][h];

			}
		}

		EigValDec(p,Eig,VC,&dtmt);


		for (i=0;i<p;i++){
			for (j=0;j<p;j++){

				if (i == j){
		
					if (Eig[i] != 0.0){
						D[i][j] = 1 / pow(Eig[i], 0.5);
					} else {
						D[i][j] = 0.0;
					}
					De[i][j] = pow(Eig[i], 0.5);	

				} else {

					D[i][j] = 0.0;		
					De[i][j] = 0.0;

				}
			}
		}



		multiplyAB(VC,p,p,D,p,p,H);		
		transpose(VC,p);
		multiplyAB(H,p,p,VC,p,p,D);


		transpose(VC,p);
		multiplyAB(VC,p,p,De,p,p,H);		
		transpose(VC,p);
		multiplyAB(H,p,p,VC,p,p,De);


		h = 0;
		for (i=0;i<p;i++){
			for (j=i;j<p;j++){

				DeSt[k][h] = De[i][j];
				h++;

			}
		}
		

		for (i=0;i<N;i++){

			if (Xic1[i] == k){

				for (j=0;j<p;j++){			

					resid[j] = E[i][j];

				}
				
				vecxmat(resid,p,D,p,p,new);


				for (j=0;j<p;j++){			

					E[i][j] = new[j];

				}

			}

		}
	


	}

	FREE_VECTOR(resid);
	FREE_VECTOR(new);
	FREE_VECTOR(Eig);
	FREE_MATRIX(VC);
	FREE_MATRIX(D);
	FREE_MATRIX(De);
	FREE_MATRIX(H);

	return;
}

void DeStandardErr(double **W, int *Xic1, int N, int p, double **DeSt, int K)
{

	double *new, **H, ***D;
	int i, j, h, k, m;

	m = p * (p+1) / 2;

	MAKE_VECTOR(new,p);

	MAKE_3ARRAY(D,K,p,p);
	MAKE_MATRIX(H,p,p);

	for (k=0;k<K;k++){
	
		i = 0;
		j = -1;

		for (h=0;h<m;h++){

			j++;

			if (j == p){
			
				i++;
				j = i;

			}

			D[k][i][j] = DeSt[k][h];

			if (i != j){
				
				D[k][j][i] = DeSt[k][h];

			}
		}

	}

	for (i=0;i<N;i++){

		vecxmat(W[i], p, D[Xic1[i]], p, p, new);

		for (j=0;j<p;j++){			

			W[i][j] = new[j];

		}
				
	}

	FREE_VECTOR(new);

	FREE_MATRIX(H);
	FREE_3ARRAY(D);

	return;
}




double TestKvsKstar_HIER(char *output_dir, int K, int Kstar, double **X, int N, int p, int RefDistN, int Nrepl, int iter, int hcrit)
{
/*
	Returns Pvalue of the test: K vs Kstar clusters
	Parameters:
		K - # of clusters under H0
		Kstar - # of parameters under Ha
		X - data
		N - # of observations
		p - dimension of data
		RefDistN - # of points for Reference distribution
		Nrepl - # if initial trials for kmeans
		iter - # of iterations for kmeans
*/

	FILE *fvc;
	char www[10];
	char fname[266];
	double **E;
	double **Y;
	double **W;
	double **c,**Xc;
	double **VarCov;
	double **DeSt;

	double *RefDist, *RefDist_;
	double *wss;
	double *u;

	int *ic1,*Xic1;
	int *nc, *Xnc;

	int *SampRes;
	
	double Wstat, Wstat_;
	double Pv, Pv_;
	double normW;
	
	int code;
	int i,j,k,sch;
	int m, randI;

	MAKE_MATRIX(E,N,p);
	MAKE_MATRIX(Y,N,p);
	MAKE_MATRIX(W,N,p);

	MAKE_VECTOR(u,N);
	MAKE_VECTOR(RefDist,RefDistN);
	MAKE_VECTOR(RefDist_,RefDistN);

	MAKE_VECTOR(ic1,N);
	MAKE_VECTOR(nc,Kstar);
	MAKE_VECTOR(Xnc,Kstar);
	MAKE_VECTOR(wss,Kstar);
	MAKE_MATRIX(c,Kstar,p);

	MAKE_VECTOR(Xic1,N);
	MAKE_MATRIX(Xc,Kstar,p);

	MAKE_VECTOR(SampRes, N);

	m = p*(p+1)/2;

	MAKE_MATRIX(DeSt,K,m);
	MAKE_MATRIX(VarCov,K,m);

	Pv = 0.0;
	Pv_ = 0.0;

/*		Find the test statistic		*/
	
	Wstat_ = HierClust(X,N,p,Xc,K,Xic1,Xnc,wss,hcrit);
	Wstat = HierClust(X,N,p,c,Kstar,ic1,nc,wss,hcrit);
	printf("WSS under H_0 = %f; ", Wstat_);
	Wstat_ -= Wstat;
	printf("WSS under H_a = %f; Statistic: %f\n", Wstat, Wstat_);
	
	sprintf(fname, "%s/", output_dir);	
	sprintf(www,"%i",Kstar);
	strcat(fname, "classes_");
	strcat(fname, www);

	fvc = fopen(fname,"w");

	for (i=0;i<N;i++){
	
		fprintf(fvc,"%d ", ic1[i]);

	}
	
	fprintf(fvc,"\n");

	fclose(fvc);	
    



/*		Find the reference distribution		*/


	SampleVarCov(X, N, p, K, Xc, Xic1, Xnc, VarCov);

	printf("Cluster sizes: ");
	for (k=0;k<K;k++){
		printf("%i ", Xnc[k]);
	}	
	printf("\n");

	

	for (i=0;i<N;i++){

		for (j=0;j<p;j++){
					
			E[i][j] = X[i][j] - Xc[Xic1[i]][j];

		}

	}


	StandardErr(E, Xic1, N, p, VarCov, K, DeSt);


	for (i=0;i<N;i++){

		u[i] = 0;

		for (j=0;j<p;j++){				
					
			u[i] = u[i] + E[i][j] * E[i][j];

		}

		u[i] = pow(u[i], 0.5);

	}
	

	for (sch=0;sch<RefDistN;sch++){
		
		code = srswor(N, N, SampRes);		
		if (code != 0) printf("Problem with function srswor\n");

		for (i=0;i<N;i++){

			randI = SampRes[i];
				
			normW = 0;
		
			for (j=0;j<p;j++){
					
				W[i][j] = rnorm(0.0, 1.0);
				normW = normW + W[i][j] * W[i][j];

			}

			normW = pow(normW, 0.5);
				
			for (j=0;j<p;j++){

				W[i][j] = W[i][j] / normW * u[randI];

			}

		}
		
		DeStandardErr(W, Xic1, N, p, DeSt, K);
		
		for (i=0;i<N;i++){
		
			if (K == 1){
			
				for (j=0;j<p;j++){
					
					Y[i][j] = W[i][j] + Xc[0][j];

				}

	
			} else {
				
				for (j=0;j<p;j++){
					
					Y[i][j] = W[i][j] + Xc[Xic1[i]][j];


				}

			}


		}
		
		RefDist[sch] = HierClust(Y, N, p, c, Kstar, ic1, nc, wss, hcrit);
		RefDist_[sch] = HierClust(Y, N, p, c, K, ic1, nc, wss, hcrit) - RefDist[sch];

		if (RefDist_[sch] > Wstat_)
			Pv_ = Pv_ + 1.0;
	}

	Pv = Pv_ / RefDistN;

	FREE_VECTOR(RefDist_);
	FREE_MATRIX(E);
	FREE_MATRIX(Y);
	FREE_MATRIX(W);

	FREE_VECTOR(u);
	FREE_VECTOR(RefDist);

	FREE_VECTOR(ic1);
	FREE_MATRIX(c);
	FREE_VECTOR(nc);
	FREE_VECTOR(Xnc);
	FREE_VECTOR(wss);

	FREE_VECTOR(Xic1);
	FREE_MATRIX(Xc);

	FREE_MATRIX(VarCov);
	FREE_MATRIX(DeSt);

	FREE_VECTOR(SampRes);

	return Pv;

}

double TestKvsKstar_SPH(char *output_dir, int K, int Kstar, double **X, int N, int p, int RefDistN, int Nrepl, int iter, int hcrit)
{

/*
	Returns Pvalue of the test: K vs Kstar clusters
	
	Parameters:
		K - # of clusters under H0
		Kstar - # of parameters under Ha
		X - data
		N - # of observations
		p - dimension of data
		RefDistN - # of points for Reference distribution
		Nrepl - # if initial trials for kmeans
		iter - # of iterations for kmeans
		hcrit - set to 1 and only used for kmeans under Ward's criterion
*/

	FILE *fvc;
	char www[10];
	char fname[266];
	double **E;
	double **Y;
	double **W;
	double **c,**Xc;
	double **VarCov;
	double **DeSt;

	double *RefDist, *RefDist_;
	double *wss;
	double *u;

	int *ic1,*Xic1;
	int *nc, *Xnc;
	int *SampRes;
	int ifault;

	double Wstat, Wstat_;
	int Pv = 0;
	double normW;
	
	int i,j,k,sch;
	int m, randI;

	hcrit = 1;
	MAKE_MATRIX(E,N,p);
	MAKE_MATRIX(Y,N,p);
	MAKE_MATRIX(W,N,p);

	MAKE_VECTOR(u,N);
	MAKE_VECTOR(RefDist,RefDistN);
	MAKE_VECTOR(RefDist_,RefDistN);

	MAKE_VECTOR(ic1,N);
	MAKE_VECTOR(nc,Kstar);
	MAKE_VECTOR(Xnc,Kstar);
	MAKE_VECTOR(wss,Kstar);
	MAKE_MATRIX(c,Kstar,p);

	MAKE_VECTOR(Xic1,N);
	MAKE_MATRIX(Xc,Kstar,p);
	
	MAKE_VECTOR(SampRes, N);


	m = p*(p+1)/2;

	MAKE_MATRIX(DeSt,K,m);
	MAKE_MATRIX(VarCov,K,m);

/*		Find the test statistic		*/
	
	/*
	  Wstat_ = RunKMeansNew(X,N,p,Xc,K,Xic1,Xnc,iter,wss,&ifault, hcrit, 1);
	  Wstat = RunKMeansNew(X,N,p,c,Kstar,ic1,nc,iter,wss,&ifault, hcrit, 1);
	*/
	
	Wstat_ = RunKMeans(X,N,p,Xc,K,Xic1,Xnc,iter,wss,&ifault, Nrepl, 1, Xc, 0);
	Wstat = RunKMeans(X,N,p,c,Kstar,ic1,nc,iter,wss,&ifault, Nrepl , 1, Xc, 0);
	printf("K under H0: %d; K under Ha: %d\n", K, Kstar);
	printf("WSS under H_0 = %f; ", Wstat_);
	Wstat_ -= Wstat;
	printf("WSS under H_a = %f; Statistic: %f\n", Wstat, Wstat_);

  	sprintf(fname, "%s", output_dir);	
	sprintf(www,"%i",Kstar);
	strcat(fname, "_");
	strcat(fname, www);

	sprintf(www,"%i",Kstar);
	  
	strcat(fname, "_");
	strcat(fname, www);

	fvc = fopen(fname,"w");

	for (i=0;i<N;i++){
	
		fprintf(fvc,"%d ", ic1[i]);

	}
	
	fprintf(fvc,"\n");

	fclose(fvc);	
    
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/*		Find the reference distribution		*/

	printf("Cluster sizes: ");
	for (k=0;k<K;k++)
		printf("%i ",Xnc[k]);
	printf("\n");

	for (i=0;i<N;i++)
		for (j=0;j<p;j++)
			E[i][j] = X[i][j] - Xc[Xic1[i]][j];

	for (i=0;i<N;i++) {
		u[i] = 0;
		for (j=0;j<p;j++)
			u[i] = u[i] + E[i][j] * E[i][j];
		u[i] = sqrt(u[i]);
	}

	for (sch=0;sch<RefDistN;sch++){
		srswor(N, N, SampRes);
		for (i=0;i<N;i++){
			randI = SampRes[i];
			normW = 0;
			for (j=0;j<p;j++){
				W[i][j] = rnorm(0.0, 1.0);
				normW = normW + W[i][j] * W[i][j];
			}
			normW = sqrt(normW);
			for (j=0;j<p;j++)
				W[i][j] = W[i][j] / normW * u[randI];
		}
		for (i=0;i<N;i++) {
			if (K == 1)
				for (j=0;j<p;j++)
					Y[i][j] = W[i][j] + Xc[0][j];
			else 
				for (j=0;j<p;j++)
					Y[i][j] = W[i][j] + Xc[Xic1[i]][j];
		}
		/*		RefDist[sch] = RunKMeansNew(Y, N, p, c, Kstar, ic1, nc, iter, wss, &ifault, hcrit, 0);
				RefDist_[sch] = RunKMeansNew(Y, N, p, c, K, ic1, nc, iter, wss, &ifault, hcrit, 0) - RefDist[sch]; */

		RefDist[sch] = RunKMeans(Y, N, p, c, Kstar, ic1, nc, iter, wss, &ifault, Nrepl, 0, Xc, 0);
		RefDist_[sch] = RunKMeans(Y, N, p, c, K, ic1, nc, iter, wss, &ifault, Nrepl, 0, Xc, 1) - RefDist[sch];
		
		/*		printf("Bootstrapped statistic = %f\n", RefDist_[sch]);*/

		if (RefDist_[sch] > Wstat_)
			Pv++;
	}
	FREE_VECTOR(RefDist_);
	FREE_VECTOR(Xnc);
	FREE_MATRIX(E);
	FREE_MATRIX(Y);
	FREE_MATRIX(W);

	FREE_VECTOR(u);
	FREE_VECTOR(RefDist);

	FREE_VECTOR(ic1);
	FREE_MATRIX(c);
	FREE_VECTOR(nc);
	FREE_VECTOR(wss);

	FREE_VECTOR(Xic1);
	FREE_MATRIX(Xc);

	FREE_MATRIX(VarCov);
	FREE_MATRIX(DeSt);
	
	FREE_VECTOR(SampRes);

	return (1.*Pv)/RefDistN;

}


int TestForward(char *output_dir, int minK, int maxK, double **X, int N, int p, int RefDistN, int Nrepl, int iter, double alpha, double *PVforw, int hcrit, short meth)
{

	int l, K, Kstar;

	K = minK;
	Kstar = minK+1;

	l = 0;

	while (Kstar <= maxK){
	
		if (meth == 1) {
			printf("testing with hierarchical clustering\n");
			PVforw[l] = TestKvsKstar_HIER(output_dir, K, Kstar, X, N, p, RefDistN, Nrepl, iter, hcrit);
		} else {
			printf("testing with random kmeans\n");
			PVforw[l] = TestKvsKstar_SPH(output_dir, K, Kstar, X, N, p, RefDistN, Nrepl, iter, hcrit);
		}	
	
		printf("%i vs %i pvalue = %f\n",K,Kstar,PVforw[l]);

		if (PVforw[l] < alpha){
		
			K = Kstar;
			Kstar = K + 1;
			printf("\n");

		} else {

			Kstar = Kstar + 1;

		}
		l++;
	}

	return K;

}

