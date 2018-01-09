#include <math.h>
#include <time.h>			
#include <stdio.h>
#include <stdlib.h>
#include <mex.h>
#include "matrix.h"

/* screening 
Author: Sen Yang
Mar 23, 2013
Arizona State University
*/
bool screening(double* v, const int K, const double lam, const double rho, const double EPS)
{
	double sumT = 0;
    double threshold = rho;
	for(int k = 0; k<K-1;k++)
	{
		sumT += v[k];
		threshold += lam;
		if (fabs(sumT) + EPS> threshold )
			return false;
	}
	sumT += v[K-1];
	threshold += lam - rho;
	if (fabs(sumT)  + EPS> threshold)
		return false;

	for(int k = 1; k < K-1; k++)
		for(int j = k; j<K-1;j++)
		{
			sumT = 0;
			threshold = (j-k+1)*lam + 2*rho;
			for(int i = k; i<j+1;i++)
				sumT += v[i];
			if (fabs(sumT) + EPS> threshold )
				return false;
		}
    
	sumT = 0;
	threshold = rho;
	for(int k = K-1; k>0;k--)
	{
		sumT += v[k];
		threshold += lam;
		if (fabs(sumT) + EPS> threshold )
			return false;
	}
	return true;

}

/* find freeset */
size_t findFreeSet(double* G, double* P, double* v,double* fx, double* fy, const int N, 
				   const int K, const double lam, const double rho, const double EPS)
{
	size_t indxNum = 0;
	bool flag = true;
	size_t indx,indxi,N2 = N*N;
	for(int i = 0; i<N; i++)
		for(int j = i+1; j<N;j++)
		{
			indx = i*N+j;
			flag = true;
			for(int k=0; k<K;k++)
			{
				indxi = indx + k*N2;
				if (P[indxi]!=0.0)
				{
					flag = false; break;
				}
				else
					v[k] = G[indxi];
			}

			if(flag && screening(v,K,lam,rho,EPS))
				continue;

			fx[indxNum] = (double) i;
			fy[indxNum] = (double) j;
			indxNum++;
		}

	return indxNum;
}

size_t globalScreening(double* G, double* adj, const int N, const int K, const double lam, const double rho)
{
	size_t indxNum = 0;
	size_t indx,N2 = N*N;
	double *v = new double[K];
	for(int i = 0; i<N; i++)
		for(int j = i+1; j<N;j++)
		{
			indx = i*N+j;
			for(int k=0; k<K;k++)
				v[k] = G[indx + k*N2];

			if(screening(v,K,lam,rho, 0.0))
				continue;

            adj[i*N+j] = adj[j*N+i] = 1;
			indxNum++;
		}
    delete [] v;
	return indxNum;
}


/* closed form solution for fused signal approximator with length 2*/
void fused_lasso2(double* v, double* z,double L, double R)
{
	if (z[0] > z[1] + 2*L)
	{
		v[0] = z[0] - L;
		v[1] = z[1] + L;
	}
	else if (z[1] > z[0] + 2*L)
	{
		v[0] = z[0] + L;
		v[1] = z[1] - L;
	}
	else
	{
		v[1] = v[0] = (z[1]+z[0])/2;
	}

	if (v[1]>R)
		v[1] -= R;
	else if(v[1]<-R)
		v[1] += R;
	else
		v[1] = 0;

	if (v[0]>R)
		v[0] -= R;
	else if(v[0]<-R)
		v[0] += R;
	else
		v[0] = 0;
}

/* fused_lasso is an efficient algorithm for fused signal approximator
Please refer to 

L. Condat, ��A direct algorithm for 1D total variation denoising,�� preprint hal-00675043, 2011. 

This function is written by Dr.Laurent Condat
*/
void fused_lasso(double* output, double* input, const size_t width, const double lambda, const double mu) {
	if (width>0) {				/*to avoid invalid memory access to input[0]*/
		int k=0, k0=0;			/*k: current sample location, k0: beginning of current segment*/
		double umin=lambda, umax=-lambda;	/*u is the dual variable*/
		double vmin=input[0]-lambda, vmax=input[0]+lambda;	/*bounds for the segment's value*/
		int kplus=0, kminus=0; 	/*last positions where umax=-lambda, umin=lambda, respectively*/
		const double twolambda=2.0*lambda;	/*auxiliary variable*/
		const double minlambda=-lambda;		/*auxiliary variable*/
		for (;;) {				/*simple loop, the exit test is inside*/
			while (k==width-1) {	/*we use the right boundary condition*/
				if (umin<0.0) {			/*vmin is too high -> negative jump necessary*/
					vmin = vmin>mu ? vmin-mu : vmin<-mu ? vmin+mu : 0.0;  
					do output[k0++]=vmin; while (k0<=kminus);
					umax=(vmin=input[kminus=k=k0])+(umin=lambda)-vmax;
				} else if (umax>0.0) {	/*vmax is too low -> positive jump necessary*/
					vmax = vmax>mu ? vmax-mu : vmax<-mu ? vmax+mu : 0.0;
					do output[k0++]=vmax; while (k0<=kplus);
					umin=(vmax=input[kplus=k=k0])+(umax=minlambda)-vmin;
				} else {
					vmin+=umin/(k-k0+1);
					vmin = vmin>mu ? vmin-mu : vmin<-mu ? vmin+mu : 0.0;  
					do output[k0++]=vmin; while(k0<=k); 
					return;
				}
			}

			if(k>width-1)
				break;

			if ((umin+=input[k+1]-vmin)<minlambda) {		/*negative jump necessary*/
				vmin = vmin>mu ? vmin-mu : vmin<-mu ? vmin+mu : 0.0; 
				do output[k0++]=vmin; while (k0<=kminus);
				vmax=(vmin=input[kplus=kminus=k=k0])+twolambda;
				umin=lambda; umax=minlambda;
			} else if ((umax+=input[k+1]-vmax)>lambda) {	/*positive jump necessary*/
				vmax = vmax>mu ? vmax-mu : vmax<-mu ? vmax+mu : 0.0;
				do output[k0++]=vmax; while (k0<=kplus);
				vmin=(vmax=input[kplus=kminus=k=k0])-twolambda;
				umin=lambda; umax=minlambda;
			} else { 	/*no jump necessary, we continue*/
				k++;
				if (umin>=lambda) {		/*update of vmin*/
					vmin+=(umin-lambda)/((kminus=k)-k0+1);
					umin=lambda;
				} 
				if (umax<=minlambda) {	/*update of vmax*/
					vmax+=(umax+lambda)/((kplus=k)-k0+1);
					umax=minlambda;
				} 	
			}
		}
	}
}

void fmgl_subfusedLasso(double* output, double* input,  const double* fx, const double* fy, double lam, double rho,
						const int N, const int K, const size_t Num)
{
	size_t N2 = N*N;
	double* sarrayi = new double[K];
	double* sarrayo = new double[K];
	size_t idx1,idx2;
	for(size_t i = 0; i<Num; i++)
	{
		double temp = fx[i]*N+fy[i];
		idx1 = (size_t) temp;
		temp = fy[i]*N+fx[i];
		idx2 = (size_t) temp;
		for(int k = 0; k<K;k++)
			sarrayi[k] = input[idx1 + k*N2];
		if(K==2)
			fused_lasso2(sarrayo, sarrayi, rho, lam);
		else
		    fused_lasso(sarrayo, sarrayi, K, rho, lam);
		for(int k = 0; k<K;k++)
		{
			output[idx1 + k*N2] = sarrayo[k];
			output[idx2 + k*N2] = sarrayo[k];
		}
	}
    
    for (size_t k = 0; k<K;k++)
    {
        size_t kN2 = k*N2;
        for(size_t i = 0; i<N; i++)
            output[i*N+ i + kN2] = input[i*N+ i + kN2];
    }
    
	delete[] sarrayi;
	delete[] sarrayo;
}