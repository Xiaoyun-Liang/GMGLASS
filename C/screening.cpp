#include "fmgl.h"
/* global screening 
Author: Sen Yang
Mar 23 2013, Arizona State University*/
void mexFunction (int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ])
{
	const mwSize *sDims;
	double *G, *adj;
	double lam, rho;
	size_t Num,M,N,K;
	sDims = mxGetDimensions(prhs[0]);

	G = mxGetPr(prhs[0]);

	lam = mxGetScalar(prhs[1]);
    rho = mxGetScalar(prhs[2]);
	
	M =(size_t) sDims[0];
	N =(size_t) sDims[1];
	K = (size_t) sDims[2];
	plhs[0] = mxCreateDoubleMatrix(M, M, mxREAL);
	adj = mxGetPr(plhs[0]);
	globalScreening(G, adj, (int)M, (int)K, lam, rho);
}
