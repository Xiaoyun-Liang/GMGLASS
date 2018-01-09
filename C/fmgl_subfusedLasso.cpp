#include "fmgl.h"
/* global screening 
Author: Sen Yang
Mar 23 2013, Arizona State University*/
void mexFunction (int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ])
{
	const mwSize *sDims;
	double *input, *output, *fx, *fy;
	double lam, rho;
	size_t Num,M,N,K;
	sDims = mxGetDimensions(prhs[0]);

	output = mxGetPr(prhs[0]);
	input = mxGetPr(prhs[1]);
	fx = mxGetPr(prhs[2]);
	fy = mxGetPr(prhs[3]);

	lam = mxGetScalar(prhs[4]);
    rho = mxGetScalar(prhs[5]);
    Num = (size_t) mxGetScalar(prhs[6]);
	M =(size_t) sDims[0];
	N =(size_t) sDims[1];
	K = (size_t) sDims[2];
	fmgl_subfusedLasso(output, input,  fx, fy, lam, rho, sDims[0], sDims[2], Num);

}