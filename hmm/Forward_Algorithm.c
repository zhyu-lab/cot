/* 
 Forward HMM algorithm


 Usage  
 -------
 //forward only
 [alpha, scale, loglik] = Forward_Algorithm(PI, A, L) 

 Inputs
 -------

 PI            Initial proabilities (N x 1) : Pr(x_1 = i) , i=1,...,N

 A             State transition probabilities matrix Pr(x_{k} = i| x_{k - 1} = j) such
               sum_{x_k}(A) = 1 => sum(A , 2) = 1, sum of row equals 1

 L        Time indexed Likelihood matrix Pr(z_k | x_k = i) (N x K), i=1,...,N, k=1,...,K. 
               Extracted from B matrix such that B = Pr(z | x) (M x N), sum(B , 1) = 1 and B(z_k , :)' = L(: , k).


 Ouputs
 -------
 alpha         alpha(i,t) = p(Q(t)=i | y(1:t)) (or p(Q(t)=i, y(1:t)) if scaled=0)
 
 scale         
 
 
 to complie 
 mex -output Forward_Algorithm.dll Forward_Algorithm.c

*/

#include <math.h>
#include <stdio.h>
#include "mex.h"

void ForwardWithScale(int , int , double *, double *, 
	double *, double *, double *, double *);
void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )
{
	
	
	double *PI , *A, *L ;
	
	double *alpha, *scale, *loglik; 
	
	const int *dimsPI , *dimsA , *dimsL ;
	
	
	int K , N , numdimsPI , numdimsA  , numdimsL;
	
	/*---------------------------------------------------------------*/
	/*---------------------- PARSE INPUT ----------------------------*/	
	/*---------------------------------------------------------------*/
	
	if( (nrhs != 3))
		
	{
		
		mexErrMsgTxt("3 inputs are requiered");
		
	}
	
	PI         = mxGetPr(prhs[0]);
    
  numdimsPI  = mxGetNumberOfDimensions(prhs[0]);
    
	dimsPI     = mxGetDimensions(prhs[0]);
	
	if ( (numdimsPI>2) & (dimsPI[1] > dimsPI[0]) )
	{
		
		mexErrMsgTxt("PI must be (N x 1)");			 
		
	}
	
	
	A         = mxGetPr(prhs[1]);
    
  numdimsA  = mxGetNumberOfDimensions(prhs[1]);
    
	dimsA     = mxGetDimensions(prhs[1]);
	
	if ( (numdimsA>2) & (dimsA[1] != dimsA[0]) )
	{
		
		mexErrMsgTxt("A must be (N x N)");			 
		
	}
	
	
	L         = mxGetPr(prhs[2]);
    
  numdimsL  = mxGetNumberOfDimensions(prhs[2]);
    
	dimsL     = mxGetDimensions(prhs[2]);
	
	if ( (numdimsL>2) & (dimsL[0] != dimsA[0]) )
	{
		
		mexErrMsgTxt("L must be (N x K)");			 
		
	}
	
  N         = dimsL[0];
	
	K         = dimsL[1];
	
	
	/*---------------------------------------------------------------*/
	/*---------------------- PARSE OUTPUT ---------------------------*/	
	/*---------------------------------------------------------------*/
	
	plhs[0]       = mxCreateDoubleMatrix(N , K , mxREAL);
	
	alpha         = mxGetPr(plhs[0]);
	
	
	plhs[1]       = mxCreateDoubleMatrix(1 , K , mxREAL);
    
	scale         = mxGetPr(plhs[1]);
    
    plhs[2]       = mxCreateDoubleMatrix(1 , 1 , mxREAL);
    
	loglik        = mxGetPr(plhs[2]);

	
	ForwardWithScale(N, K, PI, A, L, alpha, loglik, scale);
	
}

	
void ForwardWithScale(int N, int K, double *PI, double *A, double *L, 
	double *alpha, double *loglik, double *scale)
{
	int	i, j; 	/* state indices */
	int	t;	/* time index */

	double sum;	/* partial sum */

	/* 1. Initialization */

	scale[0] = 0.0;	
	for (i = 0; i < N; i++) {
		alpha[i] = PI[i]*L[i];
		scale[0] += alpha[i];
	}
	for (i = 0; i <N; i++) 
		alpha[i] /= scale[0];
	
	/* 2. Induction */

	for (t = 1; t < K; t++) {
		scale[t] = 0.0;
		for (j = 0; j <N; j++) {
			sum = 0.0;
			for (i = 0; i <N; i++)
				sum += alpha[i+(t-1)*N]*A[i+j*N]; 
				/*sum += alpha[t][i]* (phmm->A[i][j]); 
				sum   += At[j + iN]*alpha[j + kN1];*/
			
			alpha[j+t*N] = sum*L[j+t*N];	
			/*alpha[t+1][j] = sum*(phmm->B[j][O[t+1]]);*/
			scale[t] += alpha[j+t*N];
			/*scale[t+1] += alpha[t+1][j];*/
		}
        /*printf("%lf,",scale[t]);*/
		for (j = 0; j <N; j++) 
			alpha[j+t*N] /= scale[t];
			/*alpha[t+1][j] /= scale[t+1]; */
	}
    /*printf("\n");*/

	/* 3. Termination */
	*loglik = 0.0;

	for (t = 0; t < K; t++)
		*loglik += log(scale[t]);
}

