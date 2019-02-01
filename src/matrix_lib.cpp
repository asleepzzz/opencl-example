//------------------------------------------------------------------------------
//
//  PROGRAM: Matrix library for the multiplication driver
//
//  PURPOSE: This is a simple set of functions to manipulate
//           matrices used with the multiplcation driver.
//
//  USAGE:   The matrices are square and the order is
//           set as a defined constant, ORDER.
//
//  HISTORY: Written by Tim Mattson, August 2010
//           Modified by Simon McIntosh-Smith, September 2011
//           Modified by Tom Deakin and Simon McIntosh-Smith, October 2012
//           Updated to C++ Wrapper v1.2.6 by Tom Deakin, August 2013
//           Modified to assume square matrices by Simon McIntosh-Smith, Sep 2014
//
//------------------------------------------------------------------------------

#include "matmul.hpp"

//------------------------------------------------------------------------------
//
//  Function to compute the matrix product (sequential algorithm, dot prod)
//
//------------------------------------------------------------------------------

void seq_mat_mul_sdot(int M,int P,int N, std::vector<float>& A, std::vector<float>& B, std::vector<float>& C)
{
    int i, j, k;
    float tmp;

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            for (k = 0; k < P; k++) {
               C[i*N+j] += A[i*P+k] * B[k*N+j];
            }
        }
    }
}

//------------------------------------------------------------------------------
//
//  Function to initialize the input matrices A and B
//
//------------------------------------------------------------------------------
void initmat(int M ,int P,int N, std::vector<float>& A, std::vector<float>& B, std::vector<float>& C, bool isFirst)
{
    int i, j;

    /* Initialize matrices */

        if (isFirst) {
	for (i = 0; i < M; i++)
		for (j = 0; j < P; j++)
			A[i*P+j] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

	for (i = 0; i < P; i++)
		for (j = 0; j < N; j++)
			B[i*N+j] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        }
	for (i = 0; i < M; i++)
		for (j = 0; j < N; j++)
			C[i*N+j] = 0.0f;
}

//------------------------------------------------------------------------------
//
//  Function to set a matrix to zero
//
//------------------------------------------------------------------------------
void zero_mat (int M,int N, std::vector<float>& C)
{
    int i, j;

	for (i = 0; i < M; i++)
		for (j = 0; j < N; j++)
			C[i*N+j] = 0.0f;
}

//------------------------------------------------------------------------------
//
//  Function to fill Btrans(N,N) with transpose of B(N,N)
//
//------------------------------------------------------------------------------
void trans(int N, std::vector<float>& B, std::vector<float>& Btrans)
{
    int i, j;

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
		    Btrans[j*N+i] = B[i*N+j];
}

//------------------------------------------------------------------------------
//
//  Function to compute errors of the product matrix
//
//------------------------------------------------------------------------------
float error(int M,int N, std::vector<float>& C, std::vector<float>& product)
{
   int i,j;
   float cval, errsq, err;
   errsq = 0.0f;

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            err = C[i*N+j] -product[i*N+j];
            errsq += err * err;
        }
    }
 
    if (std::isnan(errsq) || errsq > TOL)
        printf("\n Errors in multiplication: %f\n",errsq);
    else
        printf("\n ============ok============\n");
    return errsq;
}

//------------------------------------------------------------------------------
//
//  Function to analyze and output results
//
//------------------------------------------------------------------------------
void results(int M,int P,int N, std::vector<float>& C, double run_time)
{

    float mflops;
    float errsq;

    mflops = 2.0 * M * P * N/(1000000.0f * run_time);
    printf(" %.2f seconds at %.1f MFLOPS \n",  run_time,mflops);
}
