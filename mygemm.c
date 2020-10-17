#include "mygemm.h"

/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/

//Register Reuse part 1
void dgemm0(const double* A, const double* B, double* C, const int n)
{
    int i;
    for(i=0; i<n; i++){
        int j;
        for(j=0; j<n; j++){
            int k;
            for(k=0; k<n; k++){
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }

}

void dgemm1(const double *A, const double *B, double *C, const int n) 
{
    int i;
    for (i = 0; i < n; i++)
    {
        int j;
        for (j = 0; j < n; j++)
        {
            register double R = C[i * n + j];
            int k;
            for (k = 0; k < n; k++)
            {
                R += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = R;
        }
    }    

}
//Register Reuse part 1 End

//Register Reuse part 2
void dgemm2(const double *A, const double *B, double *C, const int n) 
{
    int i;
    for(i = 0; i<n; i = i + 2){
        int j;
        for(j = 0; j<n; j = j + 2){
            register double c00 = C[i*n + j];
            register double c01 = C[i*n + j+1];
            register double c10 = C[(i+1)*n + j];
            register double c11 = C[(i+1)*n + j+1];

            int k;
            for(k = 0; k<n; k = k + 2){
                register double a00 = A[i*n + k];
                register double a01 = A[i*n + k+1];
                register double a10 = A[(i+1)*n + k];
                register double a11 = A[(i+1)*n + k+1];
                register double b00 = B[k*n + j];
                register double b01 = B[k*n + j+1];
                register double b10 = B[(k+1)*n + j];
                register double b11 = B[(k+1)*n + j+1];

                c00 += a00 * b00 + a01 * b10;
                c01 += a00 * b01 + a01 * b11;
                c10 += a10 * b00 + a11 * b10;
                c11 += a10 * b01 + a11 * b11;
            }
            C[i*n + j] = c00;
            C[i*n + j+1] = c01;
            C[(i+1)*n + j] = c10;
            C[(i+1)*n + j+1] = c11;
        }
    }
}
//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double *A, const double *B, double *C, const int n) 
{
    int i;
    for(i = 0; i<n; i = i + 3){
        int j;
        for(j = 0; j<n; j = j + 3){
            register double c00 = C[i*n + j];
            register double c01 = C[i*n + j+1];
            register double c02 = C[i*n + j+2];
            register double c10 = C[(i+1)*n + j];
            register double c11 = C[(i+1)*n + j+1];
            register double c12 = C[(i+1)*n + j+2];
            register double c20 = C[(i+2)*n + j];
            register double c21 = C[(i+2)*n + j+1];
            register double c22 = C[(i+2)*n + j+2];

            int k;
            for(k = 0; k<n; k++){
                register double a0 = A[i*n + k];
                register double a1 = A[(i+1)*n + k];
                register double a2 = A[(i+2)*n + k];
                register double b0 = B[k*n + j];
                register double b1 = B[k*n + j+1];
                register double b2 = B[k*n + j+2];

                c00 += a0 * b0;
                c01 += a0 * b1;
                c02 += a0 * b2;

                c10 += a1 * b0;
                c11 += a1 * b1;
                c12 += a1 * b2;

                c20 += a2 * b0;
                c21 += a2 * b1;
                c22 += a2 * b2;
            }
            C[i*n + j] = c00;
            C[i*n + j+1] = c01;
            C[i*n + j+2] = c02;

            C[(i+1)*n + j] = c10;
            C[(i+1)*n + j+1] = c11;
            C[(i+1)*n + j+2] = c12;

            C[(i+2)*n + j] = c20;
            C[(i+2)*n + j+1] = c21;
            C[(i+2)*n + j+2] = c22;
        }
    }
}
//Register Reuse part 3 End

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n) 
{
    int i;
    for (i = 0; i < n; i++)
    {
        int j;
        for (j = 0; j < n; j++)
        {
            register double r = C[i * n + j];
            int k;
            for (k = 0; k < n; k++)
            {
                r += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = r;
        }
    }    

}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i;
    for (i = 0; i < n; i += b)
    {
        int j;
        for (j = 0; j < n; j += b)
        {
            int k;
            for (k = 0; k < n; k += b)
            {
                int i1 ;
                for (i1 = i; i1 < i + b && i1 < n; i1++)
                {
                    int j1;
                    for (j1 = j; j1 < j + b && j1 < n; j1++)
                    {
                        register double r = C[i1 * n + j1];
                        int k1 = 0;
                        for (k1 = k; k1 < k + b && k1 < n; k1++)
                        {
                            r += A[i1 * n + k1] * B[k1 * n + j1];
                        }
                        C[i1 * n + j1] = r;
                    }
                }
            }
        }
    }    

}

void jik(const double *A, const double *B, double *C, const int n) 
{
    int j;
    for (j = 0; j < n; j++)
    {
        int i;
        for (i = 0; i < n; i++)
        {
            register double r = C[i * n + j];
            int k;
            for (k = 0; k < n; k++)
            {
                r += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = r;
        }
    }    

}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{
    int j;
    for (j = 0; j < n; j += b)
    {
        int i;
        for (i = 0; i < n; i += b)
        {
            int k;
            for (k = 0; k < n; k += b)
            {
                int j1;
                for (j1 = j; j1 < j + b && j1 < n; j1++)
                {
                    int i1;
                    for (i1 = i; i1 < i + b && i1 < n; i1++)
                    {
                        register double r = C[i1 * n + j1];
                        int k1;
                        for (k1 = k; k1 < k + b && k1 < n; k1++)
                        {
                            r += A[i1 * n + k1] * B[k1 * n + j1];
                        }
                        C[i1 * n + j1] = r;
                    }
                }
            }
        }
    }    

}

void kij(const double *A, const double *B, double *C, const int n) 
{
    int k;
    for (k = 0; k < n; k++)
    {
        int i;
        for (i = 0; i < n; i++)
        {
            register double r = A[i * n + k];
            int j;
            for (j = 0; j < n; j++)
            {
                C[i * n + j] += r * B[k * n + j];
            }
        }
    }

}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{
    int k;
    for (k = 0; k < n; k += b)
    {
        int i;
        for (i = 0; i < n; i += b)
        {
            int j;
            for (j = 0; j < n; j += b)
            {
                int k1;
                for (k1 = k; k1 < k + b && k1 < n; k1++)
                {
                    int i1;
                    for (i1 = i; i1 < i + b && i1 < n; i1++)
                    {
                        register double r = A[i1 * n + k1];
                        int j1;
                        for (j1 = j; j1 < j + b && j1 < n; j1++)
                        {
                            C[i1 * n + j1] += r * B[k1 * n + j1];
                        }
                    }
                }
            }
        }
    }

}


void ikj(const double *A, const double *B, double *C, const int n) 
{
    int i;
    for (i = 0; i < n; i++)
    {
        int k;
        for (k = 0; k < n; k++)
        {
            register double r = A[i * n + k];
            int j;
            for (j = 0; j < n; j++)
            {
                C[i * n + j] += r * B[k * n + j];
            }
        }
    }

}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i;
    for (i = 0; i < n; i += b)
    {
        int k;
        for (k = 0; k < n; k += b)
        {
            int j;
            for (j = 0; j < n; j += b)
            {
                int i1;
                for (i1 = i; i1 < i + b && i1 < n; i1++)
                {
                    int k1;
                    for (k1 = k; k1 < k + b && k1 < n; k1++)
                    {
                        register double r = A[i1 * n + k1];
                        int j1;
                        for (j1 = j; j1 < j + b && j1 < n; j1++)
                        {
                            C[i1 * n + j1] += r * B[k1 * n + j1];
                        }
                    }
                }
            }
        }
    }

}

void jki(const double *A, const double *B, double *C, const int n) 
{
    int j;
    for (j = 0; j < n; j++)
    {
        int k;
        for (k = 0; k < n; k++)
        {
            register double r = B[k * n + j];
            int i;
            for (i = 0; i < n; i++)
            {
                C[i * n + j] += A[i * n + k] * r;
            }
        }
    }

}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{
    int j;
    for (j = 0; j < n; j += b)
    {
        int k;
        for (k = 0; k < n; k += b)
        {
            int i;
            for (i = 0; i < n; i += b)
            {
                int j1;
                for (j1 = j; j1 < j + b && j1 < n; j1++)
                {
                    int k1;
                    for (k1 = k; k1 < k + b && k1 < n; k1++)
                    {
                        register double r = B[k1 * n + j1];
                        int i1;
                        for (i1 = i; i1 < i + b && i1 < n; i1++)
                        {
                            C[i1 * n + j1] += A[i1 * n + k1] * r;
                        }
                    }
                }
            }
        }
    }

}

void kji(const double *A, const double *B, double *C, const int n) 
{
    int k;
    for (k = 0; k < n; k++)
    {
        int j;
        for (j = 0; j < n; j++)
        {
            register double r = B[k * n + j];
            int i;
            for (i = 0; i < n; i++)
            {
                C[i * n + j] += A[i * n + k] * r;
            }
        }
    }

}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{
    int k;
    for (k = 0; k < n; k += b)
    {
        int j;
        for (j = 0; j < n; j += b)
        {
            int i;
            for (i = 0; i < n; i += b)
            {
                int k1;
                for (k1 = k; k1 < k + b && k1 < n; k1++)
                {
                    int j1;
                    for (j1 = j; j1 < j + b && j1 < n; j1++)
                    {
                        register double r = B[k1 * n + j1];
                        int i1;
                        for (i1 = i; i1 < i + b && i1 < n; i1++)
                        {
                            C[i1 * n + j1] += A[i1 * n + k1] * r;
                        }
                    }
                }
            }
        }
    }

}
//Cache Reuse part 3 End 

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{
   for (int i = 0; i < n; i += b) {
        for (int j = 0; j < n; j += b) {
            for (int k = 0; k < n; k += b) {
                int i1 = i, j1 = j, k1 = k;
                int ni = i + b > n ? n : i + b;
                int nj = j + b > n ? n : j + b;
                int nk = k + b > n ? n : k + b;

                for (i1 = i; i1 < ni; i1 += 3) {
                    for (j1 = j; j1 < nj; j1 += 3) {
                        int t = i1 * n + j1;
                        int tt = t + n;
                        int ttt = tt + n;
                        register double c00 = C[t];
                        register double c01 = C[t + 1];
                        register double c02 = C[t + 2];
                        register double c10 = C[tt];
                        register double c11 = C[tt + 1];
                        register double c12 = C[tt + 2];
                        register double c20 = C[ttt];
                        register double c21 = C[ttt + 1];
                        register double c22 = C[ttt + 2];

                        for (k1 = k; k1 < nk; k1 += 3) {
                                int ta = i1 * n + k1;
                                int tta = ta + n;
                                int ttta = tta + n;
                                int tb = k1 * n + j1 * n;
                                register double a0 = A[ta];
                                register double a1 = A[tta];
                                register double a2 = A[ttta];
                                register double b0 = B[tb];
                                register double b1 = B[tb + 1];
                                register double b2 = B[tb + 2];

                                c00 += a0 * b0;
                                c01 += a0 * b1;
                                c02 += a0 * b2;
                                c10 += a1 * b0;
                                c11 += a1 * b1;
                                c12 += a1 * b2;
                                c20 += a2 * b0;
                                c21 += a2 * b1;
                                c22 += a2 * b2;
                        }
                        C[t] = c00;
                        C[t + 1] = c01;
                        C[t + 2] = c02;
                        C[tt] = c10;
                        C[tt + 1] = c11;
                        C[tt + 2] = c12;
                        C[ttt] = c20;
                        C[ttt + 1] = c21;
                        C[ttt + 2] = c22;

                    }
                }
            }
        }
}