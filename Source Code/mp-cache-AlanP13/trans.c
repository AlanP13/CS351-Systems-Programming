/* 
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"
#include <assert.h>
#ifdef DEBUG
#define ASSERT(COND) assert(COND)
#define REQUIRES(COND) assert(COND)
#define ENSURES(COND) assert(COND)
#else
#define ASSERT(COND) ((void)0)
#define REQUIRES(COND) ((void)0)
#define ENSURES(COND) ((void)0)
#endif
int is_transpose(int M, int N, int A[N][M], int B[M][N]);
void transpose_func64(int M, int N, int A[N][M], int B[M][N]);
void transpose_func32(int M, int N, int A[N][M], int B[M][N]);
void transpose_func6167(int M, int N, int A[N][M], int B[M][N]);
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
  if(M==32)
    {
      transpose_func32(M,N,A,B);
    }
  if(M==64)
    {
      transpose_func64(M,N,A,B);
    }
  if(M==61)
    {
      transpose_func6167(M,N,A,B);
    }
}
void transpose_func32(int M, int N, int A[M][N], int B[N][M])
{
  int i, j;
  int p,q;
  int temp;
  for(i=0; i<M; i+=8)
    {
      for(j=0; j<N; j+=8)
        {
	  for(p=i+7; (p>(i+3))&&p<M; p--)
            {
	      for(q=j+3; (q>=j)&&q<N; q--)
                {
		  B[q][p]=A[p][q];
                }
            }
	  for(p=i+7; (p>(i+3))&&p<M; p--)
	    {
	      for(q=j+4; (q<=(j+7))&&q<N; q++)
		{
		  B[q-4][p-4]=A[p][q];
		}
	    }
	  for(p=i; (p<=(i+3))&&p<M; p++)
            {
	      for(q=j+4; (q<=(j+7))&&q<N; q++)
                {
		  B[q][p]=A[p][q];
                }
            }
	  for(p=i; (p<=(i+3))&&p<M; p++)
            {
	      for(q=j; (q<=(j+3))&&q<N; q++)
                {
		  B[q+4][p+4]=A[p][q];
                }
            }
	  for(p=i; (p<=(i+3))&&p<M; p++)
            {
	      for(q=j; (q<=(j+3))&&q<N; q++)
                {
		  temp=B[p+4+(j-i)][q+4+(i-j)];
		  B[p+4+(j-i)][q+4+(i-j)]=B[p+(j-i)][q+(i-j)];
		  B[p+(j-i)][q+(i-j)]=temp;
                }
            }
        }
    }
}
void transpose_func64(int M, int N, int A[N][M], int B[N][M])
{
  int i,j,p,q;
  int xDash,yDash; 
  int temp;
  for(i=0; i<M; i+=8)
    {
      for(j=0; j<N; j+=8)
        {
	  xDash=j-i;
	  yDash=i-j;
	  for(p=i; p<(i+2); p++)
            {
	      for(q=j; q<(j+8); q++)
                {
		  B[p+2+xDash][q+yDash]=A[p][q];
                }
            }
	  for(p=i+2; p<(i+4); p++)
            {
	      for(q=j; q<(j+8); q++)
                {
		  B[p-2+xDash][q+yDash]=A[p][q];
                }
            }
	  for(p=i; p<(i+2); p++)
            {
	      for(q=j; q<(j+8); q++)
                {
		  temp=B[p+xDash][q+yDash];
		  B[p+xDash][q+yDash]=B[p+2+xDash][q+yDash];
		  B[p+2+xDash][q+yDash]=temp;
                }
            }
	  for(p=i+xDash; p<(i+4+xDash); p++)
            {
	      for(q=j+yDash; q<(j+4+yDash); q++)
                {
		  if((q-(j+yDash))<(p-(i+xDash)))
                    {
		      temp=B[p][q];
		      B[p][q]=B[q+xDash][p+yDash];
		      B[q+xDash][p+yDash]=temp;
                    }
                }
            }
	  for(p=i+xDash; p<(i+4+xDash); p++)
            {
	      for(q=j+4+yDash; q<(j+8+yDash); q++)
                {
		  if((q-(j+4+yDash))<(p-(i+xDash)))
                    {
		      temp=B[p][q];
		      B[p][q]=B[q+xDash-4][p+yDash+4];
		      B[q+xDash-4][p+yDash+4]=temp;
                    }
                }
            }
	  for(p=i; p<(i+2); p++)
            {
	      for(q=j+4; q<(j+8); q++)
                {
		  temp=B[p+xDash][q+yDash];
		  B[p+xDash][q+yDash]=B[p+2+xDash][q+yDash];
		  B[p+2+xDash][q+yDash]=temp;
                }
            }
	  for(p=i+4; p<(i+6); p++)
            {
	      for(q=j; q<(j+8); q++)
                {
		  B[p+2+xDash][q+yDash]=A[p][q];
                }
            }
	  for(p=i+6; p<(i+8); p++)
            {
	      for(q=j; q<(j+8); q++)
                {
		  B[p-2+xDash][q+yDash]=A[p][q];
                }
            }
	  for(p=i+4; p<(i+6); p++)
            {
	      for(q=j; q<(j+8); q++)
                {
		  temp=B[p+xDash][q+yDash];
		  B[p+xDash][q+yDash]=B[p+2+xDash][q+yDash];
		  B[p+2+xDash][q+yDash]=temp;
                }
            }
	  for(p=i+4+xDash; p<(i+8+xDash); p++)
            {
	      for(q=j+yDash; q<(j+4+yDash); q++)
                {
		  if((q-(j+yDash))<(p-(i+4+xDash)))
                    {
		      temp=B[p][q];
		      B[p][q]=B[q+xDash+4][p+yDash-4];
		      B[q+xDash+4][p+yDash-4]=temp;
                    }
                }
            }
	  for(p=i+4+xDash; p<(i+8+xDash); p++)
            {
	      for(q=j+4+yDash; q<(j+8+yDash); q++)
                {
		  if((q-(j+4+yDash))<(p-(i+4+xDash)))
                    {
		      temp=B[p][q];
		      B[p][q]=B[q+xDash][p+yDash];
		      B[q+xDash][p+yDash]=temp;
                    }
                }
            }
	  for(p=i+4; p<(i+6); p++)
            {
	      for(q=j+0; q<(j+4); q++)
                {
		  temp=B[p+xDash][q+yDash];
		  B[p+xDash][q+yDash]=B[p+2+xDash][q+yDash];
		  B[p+2+xDash][q+yDash]=temp;
                }
            }
	  for(p=i+4; p<(i+6); p++)
            {
	      for(q=j+0; q<(j+4); q++)
                {
		  temp=B[p+xDash][q+yDash];
		  B[p+xDash][q+yDash]=B[p-2+xDash][q+4+yDash];
		  B[p-2+xDash][q+4+yDash]=temp;
                }
            }
	  for(p=i+6; p<(i+8); p++)
            {
	      for(q=j+0; q<(j+4); q++)
                {
		  temp=B[p+xDash][q+yDash];
		  B[p+xDash][q+yDash]=B[p-6+xDash][q+4+yDash];
		  B[p-6+xDash][q+4+yDash]=temp;
                }
            }
        }
    }
}
void transpose_func6167(int M, int N, int A[N][M], int B[M][N])
{
  int i,j,p,q;
  int temp;
  for(i=0; i<M; i+=18)
    {
      for(j=0; j<N; j+=18)
        {
	  for(p=i; (p<(i+18)) && (p<M); p++)
            {
	      for(q=j; (q<(j+18)) && q<N; q++)
                {
		  if((p)==(q))
                    {
		      temp=A[p][p];
                    }
		  else
                    {
		      B[p][q]=A[q][p];
                    }
                }
	      if(i==j)
                {
		  B[p][p]=temp;
                }
            }
        }
    }
}
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
  int i,j,tmp;
  REQUIRES(M > 0);
  REQUIRES(N > 0);
  for (i = 0; i < N; i++)
    {
      for (j = 0; j < M; j++)
        {
	  tmp = A[i][j];
	  B[j][i] = tmp;
        }
    }
  ENSURES(is_transpose(M, N, A, B));
}
void registerFunctions()
{
  registerTransFunction(transpose_submit, transpose_submit_desc);
  registerTransFunction(trans, trans_desc);
}
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
  int i, j;
  for (i = 0; i < N; i++)
    {
      for (j = 0; j < M; ++j)
        {
	  if (A[i][j] != B[j][i])
            {
	      return 0;
            }
        }
    }
  return 1;
}
