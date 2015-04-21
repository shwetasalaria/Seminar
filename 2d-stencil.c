#include<stdio.h>
#include "mpi.h"
#include<stdlib.h>

double **matrix, **mat;
void print(int matdim, double **matrix);
void init(int matdim)
{
  int i,j;
  matrix = (double**)calloc(matdim, sizeof(double*));
  mat = (double**)calloc(matdim, sizeof(double*));
  for(i = 0; i < matdim; i++)
  {
    matrix[i] = (double*)calloc(matdim, sizeof(double));
    mat[i] = (double*)calloc(matdim, sizeof(double));
  }
  for( i = 0; i < matdim; i++)
  for( j = 0; j < matdim; j++)
  {
     mat[i][j] = 1;
     if ( (i == 0) || (i == (matdim - 1)) || (j == 0) || (j == (matdim - 1)))
     {
	matrix[i][j] = 1;
     }
     else
     {   
     	matrix[i][j] = drand48();
     }
  }
  printf("\nThe second matrix is:\n");
  print(matdim, mat);
}

void dealloc(int matdim)
{
 int i;
 for ( i= 0; i< matdim; i++)
 {
        free(matrix[i]);
        free(mat[i]);
 }
 free(matrix);
 free(mat);
}

void swap()                           
{                                                          
    double **temp;                                           
    temp  = matrix;                                            
    matrix  = mat;                                            
    mat  = temp;                                           
}  

void print(int matdim, double **matrix)
{
 int i,j;
 for(i= 0; i < matdim; i++)
 {
        for (j= 0;j < matdim; j++)
        {
                printf("%f ", matrix[i][j]);
        }
        printf("\n");
 }
}

void update(int iteration, int matdim)
{
 int i,j, iter;
 for (iter = 1; iter <= iteration; iter ++)
 {
  printf("\nITERATION : %d", iter);
  for (i = 1; i < matdim - 1; i++)
  {
	for ( j = 1; j < matdim - 1; j++)
        {
		mat[i][j] = (matrix[i][j-1] + matrix[i-1][j] + matrix[i][j+1] + matrix[i+1][j] + matrix[i][j])/5;
	}
  }
 printf("\nThe updated second matrix is: \n");
 print(matdim,mat);
 swap(matrix,mat);
 printf("\nAfter swap function(), the original matrix is:\n");
 print(matdim,matrix); 
 }
}

int main(int argc,char *argv[])
{
 int processes, matdim,my_rank, iteration;

 MPI_Init(&argc, &argv);
 MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
 MPI_Comm_size(MPI_COMM_WORLD, &processes);

 matdim = 5;
 iteration = 2;

 if((argc > 1) && (argv[1] != NULL))
 {
 	matdim = atoi(argv[1]);
 }
 if((argc > 2) && (argv[2] != NULL))
 {
 	iteration = atoi(argv[2]);
 }
 if ( processes > matdim)
 {
        printf("\nProcesses can't be more than rows in this implementation\n");
        exit(0);
 }
 printf("\nInitializing matrix...\n");
 init(matdim);
 printf("\nInitialization completed.\n");
 printf("\nPrinting original matrix...\n");
 print(matdim,matrix);
 printf("\nUpdating elements..\n");
 update(iteration, matdim);
 printf("\nPrinting updated matrix..\n");
 print(matdim,matrix);
 printf("\nDeallocating matrix...\n");
 dealloc(matdim);
 
 MPI_Finalize();
 return 0;
}
