#include<stdio.h>
#include "mpi.h"
#include<stdlib.h>

void print(int matdim, double *matrix);

double *init(int matdim, double *matrix)
{
  int i,j;

  matrix = malloc(matdim * matdim * sizeof(double));
  
  for( i = 0; i < matdim; i++)
  {
     for ( j = 0; j< matdim; j++)
     {
    	if ( (i == 0) || (i == (matdim - 1)) || (j == 0) || (j == (matdim - 1)))
     	{
		matrix[i * matdim + j] = 1;
    	 }
     	else
     	{   
     		matrix[i * matdim + j] = drand48();
     	}
    }
  }
  printf("\nThe initialized matrix is:\n");
  print(matdim, matrix);
  return matrix; 
}

void dealloc(double *matrix)
{
 free(matrix);
 printf("\nMatrix deallocated.\n");
}

/*void swap()                           
{                                                          
    double **temp;                                           
    temp  = matrix;                                            
    matrix  = mat;                                            
    mat  = temp;                                           
} */ 

void print(int matdim, double *matrix)
{
 int i,j;
 for(i= 0; i < matdim; i++)
 {
        for (j= 0;j < matdim; j++)
        {
                printf("%f ", matrix[i * matdim + j]);
        }
        printf("\n");
 }
}


void print_rows (int nElements, double *row, int matdim) 
{
    int i;
    for (i=0; i < nElements; i++) 
    {
        printf("%f ", row[i]);
        if ((i+1) % matdim  == 0)
        	printf("\n");                                                                                                                                                                 
    }
    printf("\n");
}
/*
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
 i}
 printf("\nThe updated second matrix is: \n");
 print(matdim,mat);
 swap(matrix,mat);
 printf("\nAfter swap function(), the original matrix is:\n");
 print(matdim,matrix); 
 }
}
*/

int main(int argc,char *argv[])
{
 int matdim, processes, my_rank, iteration, destination, rows, num_elements, index, tag;
 double *matrix = NULL;
 double *recv_rows = NULL, *extra_row_left_rank = NULL, *extra_row_right_rank = NULL;
 
 MPI_Status status;
 MPI_Init(&argc, &argv);
 MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
 MPI_Comm_size(MPI_COMM_WORLD, &processes);
 
 matdim = 4;
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
    if(my_rank == 0)
    {
    	printf("\nProcesses can't be more than rows in this implementation.\n");
    }
    MPI_Finalize();
    return 0;
 }

 if(my_rank == 0)
 { 
 	printf("\nInitializing matrix on process 0...\n");
 	matrix = init(matdim, matrix);
 	printf("\nInitialization completed.\n");

 	printf("\nPrinting original matrix...\n");
 	print(matdim, matrix);

        //printf("\nDeallocating matrix...\n");
        //dealloc(matrix);
 } 
 
 //Process 0 scatters rows among other processes
 rows = matdim/processes;
 num_elements = rows * matdim;
 recv_rows = malloc(num_elements * sizeof(double));
 MPI_Scatter(matrix, num_elements, MPI_DOUBLE, recv_rows , num_elements , MPI_DOUBLE, 0, MPI_COMM_WORLD);
 printf("\n Process %d received elements: \n", my_rank);
 print_rows(num_elements, recv_rows, matdim);
 
 //Processes send elements to neighbouring processes
 extra_row_left_rank = malloc(matdim * sizeof(double));
 extra_row_right_rank = malloc(matdim * sizeof(double));
 index = num_elements - matdim;
 tag = 1;
 if (my_rank == 0)
 {
	MPI_Send(&recv_rows[index], matdim, MPI_DOUBLE, my_rank +1 , tag , MPI_COMM_WORLD);
        MPI_Recv(extra_row_right_rank, matdim, MPI_DOUBLE, my_rank + 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("\n Process 0 furthur receives elements:\n");
        print_rows(matdim, extra_row_right_rank, matdim); 
}
 
 else if (my_rank == (processes -1))
 {
	MPI_Send(recv_rows, matdim, MPI_DOUBLE, my_rank - 1, tag, MPI_COMM_WORLD);
        MPI_Recv(extra_row_left_rank, matdim, MPI_DOUBLE, my_rank - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("\n Process %d furthur receives elements:\n", my_rank);
        print_rows(matdim, extra_row_left_rank, matdim);
 }
 else
 {
	MPI_Send(recv_rows, matdim, MPI_DOUBLE, my_rank - 1, tag, MPI_COMM_WORLD);
 	MPI_Send(&recv_rows[index], matdim, MPI_DOUBLE, my_rank + 1, tag, MPI_COMM_WORLD);
 	MPI_Recv(extra_row_left_rank, matdim, MPI_DOUBLE, my_rank - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 	MPI_Recv(extra_row_right_rank, matdim, MPI_DOUBLE, my_rank + 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 
 	printf("\n Process %d receives elements from its left rank:\n", my_rank);
	print_rows(matdim, extra_row_left_rank, matdim);
 	printf("\n Process %d receives elements from its right rank:\n", my_rank);
 	print_rows(matdim, extra_row_right_rank, matdim);
 }
 
 // Stencil computation
 


 dealloc(extra_row_left_rank);
 dealloc(extra_row_right_rank);
 dealloc(recv_rows);
 if (my_rank == 0)
 { 
 	printf("\nDeallocating matrix...\n");
 	dealloc(matrix);
 }
 

 /*
 MPI_Bcast(matrix, matdim * matdim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 printf("\nUpdating elements..\n");
 update(iteration, matdim);
 printf("\nPrinting updated matrix..\n");
 print(matdim,matrix);
 printf("\nDeallocating matrix...\n");
 dealloc();
 printf("\n Value of matdim %d\n", matdim);
*/

 MPI_Finalize();
 
 return 0;
}
