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
 //printf("\nMatrix deallocated.\n");
}

void swap(double **orig_matrix, double **temp_matrix )                           
{                                                          
    double *temp;                                           
    temp  = *temp_matrix;                                            
    *temp_matrix  = *orig_matrix;                                            
    *orig_matrix  = temp;                                           
} 

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


int main(int argc,char *argv[])
{
 int matdim, processes, my_rank, iteration, destination, rows, num_elements, index, tag;
 double *matrix = NULL;
 double *result = NULL;
 double *recv_rows = NULL, *extra_row_left_rank = NULL, *extra_row_right_rank = NULL;
 int i,j,k,iter;
 
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
 
 //Iteration code
 for ( iter = 1; iter <= iteration ; iter++)
{
  
 if (my_rank == 0)
	printf("\n\n\n####################### Start of iteration %d ###################################\n\n\n", iter);
 MPI_Scatter(matrix, num_elements, MPI_DOUBLE, recv_rows , num_elements , MPI_DOUBLE, 0, MPI_COMM_WORLD);
 //printf("\n Iteration %d :Process %d received elements: \n", iter, my_rank);
 //print_rows(num_elements, recv_rows, matdim);
 
 //Processes send elements to neighbouring processes
 extra_row_left_rank = malloc(matdim * sizeof(double));
 extra_row_right_rank = malloc(matdim * sizeof(double));
 index = num_elements - matdim;
 tag = 1;
 if (my_rank == 0)
 {
	MPI_Send(&recv_rows[index], matdim, MPI_DOUBLE, my_rank +1 , tag , MPI_COMM_WORLD);
        MPI_Recv(extra_row_right_rank, matdim, MPI_DOUBLE, my_rank + 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("\n Process %d furthur receives elements from its right rank:\n", my_rank);
        //print_rows(matdim, extra_row_right_rank, matdim); 
 }
 
 else if (my_rank == (processes -1))
 {
	MPI_Send(recv_rows, matdim, MPI_DOUBLE, my_rank - 1, tag, MPI_COMM_WORLD);
        MPI_Recv(extra_row_left_rank, matdim, MPI_DOUBLE, my_rank - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("\n Process %d furthur receives elements fom its left rank:\n", my_rank);
        //print_rows(matdim, extra_row_left_rank, matdim);
 }

 else
 {
	MPI_Send(recv_rows, matdim, MPI_DOUBLE, my_rank - 1, tag, MPI_COMM_WORLD);
 	MPI_Send(&recv_rows[index], matdim, MPI_DOUBLE, my_rank + 1, tag, MPI_COMM_WORLD);
 	MPI_Recv(extra_row_left_rank, matdim, MPI_DOUBLE, my_rank - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 	MPI_Recv(extra_row_right_rank, matdim, MPI_DOUBLE, my_rank + 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 
 	//printf("\n Process %d receives elements from its left rank:\n", my_rank);
	//print_rows(matdim, extra_row_left_rank, matdim);
 	//printf("\n Process %d receives elements from its right rank:\n", my_rank);
 	//print_rows(matdim, extra_row_right_rank, matdim);
 }
 
 //Stencil computation
 //Converting 1-D matrix into 2-D matrix for ease of computation
 int rows_2D_matrix = 0; 
 if ( my_rank == 0 || my_rank == (processes - 1))
	rows_2D_matrix = rows + 1;
 else
 	rows_2D_matrix = rows +2;

 double intermediate_matrix[rows_2D_matrix][matdim];
 if (my_rank == 0)
 {
 	for (i = 0; i< rows_2D_matrix - 1; i++)
        {
         	for( j = 0; j< matdim ; j++)
                {
                	intermediate_matrix[i][j] = recv_rows[i * matdim + j];
                }
        }
       for ( k = 0; k< matdim; k++)
       	intermediate_matrix[rows_2D_matrix - 1][k] = extra_row_right_rank[k];
       /*
       i = 0, j = 0;
       printf("\nPrinting 2D-matrix for rank %d\n", my_rank);
       for(i= 0; i < rows_2D_matrix; i++)
 	{
        	for (j= 0;j < matdim; j++)
        	{
                	printf("%f ", intermediate_matrix[i][j]);
        	}
        	printf("\n");
 	}
        */
 }
 else if (my_rank == (processes - 1))
 {
        for ( k = 0; k< matdim; k++)
        intermediate_matrix[0][k] = extra_row_left_rank[k];
	//print_rows(num_elements, recv_rows, matdim);
	for (i = 0; i< rows_2D_matrix; i++)
        {
                for( j = 0; j< matdim ; j++)
                {
                        intermediate_matrix[i+1][j] = recv_rows[i * matdim + j];
                }
        }
       /*
       printf("\nPrinting 2D-matrix for rank %d\n", my_rank);
       for(i= 0; i < rows_2D_matrix; i++)
        {
                for (j= 0;j < matdim; j++)
                {
                        printf("%f ", intermediate_matrix[i][j]);
                }
                printf("\n");
        }
        */
 }      
 
 else
 {
    for ( k = 0; k< matdim; k++)
        intermediate_matrix[0][k] = extra_row_left_rank[k];
    for (i = 0; i< rows_2D_matrix; i++)
        {
                for( j = 0; j< matdim ; j++)
                {
                        intermediate_matrix[i+1][j] = recv_rows[i * matdim + j];
                }
        }
    for ( k = 0; k< matdim; k++)
        intermediate_matrix[rows_2D_matrix - 1][k] = extra_row_right_rank[k];
    /*
    printf("\nPrinting 2D-matrix for rank %d\n", my_rank);
       for(i= 0; i < rows_2D_matrix; i++)
        {
                for (j= 0;j < matdim; j++)
                {
                        printf("%f ", intermediate_matrix[i][j]);
                }
                printf("\n");
        }
   */
 }


 //Stencil code
 double result_matrix[rows_2D_matrix][matdim];
 for (i = 0; i< rows_2D_matrix; i++)
     for( j = 0; j< matdim; j++)
         result_matrix[i][j] = intermediate_matrix[i][j];
 for (i = 1;  i < rows_2D_matrix - 1; i++)
  {
        for ( j = 1; j < matdim - 1; j++)
        {
                result_matrix[i][j] = (intermediate_matrix[i][j-1] + intermediate_matrix[i-1][j] + intermediate_matrix[i][j+1] + intermediate_matrix[i+1][j] + intermediate_matrix[i][j])/5;
               
        }
  }
 /*
 printf("\nPrinting 2D-matrix computed after stencil's iteration 1 for rank %d\n", my_rank);
 for(i= 0; i < rows_2D_matrix; i++)
  {
                for (j= 0;j < matdim; j++)
                {
                        printf("%f ", result_matrix[i][j]);
                }
                printf("\n");
  }
  */
 // Converting 2D array into 1D array
 int lower_bound_x = 1, upper_bound_x = rows_2D_matrix - 1;
 double * recv_updated_rows = malloc(num_elements * sizeof(double));
 if ( my_rank == 0)
  {
	lower_bound_x = 0;
  }
 else if ( my_rank == (processes -1))
  {
        lower_bound_x = 1;
        upper_bound_x = rows_2D_matrix;
  }
 else
  {
        upper_bound_x = rows_2D_matrix - 1;
  }
 
  k = 0;
  for ( i = lower_bound_x; i < upper_bound_x; i++)
  {
	for ( j = 0; j < matdim ; j++)
	{
        	recv_updated_rows[k] = result_matrix[i][j];
  		k++;
        }
        
  }
  /*
  printf("\n Iteration %d :Updated row of rank %d is \n", iter, my_rank);
  for ( k=0; k< num_elements; k++)
  printf(" %f ", recv_updated_rows[k]);      
  */

  if (my_rank == 0)
  	result = malloc(matdim * matdim * sizeof(double));

  MPI_Gather(recv_updated_rows, num_elements, MPI_DOUBLE, result, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if (my_rank == 0)
  {
	printf("\n\n Iteration %d: The resultant matrix at process 0 is : \n\n", iter);
        print(matdim, result);
	swap(&matrix, &result);
  	printf("\n\n Iteration %d: After swap, the original matrix is: \n\n", iter);
  	print(matdim, matrix);
 //printf("\nCalling MPI_Finalize()");
        //printf("\nAfter swap, the other matrix is :\n\n");
        //print(matdim, result);
  	dealloc(result);
  }
 //Deallocation
 dealloc(extra_row_left_rank);
 dealloc(extra_row_right_rank);
 dealloc(recv_updated_rows);
 if (my_rank == 0)
 	printf("\n\n\n########################### End of iteration %d #################################\n\n\n", iter);
 MPI_Barrier(MPI_COMM_WORLD);
}

 //printing final matrix
 if(my_rank == 0)
 {
 printf("\nAll the iterations are completed.\n");
 printf("\n\nThe final matrix is : \n\n");
 print(matdim, matrix);      
 dealloc(matrix);
 }
 dealloc(recv_rows);
 printf("\n Process %d calling MPI_Finalize()\n", my_rank);
 printf("\n Status of MPI_Finalize() is %d \n", MPI_Finalize());

 
 return 0;
}
