#include<stdio.h>
#include "mpi.h"
#include<stdlib.h>


void print (int matdim, double * row)
{
    int i,nElements;
    nElements = matdim * matdim * matdim;
    for (i=0; i < nElements; i++)
    {
        printf("%f ", row[i]);                                                                                                                                           
    }
    printf("\n");
}




double *init(int matdim, double *matrix)
{
  int i,j,k;

  matrix = malloc(matdim * matdim * matdim * sizeof(double));

  for( i = 0; i < matdim; i++)
  {
     for ( j = 0; j< matdim; j++)
     {
        for( k =0; k< matdim; k++)
        {
                if ( (i == 0) || (i == (matdim - 1)) || (j == 0) || (j == (matdim - 1)) || ( k == 0) || ( k == ( matdim -1)))
                {
                        matrix[i * matdim *matdim + j * matdim + k] = 1;
                }
                else
                {
                        matrix[i * matdim * matdim + j * matdim + k] = drand48();

                }
        }
    }
  }
  return matrix;
}

void dealloc(double *matrix)
{
 free(matrix);
}

void print_rows (int nElements, double *row, int matdim)
{
    int i;
    int size_of_row = matdim * matdim;
    for (i=0; i < nElements; i++)
    {
        printf("%f ", row[i]);
        if ((i+1) % size_of_row  == 0)
                printf("\n");                                                                                                                                                       
    }
    printf("\n");
}

void swap(double **orig_matrix, double **temp_matrix )                           
{                                                          
    double *temp;                                           
    temp  = *temp_matrix;                                            
    *temp_matrix  = *orig_matrix;                                            
    *orig_matrix  = temp;                                           
} 

int main(int argc,char *argv[])
{
 int matdim, processes, my_rank, iteration, destination, rows, num_elements, index, tag, row_size;
 double *matrix = NULL;
 double *result = NULL;
 double *recv_rows = NULL, *extra_row_left_rank = NULL, *extra_row_right_rank = NULL;
 int i,j,k,iter,l;

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
 }

 //Process 0 scatters rows among other processes
 rows = matdim/processes;
 num_elements = rows * matdim * matdim;
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
 row_size = matdim * matdim;
 extra_row_left_rank = malloc(matdim * matdim * sizeof(double));
 extra_row_right_rank = malloc(matdim * matdim * sizeof(double));
 index = num_elements - row_size;
 tag = 1;
 if (my_rank == 0)
 {
        MPI_Send(&recv_rows[index], row_size, MPI_DOUBLE, my_rank +1 , tag , MPI_COMM_WORLD);
        MPI_Recv(extra_row_right_rank, row_size, MPI_DOUBLE, my_rank + 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("\n Process %d furthur receives elements from its right rank:\n", my_rank);
        //print_rows(row_size, extra_row_right_rank, matdim); 
 }
 else if (my_rank == (processes -1))
 {
        MPI_Send(recv_rows, row_size, MPI_DOUBLE, my_rank - 1, tag, MPI_COMM_WORLD);
        MPI_Recv(extra_row_left_rank, row_size, MPI_DOUBLE, my_rank - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("\n Process %d furthur receives elements fom its left rank:\n", my_rank);
        //print_rows(row_size, extra_row_left_rank, matdim);
 }
 else
 {
        MPI_Send(recv_rows, row_size, MPI_DOUBLE, my_rank - 1, tag, MPI_COMM_WORLD);
        MPI_Send(&recv_rows[index], row_size, MPI_DOUBLE, my_rank + 1, tag, MPI_COMM_WORLD);
        MPI_Recv(extra_row_left_rank, row_size, MPI_DOUBLE, my_rank - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(extra_row_right_rank, row_size, MPI_DOUBLE, my_rank + 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("\n Process %d receives elements from its left rank:\n", my_rank);
        //print_rows(row_size, extra_row_left_rank, matdim);
        //printf("\n Process %d receives elements from its right rank:\n", my_rank);
        //print_rows(row_size, extra_row_right_rank, matdim);
 }
 
 //Stencil computation
 //Converting 1-D matrix into 3-D matrix for ease of computation
 int rows_3D_matrix = 0; 
 if ( my_rank == 0 || my_rank == (processes - 1))
	rows_3D_matrix = rows + 1;
 else
 	rows_3D_matrix = rows +2;

 double intermediate_matrix[rows_3D_matrix][matdim][matdim];
 if (my_rank == 0)
 {
 	for (i = 0; i< rows_3D_matrix - 1; i++)
        {
         	for( j = 0; j< matdim ; j++)
                {
			for ( k = 0; k< matdim; k++)	
                		intermediate_matrix[i][j][k] = recv_rows[i * matdim * matdim + j * matdim + k];
                }
        }
       l = 0;
       for ( j = 0; j< matdim; j++)
	{
		for( k = 0; k< matdim; k++)
			{
       				intermediate_matrix[rows_3D_matrix - 1][j][k] = extra_row_right_rank[l];
                        	l++;
			}
        }
	/*
	printf("\nIteration %d : Printing 3D-matrix for rank %d \n", iter, my_rank);
        for(i= 0; i < rows_3D_matrix; i++)
        {
                for (j= 0;j < matdim; j++)
                {

                        for ( k = 0; k< matdim; k++)
                        {

                                printf("%f ", intermediate_matrix[i][j][k]);
                        }
                }
                printf("\n");
 
	}
	*/
 }

 else if (my_rank == (processes - 1))
 {
	l = 0;
        for ( j = 0; j< matdim; j++)
	{
		for( k = 0; k< matdim; k++)
		{
			intermediate_matrix[0][j][k] = extra_row_left_rank[l];
			l++;
		}
	}
	
	for (i = 0; i< rows_3D_matrix - 1 ;i++)
        {
                for( j = 0; j< matdim ; j++)
                {
                        for ( k = 0; k< matdim; k++)
				intermediate_matrix[i+1][j][k] = recv_rows[i * matdim * matdim + j * matdim + k];
                }
        }
	/*
	printf("\nIteration %d: Printing 3D-matrix for rank %d \n", iter, my_rank);
        for(i= 0; i < rows_3D_matrix; i++)
        {
                for (j= 0;j < matdim; j++)
                {

                        for ( k = 0; k< matdim; k++)
                        {

                                printf("%f ", intermediate_matrix[i][j][k]);
                        }
                }
                printf("\n");
 
	} 
	*/     
 }
 else
 {
    l = 0;
    for ( j = 0; j< matdim; j++)
	{
		for( k = 0; k< matdim; k++)
		{
			intermediate_matrix[0][j][k] = extra_row_left_rank[l];
			l++;
		}
	}
    for (i = 0; i< rows_3D_matrix - 2 ; i++)
        {
                for( j = 0; j< matdim ; j++)
                {
                        for( k = 0; k < matdim ; k++)
			intermediate_matrix[i+1][j][k] = recv_rows[i * matdim * matdim + j * matdim + k];
                }
        }
    l = 0;
    for ( j = 0; j< matdim; j++)
	{
		for ( k = 0; k< matdim; k++)
		{
			intermediate_matrix[rows_3D_matrix - 1][j][k] = extra_row_right_rank[l];
			l++;
		}
	}
	/*
        printf("\nIteration %d: Printing 3D-matrix for rank %d\n", iter, my_rank);
        for(i= 0; i < rows_3D_matrix; i++)
        {
                for (j= 0;j < matdim; j++)
                {

                        for ( k = 0; k< matdim; k++)
                        {

                                printf("%f ", intermediate_matrix[i][j][k]);
                        }
                }
                printf("\n");
	
 	}
	*/
 }

 //Stencil code
 double result_matrix[rows_3D_matrix][matdim][matdim];
 for (i = 0; i< rows_3D_matrix; i++)
     {
     	for( j = 0; j< matdim; j++)
	{
	  for( k = 0; k< matdim; k++)
	         result_matrix[i][j][k] = intermediate_matrix[i][j][k];
	}
      }	
 for (i = 1;  i < rows_3D_matrix - 1; i++)
  {
        for ( j = 1; j < matdim - 1; j++)
        {
	  for ( k = 1; k< matdim -1 ; k++)
                result_matrix[i][j][k] = (intermediate_matrix[i][j][k-1] + intermediate_matrix[i][j-1][k+1] + intermediate_matrix[i][j-1][k] + intermediate_matrix[i][j+1][k] + intermediate_matrix[i-1][j][k] + intermediate_matrix[i+1][j][k] + intermediate_matrix[i][j][k])/7;
               
        }
  }
 /*
 printf("\nIteration %d : Printing 3D-matrix computed after stencil's iteration 1 for rank %d\n", iter, my_rank);
 for(i= 0; i < rows_3D_matrix; i++)
  {
                for (j= 0;j < matdim; j++)
                {
                     
			for ( k = 0; k< matdim; k++)
			   printf("%f ", result_matrix[i][j][k]);
                }
                printf("\n");
  }
  */

 
 // Converting 3D array into 1D array
 int lower_bound_x = 1, upper_bound_x = rows_3D_matrix - 1;
 double * recv_updated_rows = malloc(num_elements * sizeof(double));
 if ( my_rank == 0)
  {
	lower_bound_x = 0;
  }
 else if ( my_rank == (processes -1))
  {
        lower_bound_x = 1;
        upper_bound_x = rows_3D_matrix;
  }
 else
  {
        upper_bound_x = rows_3D_matrix - 1;
  }
 
  l = 0;
  for ( i = lower_bound_x; i < upper_bound_x; i++)
  {
	for ( j = 0; j < matdim ; j++)
	{
		for ( k = 0; k< matdim; k++)
		{
        		recv_updated_rows[l] = result_matrix[i][j][k];
  			l++;
		}
        }
        
  }
  /*
  printf("\n Iteration %d :Updated row of rank %d is \n", iter, my_rank);
  for ( k=0; k< num_elements; k++)
  printf(" %f ", recv_updated_rows[k]);      
  */
  if (my_rank == 0)
  	result = malloc(matdim * matdim * matdim * sizeof(double));

  MPI_Gather(recv_updated_rows, num_elements, MPI_DOUBLE, result, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if (my_rank == 0)
  {
	//printf("\n\n Iteration %d: The resultant matrix at process 0 is : \n\n", iter);
        //print(matdim, result);
	swap(&matrix, &result);
  	printf("\n\n Iteration %d: After swap, the matrix is: \n\n", iter);
  	print(matdim, matrix);
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
 MPI_Finalize();

 
 return 0;
}
