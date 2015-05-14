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

	matrix = (double*)malloc(matdim * matdim * matdim * sizeof(double));

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
	double start, end;

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
		matrix = init(matdim, matrix);
	}

	rows = matdim/processes;
	num_elements = rows * matdim * matdim;
	recv_rows = (double*)malloc(num_elements * sizeof(double));

	start = MPI_Wtime();
	for ( iter = 1; iter <= iteration ; iter++)
	{

		MPI_Scatter(matrix, num_elements, MPI_DOUBLE, recv_rows , num_elements , MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		row_size = matdim * matdim;
		extra_row_left_rank = (double*)malloc(matdim * matdim * sizeof(double));
		extra_row_right_rank = (double*)malloc(matdim * matdim * sizeof(double));
		index = num_elements - row_size;
		tag = 1;
		if (my_rank == 0)
		{
			printf("\nProcess 0 send and receive\n");
			MPI_Sendrecv(&recv_rows[index], row_size, MPI_DOUBLE, my_rank + 1, tag, extra_row_right_rank, row_size, MPI_DOUBLE, my_rank + 1, tag, MPI_COMM_WORLD, &status);
		}
		else if (my_rank == (processes -1))
		{
			printf("\n Last process send and receive \n");
			MPI_Sendrecv(recv_rows, row_size, MPI_DOUBLE, my_rank - 1, tag, extra_row_left_rank, row_size, MPI_DOUBLE, my_rank - 1, tag, MPI_COMM_WORLD, &status);
		}
		else
		{
		//double intermediate_matrix[rows_3D_matrix][matdim][matdim];
			printf("\n Process %d send and receive \n", my_rank);
			MPI_Sendrecv(recv_rows, row_size, MPI_DOUBLE, my_rank - 1, tag, extra_row_left_rank, row_size, MPI_DOUBLE, my_rank - 1, tag, MPI_COMM_WORLD, &status);
			MPI_Sendrecv(&recv_rows[index], row_size, MPI_DOUBLE, my_rank + 1, tag, extra_row_right_rank, row_size, MPI_DOUBLE, my_rank + 1, tag, MPI_COMM_WORLD, &status);
		}



		int rows_3D_matrix = 0; 
		if ( my_rank == 0 || my_rank == (processes - 1))
			rows_3D_matrix = rows + 1;
		else
			rows_3D_matrix = rows + 2;

		double ***intermediate_matrix;
		intermediate_matrix = (double***)malloc(sizeof(double**) * rows_3D_matrix);
		for ( i = 0; i < rows_3D_matrix; i++)
		{
			intermediate_matrix[i] = (double**)malloc(sizeof(double*) * matdim);
			for ( j = 0; j < matdim; j++)
			{
				intermediate_matrix[i][j] = (double*)malloc(sizeof(double) * matdim);
			}
		}
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
		}

		printf("intermediate matrix is calculated\n");
	
		//double result_matrix[rows_3D_matrix][matdim][matdim];
		double ***result_matrix;
                result_matrix = (double***)malloc(sizeof(double**) * rows_3D_matrix);
                for ( i = 0; i < rows_3D_matrix; i++)
                {
                        result_matrix[i] = (double**)malloc(sizeof(double*) * matdim);
                        for ( j = 0; j < matdim; j++)
                        {
                                result_matrix[i][j] = (double*)malloc(sizeof(double) * matdim);
                        }
                }
		printf("\n Process %d starts to execute stencil code\n", my_rank);
		for (i = 0; i< rows_3D_matrix; i++)
		{
			for( j = 0; j< matdim; j++)
			{
				for( k = 0; k< matdim; k++)
						result_matrix[i][j][k] = intermediate_matrix[i][j][k];
			}
		}
		printf("....");	
		for (i = 1;  i < rows_3D_matrix - 1; i++)
		{
			for ( j = 1; j < matdim - 1; j++)
			{
				for ( k = 1; k< matdim - 1 ; k++)
					result_matrix[i][j][k] = (intermediate_matrix[i][j][k-1] + intermediate_matrix[i][j][k+1] + intermediate_matrix[i][j-1][k] + intermediate_matrix[i][j+1][k] + intermediate_matrix[i-1][j][k] + intermediate_matrix[i+1][j][k] + intermediate_matrix[i][j][k])/7;

			}
		}
                printf("\n Process %d computed stencil \n", my_rank);
		for ( i = 0; i< rows_3D_matrix; i++)
                {
                        for ( j = 0; j< matdim; j++)
                        {
                        free(intermediate_matrix[i][j]);
                        }
                        free(intermediate_matrix[i]);
                }
                free(intermediate_matrix);


		int lower_bound_x = 1, upper_bound_x = rows_3D_matrix - 1;
		double *recv_updated_rows = (double*)malloc(num_elements * sizeof(double));
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
					++l;
				}
			}

		}
		printf("\n Process %d got its updated rows\n", my_rank);
		if (my_rank == 0)
		{
			result = (double*)malloc(matdim * matdim * matdim * sizeof(double));
			if ( result == NULL)
				{
					printf("\nOut of memory\n");
					MPI_Abort(MPI_COMM_WORLD, 911);
				}
			printf("\nResult buffer prepared\n");
		}
		MPI_Gather(recv_updated_rows, num_elements, MPI_DOUBLE, result, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (my_rank == 0)
		{
			printf("\nswapped\n");
			swap(&matrix, &result);
			dealloc(result);
		}

		dealloc(extra_row_left_rank);
		dealloc(extra_row_right_rank);
		dealloc(recv_updated_rows);
		for ( i = 0; i< rows_3D_matrix; i++)
		{
			for ( j = 0; j< matdim; j++)
			{
			free(result_matrix[i][j]);
			}
			free(result_matrix[i]);
		}
		free(result_matrix);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	end = MPI_Wtime();
	if(my_rank == 0)
	{
		dealloc(matrix);
		printf("With matrix size= %d and %d iterations, time elapsed: %f seconds and FLOPS: .\n" , matdim, iteration, (end - start));
	}
	dealloc(recv_rows);
	//printf("With matrix size= %d and %d iterations, time elapsed: %f seconds and FLOPS: .\n" , matdim, iteration, (end - start));
	MPI_Finalize();


	return 0;
}
