/* Program to perform matrix-matrix multiplication */

#include <stdio.h>
#include <stdlib.h>

void allocate_matrix(double ***matrix, int nrows, int ncolumns);
void initialize_matrix(double **matrix, int nrows, int ncolumns);
void display_matrix(double **matrix, int nrows, int ncolumns);
void multiply_matrices(double **matrix1, int nrows_matrix1, double **matrix2, int ncolumns_matrix2, double **matrix3, int nrows_matrix3);
void free_memory(double **matrix, int nrows);

int main ()
{
	double **matrix1, **matrix2, **matrix3;
	int nrows_matrix1, nrows_matrix2, nrows_matrix3, ncolumns_matrix1, ncolumns_matrix2, ncolumns_matrix3;
	printf("\nEnter the no. of rows of first matrix: ");
        scanf("%d", &nrows_matrix1);
	printf("\nEnter the no. of columns of first matrix: ");
	scanf("%d", &ncolumns_matrix1);
        printf("\nEnter the no. of rows of second matrix: ");
        scanf("%d", &nrows_matrix2);
	printf("\nEnter the no. of columns of second matrix: ");
	scanf("%d", &ncolumns_matrix2);
	printf("\nrows and columns: %d %d %d %d", nrows_matrix1, ncolumns_matrix1, nrows_matrix2, ncolumns_matrix2); 
        if (ncolumns_matrix1 != nrows_matrix2)
        {
                printf("\nMatrices can't be multiplied!!\n");
                exit;
        }
	else
	{
		nrows_matrix3 = nrows_matrix1;
		ncolumns_matrix3 = ncolumns_matrix2;
		allocate_matrix(&matrix1, nrows_matrix1, ncolumns_matrix1);
		allocate_matrix(&matrix2, nrows_matrix2, ncolumns_matrix2);	
		allocate_matrix(&matrix3, nrows_matrix3, ncolumns_matrix3);
		initialize_matrix(matrix1, nrows_matrix1, ncolumns_matrix1);
		initialize_matrix(matrix2, nrows_matrix2, ncolumns_matrix2);
		printf("\nFirst matrix: \n");
		display_matrix(matrix1, nrows_matrix1, ncolumns_matrix1);
		printf("\nSecond matrix: \n");
		display_matrix(matrix2, nrows_matrix2, ncolumns_matrix2);
		multiply_matrices(matrix1, nrows_matrix1, matrix2, ncolumns_matrix2, matrix3, ncolumns_matrix1);
		printf("\nResult matrix: \n");
		display_matrix(matrix3, nrows_matrix3, ncolumns_matrix3);
		printf("\nFreeing memory of matrix 1");
		free_memory(matrix1, nrows_matrix1);
		printf("\nFreeing memory of matrix 2");
		free_memory(matrix2, nrows_matrix2);
		printf("\nFreeing memory of matrix 3");
		free_memory(matrix3, nrows_matrix3);
		printf("\nEnd of program");
	}
	
	return 0;
}

void allocate_matrix(double ***matrix, int nrows, int ncolumns)
{
	printf("\nStarting to allocate memory");
	int i;
	*matrix =(double **)malloc(nrows * sizeof(double *));
	for (i = 0; i< nrows; i++)
	{
		 (*matrix)[i] = (double*)malloc(ncolumns * sizeof(double));
	}
}



void initialize_matrix(double **matrix, int nrows, int ncolumns)
{
	
	int i, j;
	for (i = 0; i< nrows; i++)
	{
		for (j = 0; j < ncolumns; j++)
			matrix[i][j] = drand48();
	}
}

void multiply_matrices(double **matrix1, int nrows_matrix1, double **matrix2, int ncolumns_matrix2, double **matrix3, int ncolumns_matrix1)
{
	int i, j, k;
	for (i = 0; i< nrows_matrix1; i++)
	{
		for (j = 0; j< ncolumns_matrix2; j++)
		{
			matrix3[i][j] = 0;
			for (k = 0; k < ncolumns_matrix1; k++)
			{
				matrix3[i][j] +=  matrix1[i][k] * matrix2[k][j];
			}

		}
	}
}

void display_matrix(double **matrix, int nrows, int ncolumns)
{
	int i, j;
	for (i = 0; i < nrows; i++)
	{
		for(j = 0; j < ncolumns; j++)
		{
			printf("%f ", matrix[i][j]);
		}
		printf("\n");
	}
}

void free_memory(double **matrix, int nrows)
{
	printf("\n Freeing memory \n");
	int i;
	for (i = 0; i< nrows; i++)
		free(matrix[i]);
	free(matrix);
}

