#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mpi.h"

// Size of the matrix (NxN)
#define N 64

MPI_Status status;

// Whether to print the matrix when completed
bool printResults = true;

// Print matrix function declaration
void printMatrix(int matrix[N][N]);

// Define matrices
int matrix1[N][N];
int matrix2[N][N];
int productMatrix[N][N];

// Counter variables
int i, j, k;

int main(int argc, char **argv)
{
    int numberOfProcessors;
    int processorRank;
    int numberOfWorkers;
  
    int sourceProcessor;

    int destinationProcessor;

    int rows;

    int matrixSubset;
  
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcessors);

  
    MPI_Comm_rank(MPI_COMM_WORLD, &processorRank);

    numberOfWorkers = numberOfProcessors - 1;

   

    if (processorRank == 0)
    {
     
        clock_t begin = clock();

        printf("\nMultiplying a %dx%d matrix using %d processor(s).\n\n", N, N, numberOfProcessors);

        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                matrix1[i][j] = (rand() % 6) + 1;
                matrix2[i][j] = (rand() % 6) + 1;
            }
        }

        rows = N / numberOfWorkers;
        matrixSubset = 0;

        for (destinationProcessor = 1; destinationProcessor <= numberOfWorkers; destinationProcessor++)
        {
            
            MPI_Send(&matrixSubset, 1, MPI_INT, destinationProcessor, 1, MPI_COMM_WORLD);

           
            MPI_Send(&rows, 1, MPI_INT, destinationProcessor, 1, MPI_COMM_WORLD);

           
            MPI_Send(&matrix1[matrixSubset][0], rows * N, MPI_INT, destinationProcessor, 1, MPI_COMM_WORLD);

            
            MPI_Send(&matrix2, N * N, MPI_INT, destinationProcessor, 1, MPI_COMM_WORLD);

            
            matrixSubset = matrixSubset + rows;
        }

       
        for (i = 1; i <= numberOfWorkers; i++)
        {
            sourceProcessor = i;
            MPI_Recv(&matrixSubset, 1, MPI_INT, sourceProcessor, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, sourceProcessor, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&productMatrix[matrixSubset][0], rows * N, MPI_INT, sourceProcessor, 2, MPI_COMM_WORLD, &status);
        }

       
        clock_t end = clock();

      
        if (printResults == true)
        {
            printf("Matrix 1:\n");
            printMatrix(matrix1);
            printf("Matrix 2:\n");
            printMatrix(matrix2);
            printf("Product Matrix:\n");
            printMatrix(productMatrix);
        }
      
        double runTime = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("Runtime: %f seconds\n", runTime);
    }


    if (processorRank > 0)
    {
        sourceProcessor = 0;
        MPI_Recv(&matrixSubset, 1, MPI_INT, sourceProcessor, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, sourceProcessor, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&matrix1, rows * N, MPI_INT, sourceProcessor, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&matrix2, N * N, MPI_INT, sourceProcessor, 1, MPI_COMM_WORLD, &status);

        /* Perform matrix multiplication */
        for (k = 0; k < N; k++)
        {
            for (i = 0; i < rows; i++)
            {
                productMatrix[i][k] = 0.0;
                for (j = 0; j < N; j++)
                {
                    productMatrix[i][k] = productMatrix[i][k] + matrix1[i][j] * matrix2[j][k];
                }
            }
        }

        MPI_Send(&matrixSubset, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        MPI_Send(&productMatrix, rows * N, MPI_INT, 0, 2, MPI_COMM_WORLD);
    }

    MPI_Finalize();
}


void printMatrix(int matrix[N][N])
{
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            printf("%d\t", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
