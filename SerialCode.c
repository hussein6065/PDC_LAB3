#include <stdlib.h>
#include <stdio.h>
#include <time.h>



int size =10;

int main(int argc, char* argv[])
{
    // Matricx
    double MatrixA[size][size];
    double MatrixB[size][size];
    double MatrixC[size][size];

    int i,j,k,s =0;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            // MatrixA[i][j] = rand()%10000;
            // MatrixB[i][j] = rand()%10000;
            MatrixA[i][j] = ++s;
            MatrixB[i][j] = s;

        }
    }

    // Matrix multiplicaiton
    clock_t first_time,last_time;
    float diff;
    first_time = clock();
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
           MatrixC[i][j] = 0;
           for( k = 0; k < size; k++)
           {
               MatrixC[i][j] +=  MatrixA[i][k] * MatrixB[k][j];
               
           }
            // printf("%d \n",MatrixC[i][j]);

        }
    }
    // printf("\n I am done .....");
    // printf("\n");
    
    // int N, myrank;
    // MPI_Init(&argc,&argv);
    // MPI_Comm_size(MPI_COMM_WORLD,&N);
    // MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    // printf("From process %d out of %d, This is a new World^\n",myrank,N);
    // if(myrank !=0){
        
    // }
    // MPI_Finalize();
    // for(int i = 0; i < size; i++)
    // {
    //     for(int j = 0; j < size; j++)
    //     {
            
            
    //         printf("%d ",MatrixA[i][j]);
    //     }
    //     printf("\n");
    // }
    // for(int i = 0; i < size; i++)
    // {
    //     for(int j = 0; j < size; j++)
    //     {
            
            
    //         printf("%d ",MatrixB[i][j]);
    //     }
    //     printf("\n");
    // }
    last_time = clock();
    diff = (float)(last_time-first_time)/CLOCKS_PER_SEC;

    // Printing the execution type.
    printf("This is the time used: %fs\n",diff);
    printf("I am done ....\n");
    // for(int i = 0; i < size; i++)
    // {
    //     for(int j = 0; j < size; j++)
    //     {
            
            
    //         printf("%d ",MatrixC[i][j]);
    //     }
    //     printf("\n");
    // }
    
    
    
    
    
    
  
        
    return 0;
}

