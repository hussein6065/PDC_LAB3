#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

// Size of the square Matrix N
#define n 10


// Matrice
double matrixA[n][n];
double matrixB[n][n]; 
double matrixC[n][n];




void instantiateRandomMatrix();
void print(int l);


int main(int argc, char* argv[])
{

//    Time to keep track of the execution time;
    clock_t first_time,last_time;
    float diff;
    
    int my_rank, p;  // process rank and number of processes
    int source = 0, dest,dataCount,offset,counts;  
    int i,j,k,m,c,u,v;
    int tag = 1;  // tag for messages

   
   
    
    MPI_Status status; // stores status for MPI_Recv statements
    MPI_Datatype columntype,rowtype,resultType; // The new MPI datatypes

    
    instantiateRandomMatrix();
   
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

  
    // Partitioning and sending the sub matrix.
    offset = n/sqrt(p-1);
    dataCount = offset*n;
    // creating new MPI_Datatype

    MPI_Type_vector(dataCount ,n, n, MPI_DOUBLE, &rowtype);
    MPI_Type_commit(&rowtype);
    MPI_Type_vector(dataCount ,offset, n, MPI_DOUBLE, &columntype);
    MPI_Type_commit(&columntype);


   
    

    if(my_rank == 0)
    {

        
        first_time = clock();
        dest = 1;
        

        k = sqrt(p-1);
        // printf("The sqrt: %d\n",k);

        // Sending the data partitioned data to the process.
        for (i = 0; i < k; i++){
            for (j = 0; j < k; j++){
                
                MPI_Send( &matrixA[i*offset][0] ,1, rowtype , dest , 1 , MPI_COMM_WORLD);
                MPI_Send( &matrixB[0][j*offset] ,1, columntype , dest , 2 , MPI_COMM_WORLD);
                
                dest++;
            }
        }

       

    //    Probing to get the count and the size of the incoming data.
        MPI_Probe( 1 , 0, MPI_COMM_WORLD , &status);
        MPI_Get_count( &status , MPI_DOUBLE , &counts);
        c = status.count_lo;

       
        double result[counts];

        int start;

        // Taking the incoming data and putting it in to the MatrixC
        for (source = 1; source < p; source++){
            

            start = (source-1)*offset;
            i = ((start)/n)*offset;
            j = (start)%n;
            
            // matrixC[i][j] = source;
            MPI_Recv(&result, c , MPI_DOUBLE, source , 0 , MPI_COMM_WORLD, &status);
            for (int k = 0;k<offset*offset;k++){
             
                u = (k/offset)+i;
                v = (k%offset)+j;
                matrixC[u][v] = result[k];
            //    printf("%f\n",result[k]);
            }
            // printf("\n");
            
                
                // printf("Index: %d value: %d source: %d\n",i,result[i],source);
            //
        }
        // Printing the time it tool to execute the Multiplication.
        last_time = clock();
        diff = (float)(last_time-first_time)/CLOCKS_PER_SEC;
        printf("This is the time used: %fs\n",diff);
        // printf("\n");
        // print(2);


      

     
        

     
    
     
    }else 
    {
        // int column[d];
        // int row[d];
        // int *column;
        // int *row;
        
        
     
        // Probing to get the count and size of the incoming data.
        MPI_Probe( source , tag, MPI_COMM_WORLD , &status);
        MPI_Get_count( &status , MPI_DOUBLE , &counts);
        
        c = status.count_lo;
        
        double row[c], column[c],send[offset*offset];
        /
        
        double result[offset][offset];

        // Creating another MPI dataType.
        MPI_Type_vector(offset*offset,offset, offset, MPI_DOUBLE, &resultType);
        MPI_Type_commit(&resultType);
        

        double subColumn[n][offset],subRow[offset][n];
        
        // Receving the data from process 0
        MPI_Recv(&row , c, MPI_DOUBLE , source , 1 , MPI_COMM_WORLD, &status);
        MPI_Recv(&column ,c, MPI_DOUBLE , source , 2 , MPI_COMM_WORLD, &status);
        


        // The convertion to the row
        for (i = 0; i <  dataCount;i++ ){
            k= i/n;
            m = i%n;
            subRow[k][m] = row[i];
            // printf("i: %d, j: %d index: %d\n",m,k,i);
        }

        // The convertion to the columns
        for (i = 0; i <  offset;i++ ){
            k = 0;
            for (j = i; j < dataCount; j+=offset){
                // printf("i: %d, j: %d value: %d\n",i,j,column[j]);
                subColumn[k++][i] = column[j];
            }        
        }

        
        //   Multiplication of the ,matrix
        for (i = 0; i < offset; i++)
            {
                for (j = 0; j < offset; j++)
                    {
                   
                    result[i][j] = 0;
                    for( k = 0; k < n; k++)
                        {
                            
                            result[i][j] += subRow[i][k] * subColumn[k][j];
                            
                        }
                       
                }
               
            }
    //    sending the results to process 0;
        MPI_Send( &result[0][0] ,1, resultType, 0, 0 , MPI_COMM_WORLD);

        // Freeing the datatype;
        MPI_Type_free( &resultType);
       

    }
    
    // Freeing the datatypes;
    MPI_Type_free( &rowtype);
    MPI_Type_free( &columntype);
    
    MPI_Finalize();
    
   


}

// creating the random matrix
void instantiateRandomMatrix()
{
   
    
    int i,j;

    // for( i = 0; i < n; i++)
    // {
    //     matrixA[i] = (int *) malloc(n*sizeof(int));
    //     matrixB[i] = (int *) malloc(n*sizeof(int));
    //     matrixC[i] = (int *) malloc(n*sizeof(int));
    // }
        

    // srand(time(0));
    // subRow[offset][n];
    int num = 0;
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            // matrixA[i][j] = rand()%10000;
            // matrixB[i][j] = rand()%10000;
            
            matrixA[i][j] = ++num;
            matrixB[i][j] = num;
            // matrixC[i][j] = 0;
        }
    }
   
    
}

// function used to print the matrix
void print(int l)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            
            if(l == 1)
                printf("%f ",matrixA[i][j]);
            else if(l==2)
                printf("%f ",matrixB[i][j]);
            else
                printf("%f ",matrixC[i][j]);
        }
        printf("\n");
    }
}
