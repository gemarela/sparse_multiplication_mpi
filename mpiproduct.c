#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mmio.h"
#include "mpi.h"
#include <math.h>

//#define Numprocs 2
//#define MyRank 0
//#define ROW 6
//#define COLLUM 6


int main (int argc, char *argv[]){
  
	int ret_code;
	int M, N, nz;   //M : rows , N : collums
	int i,j,index,irow, icol, iproc;
	int Numprocs,MyRank;
	int Root = 0;
	int tag1 = 1;
	int tag2 = 2;
    	int blocks;
	int row_start;
	int row_end;
	int *i_index, *j_index;
	double *a_values;
	double *b_vector;
	double *dot_product;
	double *temp_product;
	double seconds;

	clock_t t;
    
	MM_typecode matcode;
	FILE *f,*fp;
    
	srand(time(NULL));
	
		  /* ........MPI Initialisation .......*/
		  
	MPI_Status status;
	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
	MPI_Comm_size(MPI_COMM_WORLD, &Numprocs);
		
	if (MyRank == 0){	
	
		fp = fopen("results_coo.txt","w");
			if(fp==NULL){
				printf("Error\n");
				exit(0);
			}
    
		if (argc < 2)
		{
			fprintf(stderr, "Usage: %s [matrix-market-filename]\n", argv[0]);
			MPI_Finalize();
			exit(1);
		}
		else    
		{ 
			if ((f = fopen(argv[1], "r")) == NULL) 
				exit(1);
		}

		if (mm_read_banner(f, &matcode) != 0)
		{
			printf("Could not process Matrix Market banner.\n");
			exit(1);
		}


		/*  This is how one can screen matrix types if their application */
		/*  only supports a subset of the Matrix Market data types.      */

		if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
			mm_is_sparse(matcode) )
		{
			printf("Sorry, this application does not support ");
			printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
			exit(1);
		}

		/* find out size of sparse matrix .... */

		if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
			exit(1);


		/* reseve memory for matrices */

		i_index = (int *) calloc(nz,sizeof(int));
		j_index = (int *) calloc(nz,sizeof(int));
		a_values = (double *) calloc(nz,sizeof(double));


		/* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
		/*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
		/*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

		for (i=0; i<nz; i++){
		
			fscanf(f, "%d %d %lg\n", &i_index[i], &j_index[i], &a_values[i]);
			i_index[i]--;  /* adjust from 1-based to 0-based */
			j_index[i]--;
		}

		if (f !=stdin) fclose(f);

	
		/************************/
		/* now write out matrix */
		/************************/

		mm_write_banner(stdout, matcode);
		mm_write_mtx_crd_size(stdout, M, N, nz);
		for (i=0; i<nz; i++)
			fprintf(stdout, "%d %d %20.19g\n", i_index[i]+1, j_index[i]+1, a_values[i]);
    
		/******************************************************************/
		/************MEM ALLOCATE FOR B_VECTOR AND INITIALIZE**********/
	
		fprintf(fp,"\nnz = %d\n",nz);
	
		b_vector = (double*)calloc(M,sizeof(double));
		for(i=0;i<M;i++){
			b_vector[i] = (rand()%9)+1;
		}
		//fprintf(fp,"\nB_VECTOR : \n");
		//for(i=0;i<nz;i++){
		//	fprintf(fp,"b_vector[%d] = %.2lf \n",i,b_vector[i]);
		//}
		//fprintf(fp,"\n");
	
		//fprintf(fp,"Compressed A with COO method : \n");
		//for(i=0;i<nz;i++){
		//	fprintf(fp,"a_values[%d] = %.2lf\t",i,a_values[i]);
		//	fprintf(fp,"i_index[%d] = %d\t",i,i_index[i]);
		//	fprintf(fp,"j_index[%d] = %d\n",i,j_index[i]);
		//}
		//fprintf(fp,"\n");	
		//fclose(fp);
		
		/*******************************************************************/
		
		
		if(N < Numprocs) {
			if(MyRank == 0)
				printf("No of Rows should be more than No of Processors ... \n");
			MPI_Finalize();
			exit(0);
		}  	

		if(N % Numprocs != 0) {
			if(MyRank == 0) 
				printf("Matrix Can not be Striped Evenly ..... \n");
			MPI_Finalize();
			exit(0);
		} 	
	
	}

	MPI_Barrier(MPI_COMM_WORLD); // alla processes reach a point before execution

	t = clock();
	
	MPI_Bcast(&nz, 1, MPI_INT, Root, MPI_COMM_WORLD);
	MPI_Bcast(&M, 1, MPI_INT, Root, MPI_COMM_WORLD);
	N = M;
	blocks = M / Numprocs;
	if (MyRank != 0) {
		i_index = (int *) calloc(nz,sizeof(int));
		j_index = (int *) calloc(nz,sizeof(int));
		a_values = (double *) calloc(nz,sizeof(double));
		b_vector = (double*)calloc(M,sizeof(double));
	}

	MPI_Bcast(i_index, nz, MPI_INT, Root, MPI_COMM_WORLD);
	MPI_Bcast(j_index, nz, MPI_INT, Root, MPI_COMM_WORLD);
	MPI_Bcast(a_values, nz, MPI_DOUBLE, Root, MPI_COMM_WORLD);
	MPI_Bcast(b_vector, M, MPI_DOUBLE, Root, MPI_COMM_WORLD);
	temp_product = (double*)calloc(M,sizeof(double));
	dot_product = (double*)calloc(M,sizeof(double));

	index=0;
	//i=MyRank+1;

	row_start = MyRank * blocks;
	row_end = (MyRank + 1) * blocks;
	printf("Row start: %d\n", row_start);
	printf("Row end: %d\n", row_end);
	for(i = 0; i < nz; i++){
		//fprintf(fp,"i_index[%d]: %d\n", i, i_index[i]);
	//Ypologizontai ta epimerous eswterika ginomena
		if (i_index[i] < row_start)
			continue;
		if (i_index[i] >= row_end)
			continue;
		//fprintf(fp,"a_values[%d]: %f\n", i, a_values[i]);
		//fprintf(fp,"b_vector[%d]: %f\n", j_index[i], b_vector[j_index[i]]);
		temp_product[i_index[i]] += a_values[i]*b_vector[j_index[i]];
	}

	//for(i=0;i<N;i++){
	//	printf("\ndot_product[%d] = %.2lf\n",i,temp_product[i]);
	//}

	// syllogh stoixeiwn
	if (MyRank != 0) {
		MPI_Send(temp_product, M, MPI_DOUBLE, Root, tag2, MPI_COMM_WORLD);
	}
	else {
		memcpy(dot_product, temp_product, M * sizeof(double));
		for (i = 1; i < Numprocs; i++) { 
			MPI_Recv(temp_product, M, MPI_DOUBLE, i, tag2, MPI_COMM_WORLD, &status);
			for (j = 0; j < M; j++) {
				dot_product[j] += temp_product[j];
			}
		} 
	}

	if (MyRank == 0) {
		for(i=0;i<N;i++){
			fprintf(fp,"\ndot_product[%d] = %.2lf\n",i,dot_product[i]);
		}
	}
	
	t = clock() - t;
	seconds = ((double)t)/CLOCKS_PER_SEC;
	printf("\nElapsed time is %lf seconds with %d Processors\n",seconds,Numprocs);
	
	if (MyRank == 0)
		fclose(fp);
	MPI_Finalize();
	
		
	free(a_values);
	free(i_index);
	free(j_index);
	free(b_vector);
	free(dot_product);
	
	return 0;
}
/*int proc_map(int i, int Numprocs , int M)
{
    Numprocs = Numprocs - 1;
    int r = (int) ceil( (double)M / (double)Numprocs);
    int proc = i / r;
    return proc + 1;
}*/


