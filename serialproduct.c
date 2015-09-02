#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mmio.h"

//#define ROW 6
//#define COLLUM 6


void coo_func(int *i_index, int *j_index, double *a_values, double *b_vector, int M, int nz , FILE *fp){
	
	int i,j;
	double *dot_product = (double*)calloc(M,sizeof(double));

	
	fprintf(fp,"Compressed A with COO method : \n");
	for(i=0;i<nz;i++){
		fprintf(fp,"a_values[%d] = %.2lf\t",i,a_values[i]);
		fprintf(fp,"i_index[%d] = %d\t",i,i_index[i]);
		fprintf(fp,"j_index[%d] = %d\n",i,j_index[i]);
	}
	fprintf(fp,"\n");
		
	for(i=0;i<nz; i++){

			dot_product[i_index[i]] += a_values[i] * b_vector[j_index[i]]; 
		
	}
		
	for(i=0;i<M;i++){
		fprintf(fp,"\ndot_product[%d] = %.2lf",i,dot_product[i]);
	}
		
	free(a_values);
	free(i_index);
	free(j_index);
  
}

int main (int argc, char *argv[]){
  
    int ret_code;
    int M, N, nz;   //M : rows , N : collums
    int i,j, *i_index, *j_index;
    double *a_values;
    double *b_vector;
    
    MM_typecode matcode;
    FILE *f,*fp;
    
    srand(time(NULL));
	
    fp = fopen("results serial coo.txt","w");
	if(fp==NULL){
		printf("Error\n");
		exit(0);
	}
    
    if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [matrix-market-filename]\n", argv[0]);
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

    for (i=0; i<nz; i++)
    {
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
    
    /***************************************************/

	
	fprintf(fp,"\nnz = %d\n",nz);
	
	b_vector = (double*)calloc(nz,sizeof(double));
	for(i=0;i<nz;i++){
		b_vector[i] = (rand()%9)+1;
	}
	fprintf(fp,"\nB_VECTOR : \n");
	for(i=0;i<N;i++){
		fprintf(fp,"b_vector[%d] = %.2lf \n",i,b_vector[i]);
	}
	fprintf(fp,"\n");
	
	coo_func(i_index,j_index,a_values,b_vector,M,nz,fp);
    
	
	return;
}
