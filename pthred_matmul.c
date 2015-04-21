#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<pthread.h>
#include<string.h>

int **result_sequential, **result_parallel;
int **mat1, **mat2;
long long matdim = 0;
int threadcount=0;
void **retval;
pthread_t *threads;

void *matmul(void *arg){

	long long i, j, k, beg, end;
	int **result = result_sequential;
	int tid = (intptr_t)arg;

	beg = 0;
	end = matdim;

	if(tid > -1){
		beg = ((tid * matdim)/threadcount);
		end = (((tid + 1) * matdim)/threadcount);
		result = result_parallel;
	}

	for (i = beg; i < end; i++){
		for (j = 0; j < matdim; j++){
			for (k = 0; k < matdim; k++)
				result[i][j] +=  mat1[i][k] * mat2[k][j];
		}
	}
	
	return NULL;
}

void init(){

	long long i,j;

	mat1 = (int **) calloc(matdim, sizeof(int*));
	mat2 = (int **) calloc(matdim, sizeof(int*));
	result_sequential = (int **) calloc(matdim, sizeof(int*));
	result_parallel = (int **) calloc(matdim, sizeof(int*));

	for(i = 0; i < matdim; ++i){
                mat1[i] = (int *) calloc(matdim, sizeof(int));
                mat2[i] = (int *) calloc(matdim, sizeof(int));
                result_sequential[i] = (int *) calloc(matdim, sizeof(int));
                result_parallel[i] = (int *) calloc(matdim, sizeof(int));
	}

	srand(1000);

        for(i = 0; i < matdim; ++i)
                for(j = 0; j < matdim; ++j){
                        mat1[i][j] = rand();
			mat2[i][j] = rand();
		}
}

void dealloc(){

	long long i;

	for(i = 0; i < matdim; ++i){
                free(mat1[i]);
		free(mat2[i]);
		free(result_sequential[i]);
		free(result_parallel[i]);
	}

	free(mat1);
	free(mat2);
	free(result_sequential);
	free(result_parallel);
}

int main(int argc, char *argv[]){

	int tid;
        long long i, j;
        struct timespec start, finish;
        double elapsed;

	//default thread count and matrix dimensions
	matdim = 1000;
	threadcount = 1;

	//setting thread count and matrix dimensions to that passed by user through command line
	if((argc > 1) && (argv[1] != NULL) &&(strlen(argv[1]) > 0)){
		matdim = atoll(argv[1]);
	}

	if((argc > 2) && (argv[2] != NULL) &&(strlen(argv[2]) > 0)){
		threadcount = atoi(argv[2]);
	}

	//initialization
	init();

	printf("Initialization Complete \n");

	//sequential
        clock_gettime(CLOCK_MONOTONIC, &start);

	matmul((void*)(intptr_t)-1);

        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("Sequential Elapsed Time: %f\n", elapsed);

	//parallel
	threads = (pthread_t*) calloc(threadcount, sizeof(pthread_t));

	clock_gettime(CLOCK_MONOTONIC, &start);

	for(tid = 0; tid < threadcount; ++tid){
		pthread_create(&threads[tid], NULL, matmul, (void*)(intptr_t)tid);
	}

	for(tid = 0; tid < threadcount; ++tid)
		pthread_join(threads[tid], retval);

        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("Parallel Elapsed Time: %f\n", elapsed);

	//result verification
	for(i=0; i< matdim; ++i)
		for(j=0;j< matdim;++j){
			if(result_sequential[i][j] != result_parallel[i][j]){
				printf("Invalid Result\n");
				dealloc();
				exit(0);
			}
		}


	//deallocating memory
	dealloc();

	printf("Deallocation Complete \n");
        return 0;
}
