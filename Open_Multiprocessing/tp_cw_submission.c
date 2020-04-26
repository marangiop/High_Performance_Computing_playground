//
//The objective of this piece of work is to implement the affinity scheduling algorithm (https://ieeexplore.ieee.org/document/344281) using the Open Multiprocessing (OpenMP) framework.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 729
#define reps 1000
#include <omp.h>

double a[N][N], b[N][N], c[N];
int jmax[N];


void init1(void);
void init2(void);
void runloop(int);
void loop1chunk(int, int);
void loop2chunk(int, int);
void valid1(void);
void valid2(void);
int largest(int*, int);


int main(int argc, char *argv[]) {

    double start1,start2,end1,end2;
    int r;

    init1();

    start1 = omp_get_wtime();

    for (r=0; r<reps; r++){
        runloop(1);
    }

    end1  = omp_get_wtime();

    valid1();

    printf("Total time for %d reps of loop 1 = %f\n",reps, (float)(end1-start1));


    init2();

    start2 = omp_get_wtime();

    for (r=0; r<reps; r++){
        runloop(2);
    }

    end2  = omp_get_wtime();

    valid2();

    printf("Total time for %d reps of loop 2 = %f\n",reps, (float)(end2-start2));

}

void init1(void){
    int i,j;

    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            a[i][j] = 0.0;
            b[i][j] = 3.142*(i+j);
        }
    }

}

void init2(void){
    int i,j, expr;

    for (i=0; i<N; i++){
        expr =  i%( 3*(i/30) + 1);
        if ( expr == 0) {
            jmax[i] = N;
        }
        else {
            jmax[i] = 1;
        }
        c[i] = 0.0;
    }

    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            b[i][j] = (double) (i*j+1) / (double) (N*N);
        }
    }

}


void runloop(int loopid) {
    int id_of_loaded_thread;
    int n;
    int myid;
    int i;
    int ipt;
    int lo;
    int hi;
    int inner_lo;
    int inner_hi;
    int nthreads = omp_get_max_threads();
    int remaining_iterations[16];
    int last_iteration[16];
    int end_of_programme = 0;
    int var_1;
    int var_2;


#pragma omp parallel default(none) shared(remaining_iterations, last_iteration, nthreads, ipt, end_of_programme) private(myid, inner_lo, inner_hi, lo, hi, id_of_loaded_thread, loopid, var_1, var_2)
    {
        int myid = omp_get_thread_num();

        int ipt = (int)ceil((double) N / (double) nthreads);
        int lo = myid * ipt;
        int hi = (myid + 1) * ipt;
        if (hi > N) hi = N;

        remaining_iterations[myid] = hi-lo;

        inner_lo = lo;
        inner_hi = lo + (int) ceil((double) remaining_iterations[myid]/(double) nthreads);

        var_1=inner_lo;
        var_2=inner_hi;


            #pragma omp barrier
            while(end_of_programme==0) {

            switch (loopid) {
                case 1:
                    loop1chunk(var_1, var_2);
                    break;
                case 2:
                    loop2chunk(var_1, var_2);
                    break;
            }


            #pragma omp barrier
            #pragma omp critical
                {
                    if (remaining_iterations[myid] != 0) { //If the thread has to still carry out some iterations from its local set
                        remaining_iterations[myid] = remaining_iterations[myid] - (inner_hi - inner_lo);
                        last_iteration[myid] = inner_hi;
                        inner_lo = inner_hi;
                        inner_hi = inner_lo + (int) ceil((double) remaining_iterations[myid]/(double) nthreads);
                    } else { //If the thread has no more iterations in its local set
                        id_of_loaded_thread = largest(remaining_iterations, nthreads); //Obtain the id of the most heavily loaded thread
                        if (id_of_loaded_thread == -1) { //If all other threads are also finished with their local set iterations
                            end_of_programme = 1; //Flag to end the programme when it leaves the critical region
                        } else { //Steal a chunk of iterations from the most loaded thread
                            inner_lo = last_iteration[id_of_loaded_thread];
                            inner_hi = inner_lo + (int)ceil((double)remaining_iterations[id_of_loaded_thread] /(double) nthreads);
                            remaining_iterations[id_of_loaded_thread] = remaining_iterations[id_of_loaded_thread] - (inner_hi - inner_lo);
                        }
                    }
                        var_1 = inner_lo;
                        var_2 = inner_hi;
                }
        }
        }
    }



void loop1chunk(int lo, int hi) {
    int i,j;

    for (i=lo; i<hi; i++){
        for (j=N-1; j>i; j--){
            a[i][j] += cos(b[i][j]);
        }
    }

}



void loop2chunk(int lo, int hi) {
    int i,j,k;
    double rN2;

    rN2 = 1.0 / (double) (N*N);

    for (i=lo; i<hi; i++){
        for (j=0; j < jmax[i]; j++){
            for (k=0; k<j; k++){
                c[i] += (k+1) * log (b[i][j]) * rN2;
            }
        }
    }

}

void valid1(void) {
    int i,j;
    double suma;

    suma= 0.0;
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            suma += a[i][j];
        }
    }
    printf("Loop 1 check: Sum of a is %lf\n", suma);

}


void valid2(void) {
    int i;
    double sumc;

    sumc= 0.0;
    for (i=0; i<N; i++){
        sumc += c[i];
    }
    printf("Loop 2 check: Sum of c is %f\n", sumc);
}

int largest(int* arr, int n)
{
    int i;

    // Initialize maximum element
    int max = 0;
    int max_i = 0;

    // Traverse array elements from second and
    // compare every element with current max

    for (i = 0; i < n; i++) {
        if (arr[i] > max) {
            max = arr[i];
            max_i = i;
        }
    }

    if(0 == max) {
        max_i = -1;
    }
    return max_i;
}


