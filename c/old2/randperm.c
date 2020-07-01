//Generates integers from 1 to N in random order, and outputs the first M.
//This uses method of M Knuth shuffles.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int randperm_s (float *Y, const int M, const int N);
int randperm_d (double *Y, const int M, const int N);


int randperm_s (float *Y, const int M, const int N)
{
    int m, n;
    int *X;
    struct timespec ts;

    //Checks
    if (N<1) { fprintf(stderr,"error in randperm_s: N (num elements X) must be positive\n"); return 1; }

    //Generate ints 1:N
    if (!(X=(int *)malloc((size_t)N*sizeof(int)))) { fprintf(stderr,"error in randperm_s: problem with malloc. "); perror("malloc"); return 1; }
    for (n=0; n<N; n++) { X[n] = n + 1; }

    //Seed rand
	if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in randperm_s: timespec_get.\n"); perror("timespec_get"); return 1; }
	else { srand48(ts.tv_nsec^ts.tv_sec); }

    //M Knuth shuffles
    for (m=0; m<M; m++)
	{
		n = m + (int)((N-m)*drand48());		//random index from m to N-1
        Y[m] = (float)X[n]; X[n] = X[m];    //combine swap and output
	}

    free(X);
    return 0;
}


int randperm_d (double *Y, const int M, const int N)
{
    int m, n;
    int *X;
    struct timespec ts;

    //Checks
    if (N<1) { fprintf(stderr,"error in randperm_d: N (num elements X) must be positive\n"); return 1; }

    //Generate ints 1:N
    if (!(X=(int *)malloc((size_t)N*sizeof(int)))) { fprintf(stderr,"error in randperm_d: problem with malloc. "); perror("malloc"); return 1; }
    for (n=0; n<N; n++) { X[n] = n + 1; }

    //Seed rand
	if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in randperm_d: timespec_get.\n"); perror("timespec_get"); return 1; }
	else { srand48(ts.tv_nsec^ts.tv_sec); }

    //M Knuth shuffles
    for (m=0; m<M; m++)
	{
		for (m=0; m<M; m++)
	{
		n = m + (int)((N-m)*drand48());		//random index from m to N-1
        Y[m] = (double)X[n]; X[n] = X[m];   //combine swap and output
	}
	}

    free(X);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
