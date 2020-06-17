//Generates the first N prime numbers, starting at 2.
//This uses the Sieve of Eratosthenes.
//The output Y must be pre-allocated to size N = (P-1)/2, where P is the largest prime to be output.
//The actual count of primes obtained is given in cnt.

#include <stdio.h>
#include <stdlib.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int primes_i (int *Y, int *cnt, const int N);


int primes_i (int *Y, int *cnt, const int N)
{
    int m, n;
	int8_t *sieve;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<0) { fprintf(stderr,"error in primes_i: N (num elements Y) must be nonnegative\n"); return 1; }

    //Sieve of Eratosthenes
	if (!(sieve=(int8_t*)calloc((size_t)N,1))) { fprintf(stderr,"Error in primes: problem with calloc. "); perror("caloc"); return 1; }
	for (n=0; n<N; n++)
	{
		if (sieve[n]==0) { for (m=3*n+3; m<N; m+=2*n+3) { sieve[m] = 1; } }
	}
    Y[0] = 2;
    *cnt = 1;
    for (n=0; n<N; n++)
	{
		if (sieve[n]==0)
		{
			Y[*cnt] = (int)(2*n + 3);
            *cnt = *cnt + 1;
		}
	}

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
