//Generates the first N prime numbers, starting at 2.
//This uses the Sieve of Eratosthenes.
//The output Y must be pre-allocated to size N = (P-1)/2, where P is the largest prime to be output.
//The actual count of primes obtained is given in cnt.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int primes_i (size_t *Y, size_t *cnt, const size_t N);


int primes_i (size_t *Y, size_t *cnt, const size_t N)
{
    //Sieve of Eratosthenes
    int8_t *sieve;
	if (!(sieve=(int8_t*)calloc((size_t)N,1))) { fprintf(stderr,"error in primes: problem with calloc. "); perror("caloc"); return 1; }

    for (size_t n=0; n<N; n++)
	{
		if (sieve[n]==0) { for (size_t m=3*n+3; m<N; m+=2*n+3) { sieve[m] = 1; } }
	}

    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    Y[0] = 2;
    *cnt = 1;
    for (size_t n=0; n<N; n++)
	{
        Y[*cnt] = 2*n + 3;
        *cnt = *cnt + (size_t)(sieve[n]==0);
		// if (sieve[n]==0)
		// {
		// 	Y[*cnt] = 2*n + 3;
        //  *cnt = *cnt + 1;
		// }
	}

    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
