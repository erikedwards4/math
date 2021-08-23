//Generates the first N prime numbers, starting at 2.
//This uses the Sieve of Eratosthenes.
//The output Y must be pre-allocated to size N+1,
//where N = (P-1)/2, and P is the largest prime to be output.
//The actual count of primes obtained is given in cnt.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int primes_i (int *Y, int *cnt, const size_t N);
int primes_u (size_t *Y, size_t *cnt, const size_t N);


int primes_i (int *Y, int *cnt, const size_t N)
{
    //Sieve of Eratosthenes
    int8_t *sieve;
	if (!(sieve=(int8_t*)calloc((size_t)N,1u))) { fprintf(stderr,"error in primes: problem with calloc. "); perror("caloc"); return 1; }

    for (size_t n=0u; n<N/3u; ++n, ++sieve)
	{
		if (*sieve==0)
        {
            for (size_t m=2u*n+3u; m<N-n; m+=2u*n+3u) { *(sieve+m) = 1; }
        }
	}
    sieve -= N/3u;

    int c = 1;
    Y[0] = 2;
    for (int n=0; n<(int)N; ++n, ++sieve)
	{
        Y[c] = 2*n + 3;
        c += (*sieve==0);
	}
    *cnt = c;

    sieve -= N; free(sieve);

    return 0;
}


int primes_u (size_t *Y, size_t *cnt, const size_t N)
{
    //Sieve of Eratosthenes
    int8_t *sieve;
	if (!(sieve=(int8_t*)calloc((size_t)N,1u))) { fprintf(stderr,"error in primes: problem with calloc. "); perror("caloc"); return 1; }

    for (size_t n=0u; n<N/3u; ++n, ++sieve)
	{
		if (*sieve==0)
        {
            for (size_t m=2u*n+3u; m<N-n; m+=2u*n+3u) { *(sieve+m) = 1; }
        }
	}
    sieve -= N/3u;

    size_t c = 1u;
    Y[0] = 2u;
    for (size_t n=0u; n<N; ++n, ++sieve)
	{
        Y[c] = 2u*n + 3u;
        c += (size_t)(*sieve==0);
	}
    *cnt = c;


    sieve -= N; free(sieve);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
