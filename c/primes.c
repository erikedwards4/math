//Generates the first N prime numbers, starting at 2.
//This uses the Sieve of Eratosthenes.
//The output Y must be pre-allocated to size N+1,
//where N = (P-1)/2, and P is the largest prime to be output.
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

    for (size_t n=0u; n<N/3; ++n, ++sieve)
	{
		if (*sieve==0)
        {
            for (size_t m=2*n+3; m<N-n; m+=2*n+3) { *(sieve+m) = 1; }
        }
	}
    sieve -= N/3;

    Y[0] = 2; *cnt = 1;
    for (size_t n=0u; n<N; ++n, ++sieve)
	{
        Y[*cnt] = 2*n + 3;
        *cnt = *cnt + (size_t)(*sieve==0);
	}

    sieve -= N; free(sieve);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
