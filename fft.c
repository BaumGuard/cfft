#include "fft.h"

#include <stdlib.h>
#include <math.h>

#define TWO_PI 2.0 * M_PI

/*---------------------------------------------------------------------------*/

struct fft_t {
    double complex* tmp;
    double complex* ifft_values;
    int n_samples;
};

/*---------------------------------------------------------------------------*/

void FFT_Init( fft_t** handle, int n_samples ) {
    (*handle) = (fft_t*) malloc( sizeof(fft_t) );
    (*handle)->tmp = (double complex*) calloc( n_samples, sizeof(double complex) );
    (*handle)->ifft_values = (double complex*) calloc( n_samples, sizeof(double complex) );
    (*handle)->n_samples = n_samples;
} /* FFT_Init() */

/*---------------------------------------------------------------------------*/

void FFT_iterative_inner( fft_t* handle, double complex* result ) {
    int recursion_depth = log2( handle->n_samples );
    int N = handle->n_samples;
    int M = N / 2;
    for ( int i = 0; i < recursion_depth; i++ ) {
        for ( int j = 0; j < handle->n_samples; j += N ) {
            int even_idx = 0, odd_idx = M;
            for ( int k = j; k < j+N; k++ ) {
                if ( k % 2 == 0 ) {
                    handle->tmp[even_idx++] = result[k];
                }
                else {
                    handle->tmp[odd_idx++]  = result[k];
                }
            }

            for ( int k = 0; k < N; k++ ) {
                result[k+j] = handle->tmp[k];
            }
        }

        N /= 2;
        M /= 2;
    }

    N = 2;
    M = 1;

    for ( int i = 0; i < recursion_depth; i++ ) {
        for ( int j = 0; j < handle->n_samples; j += N ) {
            for ( int k = 0; k < M; k++ ) {
                float angle = TWO_PI * ((float)k / (float)N);

                handle->tmp[k]   = result[j+k] + result[j+k+M] * cexp( angle * I );
                handle->tmp[k+M] = result[j+k] - result[j+k+M] * cexp( angle * I );
            }

            for ( int k = 0; k < N; k++ ) {
                result[k+j] = handle->tmp[k];
            }

        }

        N *= 2;
        M *= 2;
    }
} /* FFT_iterative_inner() */

void FFT( fft_t* handle, float* samples, double complex* result ) {
    for ( int i = 0; i < handle->n_samples; i++ ) {
        result[i] = creal( samples[i] );
    }

    FFT_iterative_inner( handle, result );
} /* FFT() */

/*---------------------------------------------------------------------------*/

void FFT_recursive_inner( fft_t* handle, double complex* samples, int n_samples, int start ) {

    if ( n_samples <= 1 ) {
        return;
    }

    int M = n_samples / 2;

    int even_idx = 0, odd_idx = M;
    for ( int i = start; i < start+n_samples; i++ ) {
        if ( i % 2 == 0 ) {
            handle->tmp[even_idx++] = samples[i];
        }
        else {
            handle->tmp[odd_idx++]  = samples[i];
        }
    }

    for ( int i = 0; i < n_samples; i++ ) {
        samples[i+start] = handle->tmp[i];
    }

    FFT_recursive_inner( handle, samples, M, start );
    FFT_recursive_inner( handle, samples, M, start+M );


    float angle;

    for ( int k = 0; k < M; k++ ) {
        angle = TWO_PI * ((float)k / (float)n_samples);

        handle->tmp[k]   = samples[start+k] + samples[start+k+M] * cexp( angle );
        handle->tmp[k+M] = samples[start+k] - samples[start+k+M] * cexp( angle );
    }

    for ( int i = 0; i < n_samples; i++ ) {
        samples[i+start] = handle->tmp[i];
    }
} /* FFT_recursive_inner() */

void FFT_recursive( fft_t* handle, float* samples, double complex* result ) {
    for ( int i = 0; i < handle->n_samples; i++ ) {
        result[i] = samples[i];
    }

    FFT_recursive_inner( handle, result, handle->n_samples, 0 );
} /* FFT_recursive() */

/*---------------------------------------------------------------------------*/

void FFT_manipulate( fft_t* handle, double complex* fft_result, double complex* manipulator ) {
    for ( int i = 0; i < handle->n_samples; i++ ) {
        fft_result[i] *= manipulator[i];
    }
} /* FFT_manipulate() */

/*---------------------------------------------------------------------------*/

void IFFT( fft_t* handle, double complex* fft_result, float* samples ) {
    for ( int i = 0; i < handle->n_samples; i++ ) {
        handle->ifft_values[i] = fft_result[i];
    }

    FFT_iterative_inner( handle, handle->ifft_values );

    for ( int i = 0; i < handle->n_samples; i++ ) {
        samples[i] = creal( handle->ifft_values[i] ) / (float)handle->n_samples;
    }
} /* IFFT() */

/*---------------------------------------------------------------------------*/

void FFT_magnitude( fft_t* handle, double complex* fft_result, float* magnitudes ) {
    for ( int i = 0; i < handle->n_samples; i++ ) {
        magnitudes[i] = cabs( fft_result[i] );
    }
} /* FFT_magnitude() */

/*---------------------------------------------------------------------------*/

void FFT_normalize( fft_t* handle, double complex* fft_result ) {
    float M = (float)( handle->n_samples / 2 );
    for ( int i = 0; i < handle->n_samples; i++ ) {
        fft_result[i] /= M;
    }
} /* FFT_normalize() */

/*---------------------------------------------------------------------------*/

void FFT_Close( fft_t* handle ) {
    free( handle->tmp );
    free( handle->ifft_values );
    free( handle );
} /* FFT_Close() */
