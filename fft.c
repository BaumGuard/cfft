#include "fft.h"

#include <stdlib.h>
#include <math.h>

#define TWO_PI 2.0 * M_PI

/*---------------------------------------------------------------------------*/

struct fft_t {
    Complex* tmp;
    Complex* ifft_values;
    int n_samples;
};

/*---------------------------------------------------------------------------*/

void FFT_Init( fft_t** handle, int n_samples ) {
    (*handle) = (fft_t*) malloc( sizeof(fft_t) );
    (*handle)->tmp = (Complex*) calloc( n_samples, sizeof(Complex) );
    (*handle)->ifft_values = (Complex*) calloc( n_samples, sizeof(Complex) );
    (*handle)->n_samples = n_samples;
} /* FFT_Init() */

/*---------------------------------------------------------------------------*/

void FFT_iterative_inner( fft_t* handle, Complex* result ) {
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
                    handle->tmp[odd_idx++] = result[k];
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
                float w_re = cos( angle );
                float w_im = sin( angle );

                handle->tmp[k].re =
                    result[j+k].re +
                    result[j+k+M].re * w_re -
                    result[j+k+M].im * w_im;
                handle->tmp[k].im =
                    result[j+k].im +
                    result[j+k+M].re * w_im +
                    result[j+k+M].im * w_re;

                handle->tmp[k+M].re =
                    result[j+k].re -
                    result[j+k+M].re * w_re +
                    result[j+k+M].im * w_im;
                handle->tmp[k+M].im =
                    result[j+k].im -
                    result[j+k+M].re * w_im -
                    result[j+k+M].im * w_re;

            }

            for ( int k = 0; k < N; k++ ) {
                result[k+j] = handle->tmp[k];
            }

        }

        N *= 2;
        M *= 2;
    }
} /* FFT_iterative_inner() */

void FFT( fft_t* handle, float* samples, Complex* result ) {
    for ( int i = 0; i < handle->n_samples; i++ ) {
        result[i].re = samples[i];
        result[i].im = 0.0;
    }

    FFT_iterative_inner( handle, result );
} /* FFT() */

/*---------------------------------------------------------------------------*/

void FFT_recursive_inner( fft_t* handle, Complex* samples, int n_samples, int start ) {

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
            handle->tmp[odd_idx++] = samples[i];
        }
    }

    for ( int i = 0; i < n_samples; i++ ) {
        samples[i+start] = handle->tmp[i];
    }

    FFT_recursive_inner( handle, samples, M, start );
    FFT_recursive_inner( handle, samples, M, start+M );


    float angle, w_re, w_im;

    for ( int k = 0; k < M; k++ ) {
        angle = TWO_PI * ((float)k / (float)n_samples);
        w_re = cos( angle );
        w_im = sin( angle );

        handle->tmp[k].re =
            samples[start+k].re +
            samples[start+k+M].re * w_re -
            samples[start+k+M].im * w_im;
        handle->tmp[k].im =
            samples[start+k].im +
            samples[start+k+M].re * w_im +
            samples[start+k+M].im * w_re;

        handle->tmp[k+M].re =
            samples[start+k].re -
            samples[start+k+M].re * w_re +
            samples[start+k+M].im * w_im;
        handle->tmp[k+M].im =
            samples[start+k].im -
            samples[start+k+M].re * w_im -
            samples[start+k+M].im * w_re;
    }

    for ( int i = 0; i < n_samples; i++ ) {
        samples[i+start] = handle->tmp[i];
    }
} /* FFT_recursive_inner() */

void FFT_recursive( fft_t* handle, float* samples, Complex* result ) {
    for ( int i = 0; i < handle->n_samples; i++ ) {
        result[i].re = samples[i];
        result[i].im = 0.0;
    }

    FFT_recursive_inner( handle, result, handle->n_samples, 0 );
} /* FFT_recursive() */

/*---------------------------------------------------------------------------*/

void FFT_manipulate( fft_t* handle, Complex* fft_result, Complex* manipulator ) {
    for ( int i = 0; i < handle->n_samples; i++ ) {
        fft_result[i].re *= manipulator[i].re;
        fft_result[i].im *= manipulator[i].im;
    }
} /* FFT_manipulate() */

/*---------------------------------------------------------------------------*/

void IFFT( fft_t* handle, Complex* fft_result, float* samples ) {
    for ( int i = 0; i < handle->n_samples; i++ ) {
        handle->ifft_values[i].re = fft_result[i].re;
        handle->ifft_values[i].im = -fft_result[i].im;
    }

    FFT_iterative_inner( handle, handle->ifft_values );

    for ( int i = 0; i < handle->n_samples; i++ ) {
        samples[i] = handle->ifft_values[i].re / (float)handle->n_samples;
    }
} /* IFFT() */

/*---------------------------------------------------------------------------*/

void FFT_magnitude( fft_t* handle, Complex* fft_result, float* magnitudes ) {
    for ( int i = 0; i < handle->n_samples; i++ ) {
        magnitudes[i] = sqrt(
            fft_result[i].re*fft_result[i].re +
            fft_result[i].im*fft_result[i].im
        );
    }
} /* FFT_magnitude() */

/*---------------------------------------------------------------------------*/

void FFT_normalize( fft_t* handle, Complex* fft_result ) {
    float M = (float)( handle->n_samples / 2 );
    for ( int i = 0; i < handle->n_samples; i++ ) {
        fft_result[i].re /= M;
        fft_result[i].im /= M;
    }
} /* FFT_normalize() */

/*---------------------------------------------------------------------------*/

void FFT_Close( fft_t* handle ) {
    free( handle->tmp );
    free( handle->ifft_values );
    free( handle );
} /* FFT_Close() */
