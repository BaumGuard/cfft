#ifndef FFT_H
#define FFT_H

/*---------------------------------------------------------------------------*/

// FFT handle
typedef struct fft_t fft_t;

// Complex number
typedef struct {
    float re;
    float im;
} Complex;

/*---------------------------------------------------------------------------*/

/*
Initialize the FFT processor

Args:
 - handle       : Pointer to the FFT handle to be Initialized
 - n_samples    : Number of samples per sample block (must be a power of 2!)
*/
void FFT_Init( fft_t** handle, int n_samples );



/*
Perform the Fast Fourier Transform (FFT) on a block of samples
(Convert the time domain into the frequency domain)

Args:
 - handle       : Pointer to the initialized FFT handle
 - samples      : Block of samples in float format (must be of equal length
                  as defined by n_samples in FFT_Init)
 - result       : Result array of complex numbers to write the coefficicents
                  in (Must be of equal length as the array samples)
*/
void FFT( fft_t* handle, float* samples, Complex* result );

/*
Alternative recursive FFT implementation
(Convert the time domain into the frequency domain)

Args:
 - handle       : Pointer to the initialized FFT handle
 - samples      : Block of samples in float format (must be of equal length
                  as defined by n_samples in FFT_Init)
 - result       : Result array of complex numbers to write the coefficicents
                  in (Must be of equal length as the array samples)
*/
void FFT_recursive( fft_t* handle, float* samples, Complex* result );

/*
Manipulate the FFT result (e. g. for filtering) by multiplying the FFT result
by the values in the manipulator array

Args:
 - handle       : Pointer to the initialized FFT handle
 - fft_result   : FFT result array
 - manipulator  : Values to be multiplied with the FFT array
*/
void FFT_manipulate( fft_t* handle, Complex* fft_result, Complex* manipulator );

/*
Get the magnitudes of an FFT result

Args:
 - handle       : Pointer to the initialized FFT handle
 - fft_result   : Result array of an FFT
 - magnitudes   : Pointer to a float array in which the magnitudes should be
                  stored
*/
void FFT_magnitude( fft_t* handle, Complex* fft_result, float* magnitudes );


/*
Inverse Fast Fourier Transform (converts the frequency domain into the time
domain)

Args:
 - handle       : Pointer to the initialized FFT handle
 - fft_result   : FFT result array
 - samples      : Array to store the time domain samples in
*/
void IFFT( fft_t* handle, Complex* fft_result, float* samples );


/*
Normalize the FFT

Args:
 - handle       : Pointer to the initialized FFT handle
 - fft_result   : FFT result array to be normalized
*/
void FFT_normalize( fft_t* handle, Complex* fft_result );

/*
Close the FFT processor

Args:
 - handle       : Pointer to the initialized FFT handle
*/
void FFT_Close( fft_t* handle );

#endif
