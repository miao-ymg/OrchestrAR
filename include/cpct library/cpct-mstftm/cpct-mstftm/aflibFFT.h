#ifndef _AFLIBFFT_H
#define _AFLIBFFT_H


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "aflib.h"


typedef struct {
        double re, im;
    } COMPLEX;

/*! \class aflibFFT
    \brief Performs a forward or reverse FFT.

 This class provides a FFT for other classes in this library to use. There is
 only one API for this class and it is fft_double. It will perform both a
 forward and reverse FFT. It operates on doubles.
*/

class aflibFFT {

public:

   aflibFFT();

   ~aflibFFT();

   void
   fft_double (
      unsigned  NumSamples,
      int       InverseTransform,
      const double   *RealIn,
      const double   *ImagIn,
      double   *RealOut,
      double   *ImagOut );

private:

unsigned int Nfactors;
COMPLEX      *W_factors;

   int
   fft (
      COMPLEX *in,
      unsigned  n,
      COMPLEX *out);

   int
   rft (
      COMPLEX *in,
      unsigned  n,
      COMPLEX *out);

   void
   Fourier (
      COMPLEX *in,
      unsigned  n,
      COMPLEX *out);

   unsigned
   radix (unsigned n);

   void
   split (
      register COMPLEX *in,
      register unsigned r,
      register unsigned m,
      register COMPLEX *out);

   void
   join (
      register COMPLEX *in,
      register unsigned m,
      register unsigned n,
      register COMPLEX *out);

   int
   W_init(unsigned n);

};


#endif
