
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "aflibFFT.h"


#define pi   3.1415926535897932384626434

#define		c_re(c)		((c).re)
#define		c_im(c)		((c).im)

/* C_add_mul adds product of c1 and c2 to c.  */
#define	c_add_mul(c, c1, c2)	{ COMPLEX C1, C2; C1 = (c1); C2 = (c2); \
				  c_re (c) += C1.re * C2.re - C1.im * C2.im; \
				  c_im (c) += C1.re * C2.im + C1.im * C2.re; }

/* C_conj substitutes c by its complex conjugate. */
#define c_conj(c)		{ c_im (c) = -c_im (c); }

/* C_realdiv divides complex c by real.  */
#define	c_realdiv(c, real)	{ c_re (c) /= (real); c_im (c) /= (real); }

/*
 * W gives the (already computed) Wn ^ k (= e ^ (2pi * i * k / n)).
 * Notice that the powerseries of Wn has period Nfactors.
 */
#define	W(n, k)		(W_factors [((k) * (Nfactors / (n))) % Nfactors])



/*! \brief Constructor.
*/
aflibFFT::aflibFFT()
{
   Nfactors = 0;
   W_factors = NULL;
}


/*! \brief Destructor.
*/
aflibFFT::~aflibFFT()
{
   if (W_factors != NULL)
      delete W_factors;
}


/*! \brief Performs a forward or reverse FFT.

     This is the main API is this class. It will perform either a forward or
     inverse FFT depending how InverseTransform is set. If set to FALSE then
     forward FFT will be performed, TRUE and a inverse FFT will be performed.
     The number of samlpes (NumSamples) must be a power of 2. The ImagIn
     pointer can be NULL if there are no imaginary values. The user is
     responsable for passing in pointers for RealOut and ImagOut containing
     arrays of the proper size.
*/
void
aflibFFT::fft_double (
    unsigned  NumSamples,
    int       InverseTransform,
    const double   *RealIn,
    const double   *ImagIn,
    double   *RealOut,
    double   *ImagOut )
{

   COMPLEX  in[1024];
   COMPLEX  out[1024];
   COMPLEX  * in_local = NULL;
   COMPLEX  * out_local = NULL;
   register COMPLEX  * in_ptr;
   register COMPLEX  * out_ptr;
   register unsigned int      i;


   // IF 1024 samples or less use local buffer else allocate memory
   if (NumSamples > 1024)
   {
      in_local = new COMPLEX[NumSamples];
      out_local = new COMPLEX[NumSamples];
      in_ptr = in_local;
      out_ptr = out_local;
   }
   else
   {
      in_ptr = in;
      out_ptr = out;
   }

   // Fill real and imaginary array
   for (i = 0; i < NumSamples; i++)
   {
      c_re(in_ptr[i]) = RealIn[i];
      if (ImagIn == NULL)
         c_im(in_ptr[i]) = 0.0;
      else
         c_im(in_ptr[i]) = ImagIn[i];
   }

   // Perform transform
   if (InverseTransform == TRUE)
   {
      rft(in_ptr, NumSamples, out_ptr);
   }
   else
   {
      fft(in_ptr, NumSamples, out_ptr);
   }

   // Fill real and imaginary array
   for (i = 0; i < NumSamples; i++)
   {
      RealOut[i] = c_re(out_ptr[i]);
      ImagOut[i] = c_im(out_ptr[i]);
   }

   // Free memory if local arrays were not used
   if (in_local != NULL)
      delete [] in_local;
   if (out_local != NULL)
      delete [] out_local;
}


/*
 * Forward Fast Fourier Transform on the n samples of complex array in.
 * The result is placed in out.  The number of samples, n, is arbitrary.
 * The W-factors are calculated in advance.
 */
int
aflibFFT::fft (
   COMPLEX *in,
   unsigned  n,
   COMPLEX *out)
{
	unsigned i;

	for (i = 0; i < n; i++)
		c_conj (in [i]);
	
	if (W_init (n) == -1)
		return -1;

	Fourier (in, n, out);

	for (i = 0; i < n; i++) {
		c_conj (out [i]);
		c_realdiv (out [i], n);
	}

	return 0;
}


/*
 * Reverse Fast Fourier Transform on the n complex samples of array in.
 * The result is placed in out.  The number of samples, n, is arbitrary.
 * The W-factors are calculated in advance.
 */
int
aflibFFT::rft (
   COMPLEX *in,
   unsigned  n,
   COMPLEX *out)
{
	if (W_init (n) == -1)
		return -1;

	Fourier (in, n, out);

	return 0;
}


/*
 * Recursive (reverse) complex fast Fourier transform on the n
 * complex samples of array in, with the Cooley-Tukey method.
 * The result is placed in out.  The number of samples, n, is arbitrary.
 * The algorithm costs O (n * (r1 + .. + rk)), where k is the number
 * of factors in the prime-decomposition of n (also the maximum
 * depth of the recursion), and ri is the i-th primefactor.
 */
void
aflibFFT::Fourier (
   COMPLEX *in,
   unsigned  n,
   COMPLEX *out)
{
	unsigned r;

	if ((r = radix (n)) < n)
		split (in, r, n / r, out);
	join (in, n / r, n, out);
}


/*
 * Give smallest possible radix for n samples.
 * Determines (in a rude way) the smallest primefactor of n.
 */
unsigned
aflibFFT::radix (unsigned n)
{
	unsigned r;

	if (n < 2)
		return 1;

	for (r = 2; r < n; r++)
		if (n % r == 0)
			break;
	return r;
}


/*
 * Split array in of r * m samples in r parts of each m samples,
 * such that in [i] goes to out [(i % r) * m + (i / r)].
 * Then call for each part of out Fourier, so the r recursively
 * transformed parts will go back to in.
 */
void
aflibFFT::split (
   register COMPLEX *in,
   register unsigned r,
   register unsigned m,
   register COMPLEX *out)
{
	register unsigned k, s, i, j;

	for (k = 0, j = 0; k < r; k++)
		for (s = 0, i = k; s < m; s++, i += r, j++)
			out [j] = in [i];

	for (k = 0; k < r; k++, out += m, in += m)
		Fourier (out, m, in);
}


/*
 * Sum the n / m parts of each m samples of in to n samples in out.
 * 		   r - 1
 * Out [j] becomes  sum  in [j % m] * W (j * k).  Here in is the k-th
 * 		   k = 0   k	       n		 k
 * part of in (indices k * m ... (k + 1) * m - 1), and r is the radix.
 * For k = 0, a complex multiplication with W (0) is avoided.
 */
void
aflibFFT::join (
   register COMPLEX *in,
   register unsigned m,
   register unsigned n,
   register COMPLEX *out)
{
	register unsigned i, j, jk, s;

	for (s = 0; s < m; s++)
		for (j = s; j < n; j += m) {
			out [j] = in [s];
			for (i = s + m, jk = j; i < n; i += m, jk += j)
				c_add_mul (out [j], in [i], W (n, jk));
		}
}


int
aflibFFT::W_init(unsigned n)
{
    unsigned k;
 
    if (n == Nfactors)
        return 0;
    if (Nfactors != 0 && W_factors != 0)
        delete [] W_factors;
    if ((Nfactors = n) == 0)
        return 0;
    if ((W_factors = new COMPLEX[n]) == NULL)
        return -1;
 
    for (k = 0; k < n; k++) {
        c_re (W_factors [k]) = cos (2 * pi * k / n);
        c_im (W_factors [k]) = sin (2 * pi * k / n);
    }
 
    return 0;
}



