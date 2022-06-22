/*
 * Author        : Copyright (c) Jie Deng
 * Author e-mail : ddqre@163.com
 *
 * License :
 * cpct-mstftm processing library
 * Copyright (c) Jie Deng
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details. 
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 *
 */

#include "CPCT_MSTFTM.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "aflibFFT.h"
#include "aflibConverter.h"

using namespace CPCT;

#define PI 3.1415926535897932384626434


/// --------------------constructor-------------------------
CPCT::CPCT_MSTFTM::CPCT_MSTFTM(void)
{
	this->winlen = 256;
	this->hoplen = winlen / 4;
	this->nit = 5;
	this->hamwin = new double[this->winlen];
	this->hamming(this->hamwin, this->winlen);

	this->tempo = 0;
	this->pitch = 0;

	this->dataInput = NULL;
	this->nDataInput = 0;
	this->dataOutput = NULL;
	this->nDataOutput = 0;
}

CPCT::CPCT_MSTFTM::~CPCT_MSTFTM(void)
{
	delete dataInput;
	delete dataOutput;
	delete hamwin;
}

CPCT::CPCT_MSTFTM::CPCT_MSTFTM( int winlen, int hoplen, int nit )
{
	this->winlen = winlen;
	this->hoplen = hoplen;
	this->nit = nit;
	this->hamwin = new double[this->winlen];
	this->hamming(this->hamwin, this->winlen);

	this->tempo = 0;
	this->pitch = 0;

	this->dataInput = NULL;
	this->nDataInput = 0;
	this->dataOutput = NULL;
	this->nDataOutput = 0;
}


/// ------------------------public function----------------------
void CPCT::CPCT_MSTFTM::getData( float* data, int& datalength )
{
	this->process();

	if (this->nChannels == 2)
	{
		this->nDataOutput *= 2;
		float* temp = this->dataOutput;
		this->dataOutput = new float[this->nDataOutput];
		for (int i = 0; i < this->nDataOutput; i+=2)
		{
			this->dataOutput[i] = temp[i / 2];
			this->dataOutput[i + 1] = temp[i / 2];
		}
		delete temp;
	}

	datalength = this->nDataOutput;
	for (int i = 0; i < this->nDataOutput; i++)
	{
		data[i] = this->dataOutput[i];
	}
	
}


void CPCT::CPCT_MSTFTM::setData( const float* data, int datalength, int nChannels )
{
	if (this->dataInput != NULL)
	{
		delete this->dataInput;
	}

	this->nChannels = nChannels;
	if (this->nChannels == 1)
	{
		this->nDataInput = datalength;
		this->dataInput = new float[this->nDataInput];
		for (int i = 0; i < datalength; i++)
		{
			this->dataInput[i] = data[i];
		}
	} 
	else
	{
		this->nDataInput = datalength / 2;
		this->dataInput = new float[this->nDataInput];
		for (int i = 0; i < datalength; i+=2)
		{
			this->dataInput[i / 2] = data[i];	
		}
	}
	
}


void CPCT::CPCT_MSTFTM::setParams( float tempo, float pitch )
{
	// set tempo
	if (tempo > 1)
	{
		this->tempo = 1;
	}
	else if (tempo < -1)
	{
		this->tempo = -1;
	}
	else
	{
		this->tempo = tempo;
	}
	
	// set pitch
	if (pitch > 12)
	{
		this->pitch = 12;
	} 
	else if (pitch < -12)
	{
		this->pitch = -12;
	}
	else
	{
		this->pitch = pitch;
	}
}

///------------------helper function------------------------------

double CPCT::CPCT_MSTFTM::sum( double* x, int length )
{	
	double result = 0; 
	for (int i = 0; i < length; i++)
	{
		result += x[i];
	}
	return result;
}

void CPCT::CPCT_MSTFTM::float2short( const double* f, short* s, int numElems )
{
	int iTemp;
	for (int i = 0; i < numElems; i++)
	{
		iTemp = (int)(32768.0 * f[i]);
		if (iTemp < -32768) iTemp = -32768;
		if (iTemp > 32767) iTemp = 32767;
		s[i] = (short)iTemp;
	}
}

void CPCT::CPCT_MSTFTM::short2float( const short* s, double* f, int numElems )
{
	double iTemp;
	for (int i = 0; i < numElems; i++)
	{
		iTemp = s[i] / 32768.0;
		f[i] = iTemp;
	}
}

///----------------- private function------------------------------
void CPCT::CPCT_MSTFTM::hamming( double* win, int length )
{
	for (int n = 0; n < length; n++)
	{
		win[n] = 0.54 - 0.46 * cos(2 * PI * n / (length - 1));
	}

}


bool CPCT::CPCT_MSTFTM::dft( int dir, int m, double* x1, double* y1 )
{
	long i, k;
	double arg;
	double cosarg, sinarg;
	double *x2=NULL,*y2=NULL;
	

	x2 = (double *)malloc(m * sizeof(double));
	y2 = (double *)malloc(m * sizeof(double));
	if (x2 == NULL || y2 == NULL)
		return false;

	for (i = 0; i < m; i++) 
	{
		x2[i] = 0;
		y2[i] = 0;
		arg = - dir * 2.0 * 3.141592654 * (double)i / (double)m;
		for (k=0;k<m;k++) 
		{
			cosarg = cos(k * arg);
			sinarg = sin(k * arg);
			x2[i] += (x1[k] * cosarg - y1[k] * sinarg);
			y2[i] += (x1[k] * sinarg + y1[k] * cosarg);
		}
	}

	/* Copy the data back */
	if (dir == 1) 
	{
		for (i=0;i<m;i++) 
		{
			x1[i] = x2[i] / (double)m;
			y1[i] = y2[i] / (double)m;
		}
	}
	else 
	{
		for (i=0;i<m;i++) 
		{
			x1[i] = x2[i];
			y1[i] = y2[i];
		}
	}

	free(x2);
	free(y2);
	return true;
}


bool CPCT::CPCT_MSTFTM::fft( int dir, int m, const double* xi, const double* yi, double* xo, double* yo )
{
	long n,i,i1,j,k,i2,l,l1,l2;
	double c1,c2,tx,ty,t1,t2,u1,u2,z;
	for (int num = 0; num < m; num++)
	{
		xo[num] = xi[num];
		yo[num] = yi[num];
	}

	/* Calculate the number of points */
	n = 1;
	for (i=0;i<m;i++) 
		n *= 2;

	/* Do the bit reversal */
	i2 = n >> 1;
	j = 0;
	for (i=0;i<n-1;i++) {
		if (i < j) {
			tx = xo[i];
			ty = yo[i];
			xo[i] = xo[j];
			yo[i] = yo[j];
			xo[j] = tx;
			yo[j] = ty;
		}
		k = i2;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}

	/* Compute the FFT */
	c1 = -1.0; 
	c2 = 0.0;
	l2 = 1;
	for (l=0;l<m;l++) {
		l1 = l2;
		l2 <<= 1;
		u1 = 1.0; 
		u2 = 0.0;
		for (j=0;j<l1;j++) {
			for (i=j;i<n;i+=l2) {
				i1 = i + l1;
				t1 = u1 * xo[i1] - u2 * yo[i1];
				t2 = u1 * yo[i1] + u2 * xo[i1];
				xo[i1] = xo[i] - t1; 
				yo[i1] = yo[i] - t2;
				xo[i] += t1;
				yo[i] += t2;
			}
			z =  u1 * c1 - u2 * c2;
			u2 = u1 * c2 + u2 * c1;
			u1 = z;
		}
		c2 = sqrt((1.0 - c1) / 2.0);
		if (dir == 1) 
			c2 = -c2;
		c1 = sqrt((1.0 + c1) / 2.0);
	}

	/* Scaling for forward transform */
	if (dir == 1) {
		for (i=0;i<n;i++) {
			xo[i] /= n;
			yo[i] /= n;
		}
	}

	return true;
}


void CPCT::CPCT_MSTFTM::fft1( unsigned NumSamples, 
							 int InverseTransform, 
							 const double *RealIn, 
							 const double *ImagIn, 
							 double *RealOut, 
							 double *ImagOut )
{
	aflibFFT *f = new aflibFFT();
	f->fft_double(NumSamples, InverseTransform, RealIn, ImagIn, RealOut, ImagOut);
	delete f;
}

int CPCT::CPCT_MSTFTM::resample( double factor, 
								int channels, 
								int &inCount, 
								int outCount, 
								short inArray[], 
								short outArray[] )
{
	aflibConverter *c = new aflibConverter(true, true, true);
	int a;
	c->initialize(factor, channels, 1.000000);
	a = c->resample(inCount, outCount, inArray, outArray);
	delete c;
	return a;
}

void CPCT::CPCT_MSTFTM::recon( const double* xSTFTM, 
							  double* x_res, 
							  int nit, 
							  const double* win,
							  int datalength )
{
	double *phi = new double[datalength];
	double *temp = new double[datalength];
	//double *zero = new double[datalength];
	//for (int i = 0; i < datalength; i++)
	//{
	//	zero[i] = 0;
	//}
	double *tempRealOut = new double[datalength];
	double *tempImagOut = new double[datalength];
	double *tempRealOut2 = new double[datalength];
	double *tempImagOut2 = new double[datalength];
	double *rubbish =new double[datalength];

	for (int i = 1; i <= nit; i++)
	{
		for (int j = 0; j < datalength; j++)
		{
			temp[j] = win[j] * x_res[j];
		}

		this->fft1(datalength, FALSE, temp, NULL, tempRealOut, tempImagOut);
		
		for (int j = 0; j < datalength; j++)
		{
			if (tempRealOut[j] == 0)
			{
				if (tempImagOut[j] > 0) 
				{
					phi[j] = PI / 2;
				}
				else
				{
					phi[j] = PI / -2;
				}
				
			}
			if (tempRealOut[j] > 0)
			{
				phi[j] = atan(tempImagOut[j] / tempRealOut[j]);
			}
			if (tempRealOut[j] < 0)
			{
				if (tempImagOut[j] < 0)
				{
					phi[j] = atan(tempImagOut[j] / tempRealOut[j]) - PI;
				}
				else
				{
					phi[j] = PI + atan(tempImagOut[j] / tempRealOut[j]);
				}
			}
			 
		}

		for (int j = 0; j < datalength; j++)
		{
			tempRealOut2[j] = xSTFTM[j] * cos(phi[j]);
			tempImagOut2[j] = xSTFTM[j] * sin(phi[j]);
		}

		this->fft1(datalength, TRUE, tempRealOut2, tempImagOut2, x_res, rubbish);
	}

	delete phi;
	delete temp;
	//delete zero;
	delete tempRealOut;
	delete tempImagOut;
	delete tempRealOut2;
	delete tempImagOut2;
	delete rubbish;
}




void CPCT::CPCT_MSTFTM::tsm()
{
	if (this->dataOutput != NULL)
	{
		delete this->dataOutput;
	}
	// initialize the output data
	this->nDataOutput = (int)(this->nDataInput / pow(2, tempo));
	this->dataOutput = new float[this->nDataOutput];
	for (int i = 0; i < this->nDataOutput; i++)
	{
		this->dataOutput[i] = 0;
	}

	int m_S = (int)(this->hoplen / pow(2,tempo));
	int overlap = this->winlen - this->hoplen;							// length of overlap
	int nframe = (int)((this->nDataInput - overlap) / this->hoplen);	// numbers of frames
	double U = this->sum(this->hamwin, this->winlen) / m_S;

	int k = 0, kk = 0;
	//double *zero = new double[this->winlen];
	//for (int i = 0; i < this->winlen; i++)
	//{
	//	zero[i] = 0;
	//}
	double *frm = new double[this->winlen];	//	data of each frame
	double *RealOut = new double[this->winlen]; 	//  Fourier Transform output, real part
	double *ImagOut = new double[this->winlen];		//  Fourier Transform output, image part
	double *xSTFTM = new double[this->winlen];
	double *res =new double[this->winlen];
	

	for (int nfrm = 1; nfrm <= nframe; nfrm++)
	{		
		for (int i = 0; i < this->winlen; i++)
		{
			frm[i] = this->hamwin[i] * this->dataInput[k + i] / U; 
		}
		
		this->fft1(this->winlen, FALSE, frm, NULL, RealOut, ImagOut);

		for (int i = 0; i < this->winlen; i++)
		{
			xSTFTM[i] = sqrt(pow(RealOut[i], 2) + pow(ImagOut[i], 2));
		}

		if (kk + this->winlen - 1< this->nDataOutput)
		{
			for (int i = 0; i < this->winlen; i++)
			{
				res[i] = dataOutput[kk + i];
			}
		}
		else
		{
			for (int i = 0; i < this->nDataOutput - kk; i++)
			{
				res[i] = dataOutput[kk +i];
			}
			for (int i = this->nDataOutput - kk; i < this->winlen; i++)
			{
				res[i] = 0;
			}
		}

		this->recon(xSTFTM, res, this->nit, this->hamwin, this->winlen);

		if (kk + this->winlen -1 < this->nDataOutput)
		{
			for (int i = 0; i < this->winlen; i++)
			{
				dataOutput[kk + i] += (float)res[i];
			}
		} 
		else
		{
			for (int i = 0; i < this->nDataOutput - kk; i++)
			{
				dataOutput[kk + i] += (float)res[i]; 
			}
		}
		k += this->hoplen;
		kk += m_S;		
	}
	/*
	
	// copy dataOutput to dataInput for later processing
	this->nDataInput = this->nDataOutput;
	for (int i = 0;i < this->nDataInput; i++)
	{
		this->dataInput[i] = this->dataOutput[i];
	}
	// ---------------------------------------------
	*/
	//delete zero;
	delete frm;
	delete RealOut;
	delete ImagOut;
	delete xSTFTM;
	delete res;

}

void CPCT::CPCT_MSTFTM::pm()
{
	if (this->dataOutput != NULL)
	{
		delete this->dataOutput;
	}
	this->nDataOutput = this->nDataInput;
	this->dataOutput = new float[this->nDataOutput];
	for (int i = 0; i < this->nDataOutput; i++)
	{
		this->dataOutput[i] = 0;
	}

	double scale = pow(2, this->pitch / 12);
	int overlap = this->winlen - this->hoplen;
	int nframe = (int) ((this->nDataInput - overlap) / this->hoplen);
	int Lq = (int) (this->winlen * scale);
	int Lq1 = Lq;
	// Lq = Lq + 50;
	double factor = (double)this->winlen / (double)Lq1;
	double U = this->sum(this->hamwin, this->winlen) / this->hoplen;

	double *winq = new double[Lq];
	this->hamming(winq, Lq);
	double *frm = new double[Lq];
	double *frm_resamp = new double[this->winlen];
	int inCount = Lq;
	short *inArray = new short[Lq];
	int outCount = this->winlen;
	short *outArray = new short[this->winlen];
	double *RealOut = new double[this->winlen]; 	//  Fourier Transform output, real part
	double *ImagOut = new double[this->winlen];		//  Fourier Transform output, image part
	double *xSTFTM = new double[this->winlen];
	double *res =new double[this->winlen];


	int k = 0;
	for (int nfrm = 1; nfrm <= nframe; nfrm++)
	{
		if (k + Lq -1 < this->nDataInput)
		{
			for (int i = 0; i < Lq; i++)
			{
				frm[i] = winq[i] * this->dataInput[k + i] / U;
			}
		} 
		else
		{
			for (int i = 0; i < this->nDataInput - k; i++)
			{
				frm[i] = winq[i] * this->dataInput[k + i] / U;
			}
			for (int i = this->nDataInput - k; i < Lq; i++)
			{
				frm[i] = 0;
			}
		}

		this->float2short(frm, inArray, Lq);
		this->resample(factor, this->nChannels, inCount, outCount, inArray, outArray);
		this->short2float(outArray, frm_resamp, this->winlen);
		
		this->fft1(this->winlen, false, frm_resamp, NULL, RealOut, ImagOut);
		for (int i = 0; i < this->winlen; i++ )
		{
			xSTFTM[i] = sqrt(pow(RealOut[i], 2) + pow(ImagOut[i], 2));
		}

		if (k + this->winlen - 1< this->nDataOutput)
		{
			for (int i = 0; i < this->winlen; i++)
			{
				res[i] = dataOutput[k + i];
			}
		}
		else
		{
			for (int i = 0; i < this->nDataOutput - k; i++)
			{
				res[i] = dataOutput[k +i];
			}
			for (int i = this->nDataOutput - k; i < this->winlen; i++)
			{
				res[i] = 0;
			}
		}
		
		this->recon(xSTFTM, res, this->nit, this->hamwin, this->winlen);

		if (k + this->winlen -1 < this->nDataOutput)
		{
			for (int i = 0; i < this->winlen; i++)
			{
				dataOutput[k + i] += (float)res[i];
			}
		} 
		else
		{
			for (int i = 0; i < this->nDataOutput - k; i++)
			{
				dataOutput[k + i] += (float)res[i]; 
			}
		}

		k += this->hoplen;
	}
	/*
	
	// copy dataOutput to dataInput for later processing
	this->nDataInput = this->nDataOutput;
	for (int i = 0;i < this->nDataInput; i++)
	{
		this->dataInput[i] = this->dataOutput[i];
	}
	// ---------------------------------------------
	*/
	delete winq;
	delete frm;
	delete frm_resamp;
	delete inArray;
	delete outArray;
	delete RealOut;
	delete ImagOut;
	delete xSTFTM;
	delete res;

}

void CPCT::CPCT_MSTFTM::process()
{	
	if (this->tempo != 0 && this->pitch == 0)
		this->tsm();
	
	if (this->pitch != 0 && this->tempo == 0)		
		this->pm();	

	if (this->tempo == 0 && this->pitch == 0)
	{
		this->nDataOutput = this->nDataInput;
		this->dataOutput = new float[this->nDataOutput];
		for (int i = 0; i < this->nDataOutput; i++)
		{
			this->dataOutput[i] = this->dataInput[i];
		}
	}

	if (this->tempo != 0 && this->pitch != 0)
	{
		this->tsm();
		delete this->dataInput;
		this->nDataInput = this->nDataOutput;
		this->dataInput = new float[this->nDataInput];
		for (int i = 0; i< this->nDataInput; i++)
		{
			this->dataInput[i] = this->dataOutput[i];
		}
		this->pm();
	}
}





