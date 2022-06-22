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

#ifndef CPCT_MSTFTM_H
#define CPCT_MSTFTM_H

namespace CPCT
{
	class CPCT_MSTFTM
	{
	private:
/// input parameters
		float* dataInput;	// input data
		int nDataInput;		// length of input data	
		int nChannels;		// number of channels

/// control parameters, default: tempo = 0, pitch = 0;
		float tempo;		// used to control the tempo of the audio [-1,1],1 means 2* faster,-1 means 0.5* slower
		float pitch;		// used to control the pitch of the audio [-12,12] + means higher,- means lower

/// output parameters
		float* dataOutput;	// output data
		int nDataOutput;	// length of output data

/// processing parameters, initialized by constructor function
		int winlen;		// length of processing window
		int hoplen;		// length of hop size, overlap size = winlen - hoplen;
		int nit;		// times of iteration in mstftm based signal estimation
		double *hamwin;	// hamming window

	private:
/// helper function
		// function for sum the elements in the array x
		double sum(double* x, int length);
		// change float array[-1, 1] to short array[-32768,32767]
		void float2short(const double* f, short* s, int numElems);
		// change short array[-32768,32767] to float array[-1,1]
		void short2float(const short* s, double* f, int numElems);

/// private function

		// function for creating hamming window
		void hamming(double* win, int length);

		/* 
		//dft function
		// x and y are the real and imaginary arrays
		// dir =  1 gives forward transform
		// dir = -1 gives reverse transform 
		*/
		bool dft(int dir, int m, double* x1, double* y1);

		/* 
		//fft function
		// x and y are the real and imaginary arrays of 2^m points.
		// dir =  1 gives forward transform
		// dir = -1 gives reverse transform 
		*/
		bool fft(int dir, int m, const double* xi, const double* yi, double* xo, double* yo);
		
		// fft function, fft1 is based on class aflibFFT
		void fft1(unsigned NumSamples,
			int InverseTransform,
			const double   *RealIn,
			const double   *ImagIn,
			double   *RealOut,
			double   *ImagOut );

		// resampling function, resample is based on class aflibConverter
		int resample(double factor, 
			int channels, 
			int &inCount, 
			int outCount, 
			short inArray[], 
			short outArray[] );

		//
		void recon(const double* xSTFTM,	// magnitude of processing data
			double* x_res,			// update in each iteration, when iteration ends, it saves the estimation results
			int nit,				// times of iteration	
			const double* win,				// hamming window
			int datalength);		// processing data length, length(xSTFTM) = length(x_res) = length(win) = datalength



		// function for doing time scale modification
		// input params: (class members) dataInput, nDataInput, dataOutput, nDataOutput,tempo
		void tsm();

		// function for doing pitch modification
		// input params: (class members) dataInput, nDataInput, dataOutput, nDataOutput,pitch
		void pm();
		
		void process();

	public:
/// public function
		CPCT_MSTFTM(void);
		~CPCT_MSTFTM(void);
		CPCT_MSTFTM(int winlen, int hoplen, int nit);

		// float* data is the processed sound data
		// int& datalength return the processed data length
		void getData(float* data, int& datalength);
		
		// const float* data is the unprocessed data
		// int datalength is the unprocessed data length
		// int nChannels is the number of channels
		void setData(const float* data, int datalength, int nChannels);
		
		// set the float tempo, float pitch
		void setParams(float tempo, float pitch);

	};
	
}


#endif