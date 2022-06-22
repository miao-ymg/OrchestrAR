/*
 * Author        : Copyright (c) Jie Deng
 * Author e-mail : ddqre@163.com
 *
 * License :
 * cpct_dll processing library
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

#ifndef CPCT_DLL_H
#define CPCT_DLL_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef DLL_EXPORTS
#define API_CPCT _declspec(dllexport)
#else 
#define API_CPCT _declspec(dllimport)
#endif

typedef void* HANDLE;

// create cpct-mstftm by default parameters
API_CPCT HANDLE createCpctByDefault();

// create cpct-mstftm by specific parameters
API_CPCT HANDLE createCpctByParams(int winlen, int hoplen, int nit);

// const float* data is the unprocessed data
// int datalength is the unprocessed data length
// int nChannels is the number of channels
API_CPCT void setData(HANDLE h, const float* data, int datalength, int nChannels);

// // set the tempo and pitch
API_CPCT void setParams(HANDLE h, float tempo, float pitch);

// float* data is the processed sound data
// int& datalength return the processed data length
API_CPCT void getData(HANDLE h, float* data, int& datalength);

API_CPCT void destroyCpct(HANDLE h);

#ifdef __cplusplus
}
#endif

#endif