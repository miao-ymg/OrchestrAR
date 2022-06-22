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

#include "stdafx.h"
#include "cpct_dll.h"

#include "../../cpct-mstftm/cpct-mstftm/CPCT_MSTFTM.h"
#pragma comment(lib, "../Debug/cpct-mstftm.lib")

using namespace CPCT;


API_CPCT HANDLE createCpctByDefault()
{
	CPCT_MSTFTM *cpct = new CPCT_MSTFTM();
	return (HANDLE)cpct;
}

API_CPCT HANDLE createCpctByParams( int winlen, int hoplen, int nit )
{
	CPCT_MSTFTM *cpct = new CPCT_MSTFTM(winlen, hoplen, nit);
	return (HANDLE)cpct;
}

API_CPCT void setData( HANDLE h, const float* data, int datalength, int nChannels )
{
	CPCT_MSTFTM *cpct = (CPCT_MSTFTM*)h;
	cpct->setData(data, datalength, nChannels);
}

API_CPCT void setParams( HANDLE h, float tempo, float pitch )
{
	CPCT_MSTFTM *cpct = (CPCT_MSTFTM*)h;
	cpct->setParams(tempo, pitch);
}

API_CPCT void getData( HANDLE h, float* data, int& datalength )
{
	CPCT_MSTFTM *cpct = (CPCT_MSTFTM*)h;
	cpct->getData(data, datalength);
}

API_CPCT void destroyCpct( HANDLE h )
{
	CPCT_MSTFTM *cpct = (CPCT_MSTFTM*)h;
	delete cpct;
}




