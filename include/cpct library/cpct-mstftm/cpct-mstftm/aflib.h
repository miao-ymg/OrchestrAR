/*
 * Copyright: (C) 1999-2001 Bruce W. Forsberg
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
 *   Bruce Forsberg  forsberg@tns.net
 *
 */

// aflib
// This file contains all of the common definitions used throughout the
// library.

#ifndef _AFLIB_H
#define _AFLIB_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define FALSE 0
#define TRUE  1


enum aflib_data_size
{
   AFLIB_SIZE_UNDEFINED,
   AFLIB_DATA_8S,
   AFLIB_DATA_8U,
   AFLIB_DATA_16S,
   AFLIB_DATA_16U,
   AFLIB_DATA_32S
};

enum aflib_data_endian
{
   AFLIB_ENDIAN_UNDEFINED,
   AFLIB_ENDIAN_LITTLE,
   AFLIB_ENDIAN_BIG
};

enum aflib_data_orientation
{
   AFLIB_ORIENTATION_UNDEFINED,
   AFLIB_SEQUENTIAL,
   AFLIB_INTERLEAVE
};

enum aflibStatus
{
   AFLIB_SUCCESS = 0,
   AFLIB_ERROR_OPEN = 1,
   AFLIB_ERROR_UNSUPPORTED = 2,
   AFLIB_ERROR_INITIALIZATION_FAILURE = 3,
   AFLIB_NOT_FOUND = 4,
   AFLIB_END_OF_FILE = 5,
   AFLIB_NO_DATA = 6,
   AFLIB_CANT_AUTO_TYPE = 7
};

// Currently supported types
enum aflibFileType
{
   AFLIB_AUTO_TYPE,
   AFLIB_DEV_TYPE,
   AFLIB_MPEG_TYPE,
   AFLIB_WAV_TYPE,
   AFLIB_AU_TYPE
};

enum aflibUndoRedo
{
   AFLIB_UNDO_MODE,
   AFLIB_REDO_MODE,
   AFLIB_NONE_MODE
};



#endif

