/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the COPYING file, which can be found at the root of the source code       *
 * distribution tree, or in https://support.hdfgroup.org/ftp/HDF5/releases.  *
 * If you do not have access to either file, you may request a copy from     *
 * help@hdfgroup.org.                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * Programmer:  Quincey Koziol <koziol@hdfgroup.org>
 *              Wednesday, March 17, 2010
 *
 * Purpose:     srcdir querying support.
 */
#ifndef _H5SRCDIR_H
#define _H5SRCDIR_H

/* Include the header file with the correct relative path for the srcdir string */
#include "H5srcdir_str.h"

/* Buffer to construct path in and return pointer to */
extern char srcdir_path[1024];

/* Buffer to construct file in and return pointer to */
extern char srcdir_testpath[1024];

/* Just return the srcdir path */
const char *H5_get_srcdir(void);

/* Append the test file name to the srcdir path and return the whole string */
const char *H5_get_srcdir_filename(const char *);

#endif /* _H5SRCDIR_H */

