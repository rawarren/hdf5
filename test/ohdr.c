/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdf.ncsa.uiuc.edu/HDF5/doc/Copyright.html.  If you do not have     *
 * access to either file, you may request a copy from hdfhelp@ncsa.uiuc.edu. *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Programmer:  Robb Matzke <matzke@llnl.gov>
 *              Tuesday, November 24, 1998
 */
#include "h5test.h"
#include "H5Iprivate.h"
/*
 * This file needs to access private datatypes from the H5O package.
 */
#define H5O_PACKAGE
#include "H5Opkg.h"

/*
 * This file needs to access private datatypes from the H5G package.
 */
#define H5G_PACKAGE
#include "H5Gpkg.h"

const char *FILENAME[] = {
    "ohdr",
    NULL
};

/* The tbogus.h5 is generated from gen_bogus.c in HDF5 'test' directory.
 * To get this data file, define H5O_ENABLE_BOGUS in src/H5Oprivate, rebuild
 * the library and simply compile gen_bogus.c with that HDF5 library and run it. */
#define FILE_BOGUS "tbogus.h5"


/*-------------------------------------------------------------------------
 * Function:	main
 *
 * Purpose:
 *
 * Return:	Success:
 *
 *		Failure:
 *
 * Programmer:	Robb Matzke
 *              Tuesday, November 24, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
int
main(void)
{
    hid_t	fapl=-1, file=-1;
    hid_t	dset=-1;
    H5F_t	*f=NULL;
    char	filename[1024];
    H5O_loc_t	oh_loc;
    time_t	time_new, ro;
    int		i;
    const char  *envval = NULL;
 
    /* Reset library */
    h5_reset();
    fapl = h5_fileaccess();
    h5_fixname(FILENAME[0], fapl, filename, sizeof filename);
    if ((file=H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl))<0)
	goto error;
    if (NULL==(f=H5I_object(file))) {
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }

    /*
     * Test object header creation
     * (using default group creation property list only because it's convenient)
     */
    TESTING("object header creation");
    HDmemset(&oh_loc, 0, sizeof(oh_loc));
    if(H5O_create(f, H5P_DATASET_XFER_DEFAULT, (size_t)64, H5P_GROUP_CREATE_DEFAULT, &oh_loc/*out*/)<0) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    PASSED();

    /* create a new message */
    TESTING("message creation");
    time_new = 11111111;
    if(H5O_msg_create(&oh_loc, H5O_MTIME_NEW_ID, 0, 0, &time_new, H5P_DATASET_XFER_DEFAULT) < 0) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    if (H5AC_flush(f, H5P_DATASET_XFER_DEFAULT, TRUE)<0) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    if (NULL==H5O_msg_read(&oh_loc, H5O_MTIME_NEW_ID, 0, &ro, H5P_DATASET_XFER_DEFAULT)) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    if (ro!=time_new) {
	H5_FAILED();
	HDfprintf(stdout, "    got: {%ld}\n", (long)ro);
	HDfprintf(stdout, "    ans: {%ld}\n", (long)time_new);
	goto error;
    }
    PASSED();

    /*
     * Test modification of an existing message.
     */
    TESTING("message modification");
    time_new = 33333333;
    if (H5O_msg_write(&oh_loc, H5O_MTIME_NEW_ID, 0, 0, 0, &time_new, H5P_DATASET_XFER_DEFAULT)<0) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    if (H5AC_flush(f, H5P_DATASET_XFER_DEFAULT, TRUE)<0) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    if (NULL==H5O_msg_read(&oh_loc, H5O_MTIME_NEW_ID, 0, &ro, H5P_DATASET_XFER_DEFAULT)) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    if (ro!=time_new) {
	H5_FAILED();
	HDfprintf(stdout, "    got: {%ld}\n", (long)ro);
	HDfprintf(stdout, "    ans: {%ld}\n", (long)time_new);
	goto error;
    }
    PASSED();


    /*
     * Test creation of a second message of the same type.
     */
    TESTING("duplicate message creation");
    time_new = 55555555;
    if(H5O_msg_create(&oh_loc, H5O_MTIME_NEW_ID, 0, 0, &time_new, H5P_DATASET_XFER_DEFAULT) < 0) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    if (H5AC_flush(f, H5P_DATASET_XFER_DEFAULT, TRUE)<0) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    if (NULL==H5O_msg_read(&oh_loc, H5O_MTIME_NEW_ID, 1, &ro, H5P_DATASET_XFER_DEFAULT)) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    if (ro!=time_new) {
	H5_FAILED();
	HDfprintf(stdout, "    got: {%ld}\n", (long)ro);
	HDfprintf(stdout, "    ans: {%ld}\n", (long)time_new);
	goto error;
    }
    PASSED();

    /*
     * Test modification of the second message with a symbol table.
     */
    TESTING("duplicate message modification");
    time_new = 77777777;
    if (H5O_msg_write(&oh_loc, H5O_MTIME_NEW_ID, 1, 0, 0, &time_new, H5P_DATASET_XFER_DEFAULT)<0) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    if (H5AC_flush(f, H5P_DATASET_XFER_DEFAULT, TRUE)<0) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    if (NULL==H5O_msg_read(&oh_loc, H5O_MTIME_NEW_ID, 1, &ro, H5P_DATASET_XFER_DEFAULT)) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    if (ro!=time_new) {
	H5_FAILED();
	HDfprintf(stdout, "    got: {%ld}\n", (long)ro);
	HDfprintf(stdout, "    ans: {%ld}\n", (long)time_new);
	goto error;
    }
    PASSED();

    /*
     * Test creation of a bunch of messages one after another to see
     * what happens when the object header overflows in core.
     * (Use 'old' MTIME message here, because it is large enough to be
     *  replaced with a continuation message (the new one is too small)
     *  and the library doesn't understand how to migrate more than one
     *  message from an object header currently - QAK - 10/8/03)
     */
    TESTING("object header overflow in memory");
    for (i=0; i<40; i++) {
        time_new = (i+1)*1000+1;
        if(H5O_msg_create(&oh_loc, H5O_MTIME_ID, 0, 0, &time_new, H5P_DATASET_XFER_DEFAULT) < 0) {
	    H5_FAILED();
            H5Eprint_stack(H5E_DEFAULT, stdout);
	    goto error;
	}
    }
    if (H5AC_flush(f, H5P_DATASET_XFER_DEFAULT, TRUE)<0) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    PASSED();

    /*
     * Test creation of a bunch of messages one after another to see
     * what happens when the object header overflows on disk.
     */
    TESTING("object header overflow on disk");
    for (i=0; i<10; i++) {
        time_new = (i + 1) * 1000 + 10;
        if(H5O_msg_create(&oh_loc, H5O_MTIME_NEW_ID, 0, 0, &time_new, H5P_DATASET_XFER_DEFAULT) < 0) {
	    H5_FAILED();
            H5Eprint_stack(H5E_DEFAULT, stdout);
	    goto error;
	}
        if (H5AC_flush(f, H5P_DATASET_XFER_DEFAULT, TRUE)<0) {
	    H5_FAILED();
            H5Eprint_stack(H5E_DEFAULT, stdout);
	    goto error;
	}
    }
    PASSED();

    /*
     * Delete all time messages.
     */
    TESTING("message deletion");
    if (H5O_msg_remove(&oh_loc, H5O_MTIME_NEW_ID, H5O_ALL, TRUE, H5P_DATASET_XFER_DEFAULT)<0) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    if (H5O_msg_remove(&oh_loc, H5O_MTIME_ID, H5O_ALL, TRUE, H5P_DATASET_XFER_DEFAULT)<0) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    if (H5O_msg_read(&oh_loc, H5O_MTIME_NEW_ID, 0, &ro, H5P_DATASET_XFER_DEFAULT)) {
	H5_FAILED();
	puts("    H5O_msg_read() should have failed but didn't");
	H5Eclear_stack(H5E_DEFAULT);
	goto error;
    }
    if (H5O_msg_read(&oh_loc, H5O_MTIME_ID, 0, &ro, H5P_DATASET_XFER_DEFAULT)) {
	H5_FAILED();
	puts("    H5O_msg_read() should have failed but didn't");
	H5Eclear_stack(H5E_DEFAULT);
	goto error;
    }
    PASSED();


    /* release resources */
    TESTING("object header closing");
    if (H5O_close(&oh_loc)<0) {
	H5_FAILED();
	H5Eprint_stack(H5E_DEFAULT, stdout);
	goto error;
    }
    if (H5Fclose(file)<0) goto error;
    PASSED();

    /* Test reading dataset with undefined object header message */
    TESTING("reading object with unknown header message");
    envval = HDgetenv("HDF5_DRIVER");
    if (envval == NULL) 
        envval = "nomatch";
    if (HDstrcmp(envval, "core") && HDstrcmp(envval, "multi") && HDstrcmp(envval, "split") && HDstrcmp(envval, "family")) {
        {
	    char testfile[512]="";
	    char *srcdir = getenv("srcdir");

	    /* Build path to test file */
	    if (srcdir && ((HDstrlen(srcdir) + HDstrlen(FILE_BOGUS) + 1) < sizeof(testfile))){
		HDstrcpy(testfile, srcdir);
		HDstrcat(testfile, "/");
	    }
	    HDstrcat(testfile, FILE_BOGUS);

	    if ((file=H5Fopen(testfile, H5F_ACC_RDONLY, fapl))<0)
		goto error;

	    /* Open the dataset with the unknown header message (generated with gen_bogus.c) */
	    if((dset=H5Dopen(file,"/Dataset1"))<0)
		goto error;
	    if (H5Dclose(dset)<0) goto error;

	    if (H5Fclose(file)<0) goto error;
	}
	PASSED();
    }
    else 
    {
        SKIPPED();
        puts("   Test not compatible with current Virtual File Driver");
    }

    puts("All object header tests passed.");
    h5_cleanup(FILENAME, fapl);
    return 0;

    error:
        puts("*** TESTS FAILED ***");
        H5E_BEGIN_TRY {
            H5Fclose(file);
        } H5E_END_TRY;
        return 1;
}
