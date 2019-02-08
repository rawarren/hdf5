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
 * Programmer:  Robb Matzke <matzke@llnl.gov>
 *              Thursday, July 29, 1999
 *
 * Purpose:	This is the MPI-2 I/O driver.
 *
 */

#include "H5FDdrvr_module.h" /* This source code file is part of the H5FD driver module */


#include "H5private.h"		/* Generic Functions			*/
#include "H5CXprivate.h"        /* API Contexts                         */
#include "H5Dprivate.h"		/* Dataset functions			*/
#include "H5Eprivate.h"		/* Error handling		  	*/
#include "H5Fprivate.h"		/* File access				*/
#include "H5FDprivate.h"	/* File drivers				*/
#include "H5FDmpi.h"            /* MPI-based file drivers		*/
#include "H5Iprivate.h"		/* IDs			  		*/
#include "H5MMprivate.h"	/* Memory management			*/
#include "H5Pprivate.h"         /* Property lists                       */

#ifdef H5_HAVE_PARALLEL

/*
 * The driver identification number, initialized at runtime if H5_HAVE_PARALLEL
 * is defined. This allows applications to still have the H5FD_MPIO
 * "constants" in their source code.
 */
static hid_t H5FD_MPIO_g = 0;

/* Whether to allow collective I/O operations */
/* (Value can be set from environment variable also) */
hbool_t H5FD_mpi_opt_types_g = TRUE;

/*
 * The view is set to this value
 */
static char H5FD_mpi_native_g[] = "native";

/*
 * The description of a file belonging to this driver.
 * The EOF value is only used just after the file is opened in order for the
 * library to determine whether the file is empty, truncated, or okay. The MPIO
 * driver doesn't bother to keep it updated since it's an expensive operation.
 */
typedef struct H5FD_mpio_t {
    H5FD_t	pub;		/*public stuff, must be first		*/

/* For comparisons */
#ifndef H5_HAVE_WIN32_API
    /* On most systems the combination of device and i-node number uniquely
     * identify a file.  Note that Cygwin, MinGW and other Windows POSIX
     * environments have the stat function (which fakes inodes)
     * and will use the 'device + inodes' scheme as opposed to the
     * Windows code further below.
     */
    dev_t       device;     /* file device number   */
    ino_t       inode;      /* file i-node number   */

#else
    DWORD       nFileIndexLow;
    DWORD       nFileIndexHigh;
    DWORD       dwVolumeSerialNumber;

#endif
    MPI_File	f;		/*MPIO file handle			*/
    MPI_Comm	comm;		/*communicator				*/
    MPI_Info	info;		/*file information			*/
    int         mpi_rank;       /* This process's rank                  */
    int         mpi_size;       /* Total number of processes            */
    haddr_t	eof;		/*end-of-file marker			*/
    haddr_t	eoa;		/*end-of-address marker			*/
    haddr_t	last_eoa;	/* Last known end-of-address marker	*/
    haddr_t	local_eof;	/* Local end-of-file address for each process */
} H5FD_mpio_t;

/* Private Prototypes */

/* Callbacks */
static herr_t H5FD__mpio_term(void);
static void *H5FD__mpio_fapl_get(H5FD_t *_file);
static void *H5FD__mpio_fapl_copy(const void *_old_fa);
static herr_t H5FD__mpio_fapl_free(void *_fa);
static H5FD_t *H5FD__mpio_open(const char *name, unsigned flags, hid_t fapl_id,
                  haddr_t maxaddr);
static herr_t H5FD__mpio_close(H5FD_t *_file);
static int H5FD__mpio_cmp(const H5FD_t *_f1, const H5FD_t *_f2);
static herr_t H5FD__mpio_query(const H5FD_t *_f1, unsigned long *flags);
static haddr_t H5FD__mpio_get_eoa(const H5FD_t *_file, H5FD_mem_t type);
static herr_t H5FD__mpio_set_eoa(H5FD_t *_file, H5FD_mem_t type, haddr_t addr);
static haddr_t H5FD__mpio_get_eof(const H5FD_t *_file, H5FD_mem_t type);
static herr_t H5FD__mpio_get_handle(H5FD_t *_file, hid_t fapl, void** file_handle);
static herr_t H5FD__mpio_read(H5FD_t *_file, H5FD_mem_t type, hid_t dxpl_id, haddr_t addr,
            size_t size, void *buf);
static herr_t H5FD__mpio_write(H5FD_t *_file, H5FD_mem_t type, hid_t dxpl_id, haddr_t addr,
            size_t size, const void *buf);
static herr_t H5FD__mpio_flush(H5FD_t *_file, hid_t dxpl_id, hbool_t closing);
static herr_t H5FD__mpio_truncate(H5FD_t *_file, hid_t dxpl_id, hbool_t closing);
static int H5FD__mpio_mpi_rank(const H5FD_t *_file);
static int H5FD__mpio_mpi_size(const H5FD_t *_file);
static MPI_Comm H5FD__mpio_communicator(const H5FD_t *_file);
static herr_t H5FD__mpio_get_info(H5FD_t *_file, void** mpi_info);

/* The MPIO file driver information */
static const H5FD_class_mpi_t H5FD_mpio_g = {
    {   /* Start of superclass information */
    "mpio",					/*name			*/
    HADDR_MAX,					/*maxaddr		*/
    H5F_CLOSE_SEMI,				/*fc_degree		*/
    H5FD__mpio_term,                            /*terminate             */
    NULL,					/*sb_size		*/
    NULL,					/*sb_encode		*/
    NULL,					/*sb_decode		*/
    sizeof(H5FD_mpio_fapl_t),			/*fapl_size		*/
    H5FD__mpio_fapl_get,			/*fapl_get		*/
    H5FD__mpio_fapl_copy,			/*fapl_copy		*/
    H5FD__mpio_fapl_free, 			/*fapl_free		*/
    0,		                		/*dxpl_size		*/
    NULL,					/*dxpl_copy		*/
    NULL,					/*dxpl_free		*/
    H5FD__mpio_open,				/*open			*/
    H5FD__mpio_close,				/*close			*/
    H5FD__mpio_cmp,                             /*cmp                   */
    H5FD__mpio_query,		                /*query			*/
    NULL,					/*get_type_map		*/
    NULL,					/*alloc			*/
    NULL,					/*free			*/
    H5FD__mpio_get_eoa,				/*get_eoa		*/
    H5FD__mpio_set_eoa, 			/*set_eoa		*/
    H5FD__mpio_get_eof,				/*get_eof		*/
    H5FD__mpio_get_handle,                      /*get_handle            */
    H5FD__mpio_read,				/*read			*/
    H5FD__mpio_write,				/*write			*/
    H5FD__mpio_flush,				/*flush			*/
    H5FD__mpio_truncate,			/*truncate		*/
    NULL,                                       /*lock                  */
    NULL,                                       /*unlock                */
    H5FD_FLMAP_DICHOTOMY                        /*fl_map                */
    },  /* End of superclass information */
    H5FD__mpio_mpi_rank,                        /*get_rank              */
    H5FD__mpio_mpi_size,                        /*get_size              */
    H5FD__mpio_communicator,                    /*get_comm              */
    H5FD__mpio_get_info                         /*get_info              */
};

#ifdef H5FDmpio_DEBUG
/* Flags to control debug actions in H5Fmpio.
 * Meant to be indexed by characters.
 *
 * 'c' show result of MPI_Get_count after read
 * 'r' show read offset and size
 * 't' trace function entry and exit
 * 'w' show write offset and size
 */
static int H5FD_mpio_Debug[256] =
        { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
#endif


/*--------------------------------------------------------------------------
NAME
   H5FD__init_package -- Initialize interface-specific information

USAGE
    herr_t H5FD__init_package()

RETURNS
    Non-negative on success/Negative on failure

DESCRIPTION
    Initializes any interface-specific data or routines.  (Just calls
    H5FD_mpio_init currently).

--------------------------------------------------------------------------*/
static herr_t
H5FD__init_package(void)
{
    herr_t ret_value = SUCCEED;

    FUNC_ENTER_STATIC

    if(H5FD_mpio_init() < 0)
        HGOTO_ERROR(H5E_VFL, H5E_CANTINIT, FAIL, "unable to initialize mpio VFD")

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* H5FD__init_package() */


/*-------------------------------------------------------------------------
 * Function:    H5FD_mpio_init
 *
 * Purpose:     Initialize this driver by registering the driver with the
 *              library.
 *
 * Return:      Success:    The driver ID for the mpio driver
 *              Failure:    H5I_INVALID_HID
 *
 * Programmer:  Robb Matzke
 *              Thursday, August 5, 1999
 *
 *-------------------------------------------------------------------------
 */
hid_t
H5FD_mpio_init(void)
{
#ifdef H5FDmpio_DEBUG
    static int H5FD_mpio_Debug_inited = 0;
#endif /* H5FDmpio_DEBUG */
    const char *s;                          /* String for environment variables */
    hid_t ret_value = H5I_INVALID_HID;      /* Return value */

    FUNC_ENTER_NOAPI(H5I_INVALID_HID)

    /* Register the MPI-IO VFD, if it isn't already */
    if(H5I_VFL != H5I_get_type(H5FD_MPIO_g))
        H5FD_MPIO_g = H5FD_register((const H5FD_class_t *)&H5FD_mpio_g, sizeof(H5FD_class_mpi_t), FALSE);

    /* Allow MPI buf-and-file-type optimizations? */
    s = HDgetenv("HDF5_MPI_OPT_TYPES");
    if(s && HDisdigit(*s)) {
        long env_val = HDstrtol(s, NULL, 0);
        H5FD_mpi_opt_types_g = (0 == env_val) ? FALSE : TRUE;
    }

#ifdef H5FDmpio_DEBUG
    if(!H5FD_mpio_Debug_inited) {
        /* Retrieve MPI-IO debugging environment variable */
        s = HDgetenv("H5FD_mpio_Debug");
        if(s) {
            /* Set debug mask */
            while(*s) {
                H5FD_mpio_Debug[(int)*s]++;
                s++;
            } /* end while */
        } /* end if */
        H5FD_mpio_Debug_inited++;
    } /* end if */
#endif /* H5FDmpio_DEBUG */

    /* Set return value */
    ret_value = H5FD_MPIO_g;

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5FD_mpio_init() */


/*---------------------------------------------------------------------------
 * Function:	H5FD__mpio_term
 *
 * Purpose:	Shut down the VFD
 *
 * Returns:     Non-negative on success or negative on failure
 *
 * Programmer:  Quincey Koziol
 *              Friday, Jan 30, 2004
 *
 *---------------------------------------------------------------------------
 */
static herr_t
H5FD__mpio_term(void)
{
    FUNC_ENTER_STATIC

    /* Reset VFL ID */
    H5FD_MPIO_g = 0;

    FUNC_LEAVE_NOAPI(SUCCEED)
} /* end H5FD__mpio_term() */


/*-------------------------------------------------------------------------
 * Function:	H5Pset_fapl_mpio
 *
 * Purpose:	Store the user supplied MPIO communicator comm and info in
 *		the file access property list FAPL_ID which can then be used
 *		to create and/or open the file.  This function is available
 *		only in the parallel HDF5 library and is not collective.
 *
 *		comm is the MPI communicator to be used for file open as
 *		defined in MPI_FILE_OPEN of MPI-2. This function makes a
 *		duplicate of comm. Any modification to comm after this function
 *		call returns has no effect on the access property list.
 *
 *		info is the MPI Info object to be used for file open as
 *		defined in MPI_FILE_OPEN of MPI-2. This function makes a
 *		duplicate of info. Any modification to info after this
 *		function call returns has no effect on the access property
 *		list.
 *
 *              If fapl_id has previously set comm and info values, they
 *              will be replaced and the old communicator and Info object
 *              are freed.
 *
 * Return:	Success:	Non-negative
 * 		Failure:	Negative
 *
 * Programmer:	Albert Cheng
 *		Feb 3, 1998
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_fapl_mpio(hid_t fapl_id, MPI_Comm comm, MPI_Info info)
{
    H5FD_mpio_fapl_t	fa;
    H5P_genplist_t *plist;      /* Property list pointer */
    herr_t ret_value;

    FUNC_ENTER_API(FAIL)
    H5TRACE3("e", "iMcMi", fapl_id, comm, info);

    /* Check arguments */
    if(fapl_id == H5P_DEFAULT)
        HGOTO_ERROR(H5E_PLIST, H5E_BADVALUE, FAIL, "can't set values in default property list")
    if(NULL == (plist = H5P_object_verify(fapl_id, H5P_FILE_ACCESS)))
        HGOTO_ERROR(H5E_PLIST, H5E_BADTYPE, FAIL, "not a file access list")
    if(MPI_COMM_NULL == comm)
	HGOTO_ERROR(H5E_PLIST, H5E_BADTYPE, FAIL, "not a valid communicator")

    /* Initialize driver specific properties */
    HDmemset(&fa, 0, sizeof(H5FD_mpio_fapl_t));
    fa.comm = comm;
    fa.info = info;

    /* duplication is done during driver setting. */
    ret_value = H5P_set_driver(plist, H5FD_MPIO, &fa);

done:
    FUNC_LEAVE_API(ret_value)
} /* H5Pset_fapl_mpio() */


/*-------------------------------------------------------------------------
 * Function:	H5Pget_fapl_mpio
 *
 * Purpose:	If the file access property list is set to the H5FD_MPIO
 *		driver then this function returns duplicates of the MPI
 *		communicator and Info object stored through the comm and
 *		info pointers.  It is the responsibility of the application
 *		to free the returned communicator and Info object.
 *
 * Return:	Success:	Non-negative with the communicator and
 *				Info object returned through the comm and
 *				info arguments if non-null. Since they are
 *				duplicates of the stored objects, future
 *				modifications to the access property list do
 *				not affect them and it is the responsibility
 *				of the application to free them.
 * 		Failure:	Negative
 *
 * Programmer:	Robb Matzke
 *		Thursday, February 26, 1998
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_fapl_mpio(hid_t fapl_id, MPI_Comm *comm/*out*/, MPI_Info *info/*out*/)
{
    H5P_genplist_t *plist;      /* Property list pointer */
    const H5FD_mpio_fapl_t *fa;       /* MPIO fapl info */
    MPI_Comm	comm_tmp = MPI_COMM_NULL;
    hbool_t     comm_copied = FALSE;    /* MPI Comm has been duplicated */
    int		mpi_code;		/* MPI return code */
    herr_t      ret_value = SUCCEED;    /* Return value */

    FUNC_ENTER_API(FAIL)
    H5TRACE3("e", "ixx", fapl_id, comm, info);

    /* Check arguments */
    if(NULL == (plist = H5P_object_verify(fapl_id, H5P_FILE_ACCESS)))
        HGOTO_ERROR(H5E_PLIST, H5E_BADTYPE, FAIL, "not a file access list")
    if(H5FD_MPIO != H5P_peek_driver(plist))
        HGOTO_ERROR(H5E_PLIST, H5E_BADVALUE, FAIL, "incorrect VFL driver")
    if(NULL == (fa = (const H5FD_mpio_fapl_t *)H5P_peek_driver_info(plist)))
        HGOTO_ERROR(H5E_PLIST, H5E_BADVALUE, FAIL, "bad VFL driver info")

    /* Store the duplicated communicator in a temporary variable for error */
    /* recovery in case the INFO duplication fails. */
    if(comm) {
	if(MPI_SUCCESS != (mpi_code = MPI_Comm_dup(fa->comm, &comm_tmp)))
	    HMPI_GOTO_ERROR(FAIL, "MPI_Comm_dup failed", mpi_code)
        comm_copied = TRUE;
    } /* end if */

    if(info) {
	if(MPI_INFO_NULL != fa->info) {
	    if(MPI_SUCCESS != (mpi_code = MPI_Info_dup(fa->info, info)))
		HMPI_GOTO_ERROR(FAIL, "MPI_Info_dup failed", mpi_code)
	} /* end if */
        else
	    /* do not dup it */
	    *info = MPI_INFO_NULL;
    } /* end if */

    /* Store the copied communicator, now that the Info object has been
     *  successfully copied.
     */
    if(comm)
        *comm = comm_tmp;

done:
    if(ret_value < 0)
	/* need to free anything created here */
	if(comm_copied)
	    MPI_Comm_free(&comm_tmp);

    FUNC_LEAVE_API(ret_value)
} /* end H5Pget_fapl_mpio() */


/*-------------------------------------------------------------------------
 * Function:	H5Pset_dxpl_mpio
 *
 * Purpose:	Set the data transfer property list DXPL_ID to use transfer
 *		mode XFER_MODE. The property list can then be used to control
 *		the I/O transfer mode during data I/O operations. The valid
 *		transfer modes are:
 *
 * 		H5FD_MPIO_INDEPENDENT:
 *			Use independent I/O access (the default).
 *
 * 		H5FD_MPIO_COLLECTIVE:
 *			Use collective I/O access.
 *
 * Return:	Success:	Non-negative
 * 		Failure:	Negative
 *
 * Programmer:	Albert Cheng
 *		April 2, 1998
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_dxpl_mpio(hid_t dxpl_id, H5FD_mpio_xfer_t xfer_mode)
{
    H5P_genplist_t *plist;      /* Property list pointer */
    herr_t ret_value = SUCCEED; /* Return value */

    FUNC_ENTER_API(FAIL)
    H5TRACE2("e", "iDt", dxpl_id, xfer_mode);

    /* Check arguments */
    if(dxpl_id == H5P_DEFAULT)
        HGOTO_ERROR(H5E_PLIST, H5E_BADVALUE, FAIL, "can't set values in default property list")
    if(NULL == (plist = H5P_object_verify(dxpl_id, H5P_DATASET_XFER)))
        HGOTO_ERROR(H5E_PLIST, H5E_BADTYPE, FAIL, "not a dxpl")
    if(H5FD_MPIO_INDEPENDENT != xfer_mode && H5FD_MPIO_COLLECTIVE != xfer_mode)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "incorrect xfer_mode")

    /* Set the transfer mode */
    if(H5P_set(plist, H5D_XFER_IO_XFER_MODE_NAME, &xfer_mode) < 0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTSET, FAIL, "unable to set value")

done:
    FUNC_LEAVE_API(ret_value)
} /* end H5Pset_dxpl_mpio() */


/*-------------------------------------------------------------------------
 * Function:	H5Pget_dxpl_mpio
 *
 * Purpose:	Queries the transfer mode current set in the data transfer
 *		property list DXPL_ID. This is not collective.
 *
 * Return:	Success:	Non-negative, with the transfer mode returned
 *				through the XFER_MODE argument if it is
 *				non-null.
 * 		Failure:	Negative
 *
 * Programmer:	Albert Cheng
 *		April 2, 1998
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_dxpl_mpio(hid_t dxpl_id, H5FD_mpio_xfer_t *xfer_mode/*out*/)
{
    H5P_genplist_t *plist;              /* Property list pointer */
    herr_t      ret_value = SUCCEED;    /* Return value */

    FUNC_ENTER_API(FAIL)
    H5TRACE2("e", "ix", dxpl_id, xfer_mode);

    /* Check arguments */
    if(NULL == (plist = H5P_object_verify(dxpl_id, H5P_DATASET_XFER)))
        HGOTO_ERROR(H5E_PLIST, H5E_BADTYPE, FAIL, "not a dxpl")

    /* Get the transfer mode */
    if(xfer_mode)
        if(H5P_get(plist, H5D_XFER_IO_XFER_MODE_NAME, xfer_mode) < 0)
            HGOTO_ERROR(H5E_PLIST, H5E_CANTSET, FAIL, "unable to get value")

done:
    FUNC_LEAVE_API(ret_value)
} /* end H5Pget_dxpl_mpio() */


/*-------------------------------------------------------------------------
 * Function:	H5Pset_dxpl_mpio_collective_opt
 *
 * Purpose:	To set a flag to choose linked chunk I/O or multi-chunk I/O
 *		without involving decision-making inside HDF5
 *
 * Note:	The library will do linked chunk I/O or multi-chunk I/O without
 *		involving communications for decision-making process.
 *		The library won't behave as it asks for only when we find
 *		that the low-level MPI-IO package doesn't support this.
 *
 * Return:	Success:	Non-negative
 * 		Failure:	Negative
 *
 * Programmer:	Kent Yang
 *		? ?, ?
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_dxpl_mpio_collective_opt(hid_t dxpl_id, H5FD_mpio_collective_opt_t opt_mode)
{
    H5P_genplist_t *plist;      /* Property list pointer */
    herr_t ret_value = SUCCEED; /* Return value */

    FUNC_ENTER_API(FAIL)
    H5TRACE2("e", "iDc", dxpl_id, opt_mode);

    /* Check arguments */
    if(dxpl_id == H5P_DEFAULT)
        HGOTO_ERROR(H5E_PLIST, H5E_BADVALUE, FAIL, "can't set values in default property list")
    if(NULL == (plist = H5P_object_verify(dxpl_id, H5P_DATASET_XFER)))
        HGOTO_ERROR(H5E_PLIST, H5E_BADTYPE, FAIL, "not a dxpl")

    /* Set the transfer mode */
    if(H5P_set(plist, H5D_XFER_MPIO_COLLECTIVE_OPT_NAME, &opt_mode) < 0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTSET, FAIL, "unable to set value")

done:
    FUNC_LEAVE_API(ret_value)
} /* end H5Pset_dxpl_mpio_collective_opt() */


/*-------------------------------------------------------------------------
 * Function:	H5Pset_dxpl_mpio_chunk_opt
 *
 * Purpose:	To set a flag to choose linked chunk I/O or multi-chunk I/O
 *		without involving decision-making inside HDF5
 *
 * Note:	The library will do linked chunk I/O or multi-chunk I/O without
 *		involving communications for decision-making process.
 *		The library won't behave as it asks for only when we find
 *		that the low-level MPI-IO package doesn't support this.
 *
 * Return:	Success:	Non-negative
 * 		Failure:	Negative
 *
 * Programmer:	Kent Yang
 *		? ?, ?
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_dxpl_mpio_chunk_opt(hid_t dxpl_id, H5FD_mpio_chunk_opt_t opt_mode)
{
    H5P_genplist_t *plist;      /* Property list pointer */
    herr_t ret_value = SUCCEED; /* Return value */

    FUNC_ENTER_API(FAIL)
    H5TRACE2("e", "iDh", dxpl_id, opt_mode);

    /* Check arguments */
    if(dxpl_id == H5P_DEFAULT)
        HGOTO_ERROR(H5E_PLIST, H5E_BADVALUE, FAIL, "can't set values in default property list")
    if(NULL == (plist = H5P_object_verify(dxpl_id, H5P_DATASET_XFER)))
        HGOTO_ERROR(H5E_PLIST, H5E_BADTYPE, FAIL, "not a dxpl")

    /* Set the transfer mode */
    if(H5P_set(plist, H5D_XFER_MPIO_CHUNK_OPT_HARD_NAME, &opt_mode) < 0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTSET, FAIL, "unable to set value")

done:
    FUNC_LEAVE_API(ret_value)
} /* end H5Pset_dxpl_mpio_chunk_opt() */


/*-------------------------------------------------------------------------
 * Function:	H5Pset_dxpl_mpio_chunk_opt_num
 *
 * Purpose:	To set a threshold for doing linked chunk IO
 *
 * Note:	If the number is greater than the threshold set by the user,
 *		the library will do linked chunk I/O; otherwise, I/O will be
 *		done for every chunk.
 *
 * Return:	Success:	Non-negative
 * 		Failure:	Negative
 *
 * Programmer:	Kent Yang
 *		? ?, ?
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_dxpl_mpio_chunk_opt_num(hid_t dxpl_id, unsigned num_chunk_per_proc)
{
    H5P_genplist_t *plist;      /* Property list pointer */
    herr_t ret_value = SUCCEED; /* Return value */

    FUNC_ENTER_API(FAIL)
    H5TRACE2("e", "iIu", dxpl_id, num_chunk_per_proc);

    /* Check arguments */
    if(dxpl_id == H5P_DEFAULT)
        HGOTO_ERROR(H5E_PLIST, H5E_BADVALUE, FAIL, "can't set values in default property list")
    if(NULL == (plist = H5P_object_verify(dxpl_id, H5P_DATASET_XFER)))
        HGOTO_ERROR(H5E_PLIST, H5E_BADTYPE, FAIL, "not a dxpl")

    /* Set the transfer mode */
    if(H5P_set(plist, H5D_XFER_MPIO_CHUNK_OPT_NUM_NAME, &num_chunk_per_proc) < 0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTSET, FAIL, "unable to set value")

done:
    FUNC_LEAVE_API(ret_value)
} /* end H5Pset_dxpl_mpio_chunk_opt_num() */


/*-------------------------------------------------------------------------
 * Function:	H5Pset_dxpl_mpio_chunk_opt_ratio
 *
 * Purpose:	To set a threshold for doing collective I/O for each chunk
 *
 * Note:	The library will calculate the percentage of the number of
 *		process holding selections at each chunk. If that percentage
 *		of number of process in the individual chunk is greater than
 *		the threshold set by the user, the library will do collective
 *		chunk I/O for this chunk; otherwise, independent I/O will be
 *		done for this chunk.
 *
 * Return:	Success:	Non-negative
 * 		Failure:	Negative
 *
 * Programmer:	Kent Yang
 *		? ?, ?
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_dxpl_mpio_chunk_opt_ratio(hid_t dxpl_id, unsigned percent_num_proc_per_chunk)
{
    H5P_genplist_t *plist;      /* Property list pointer */
    herr_t ret_value = SUCCEED; /* Return value */

    FUNC_ENTER_API(FAIL)
    H5TRACE2("e", "iIu", dxpl_id, percent_num_proc_per_chunk);

    /* Check arguments */
    if(dxpl_id == H5P_DEFAULT)
        HGOTO_ERROR(H5E_PLIST, H5E_BADVALUE, FAIL, "can't set values in default property list")
    if(NULL == (plist = H5P_object_verify(dxpl_id, H5P_DATASET_XFER)))
        HGOTO_ERROR(H5E_PLIST, H5E_BADTYPE, FAIL, "not a dxpl")

    /* Set the transfer mode */
    if(H5P_set(plist, H5D_XFER_MPIO_CHUNK_OPT_RATIO_NAME, &percent_num_proc_per_chunk) < 0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTSET, FAIL, "unable to set value")

done:
    FUNC_LEAVE_API(ret_value)
} /* end H5Pset_dxpl_mpio_chunk_opt_ratio() */


/*-------------------------------------------------------------------------
 * Function:	H5FD__mpio_fapl_get
 *
 * Purpose:	Returns a file access property list which could be used to
 *		create another file the same as this one.
 *
 * Return:	Success:	Ptr to new file access property list with all
 *				fields copied from the file pointer.
 *		Failure:	NULL
 *
 * Programmer:	Robb Matzke
 *              Friday, August 13, 1999
 *
 *-------------------------------------------------------------------------
 */
static void *
H5FD__mpio_fapl_get(H5FD_t *_file)
{
    H5FD_mpio_t		*file = (H5FD_mpio_t*)_file;
    H5FD_mpio_fapl_t	*fa = NULL;
    void      *ret_value;       /* Return value */

    FUNC_ENTER_STATIC

    /* Sanity checks */
    HDassert(file);
    HDassert(H5FD_MPIO == file->pub.driver_id);

    /* Check arguments */
    if(NULL == (fa = (H5FD_mpio_fapl_t *)H5MM_calloc(sizeof(H5FD_mpio_fapl_t))))
        HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed")

    /* Duplicate communicator and Info object. */
    if(FAIL == H5FD_mpi_comm_info_dup(file->comm, file->info, &fa->comm, &fa->info))
	HGOTO_ERROR(H5E_INTERNAL, H5E_CANTCOPY, NULL, "Communicator/Info duplicate failed")

    /* Set return value */
    ret_value = fa;

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5FD__mpio_fapl_get() */


/*-------------------------------------------------------------------------
 * Function:	H5FD__mpio_fapl_copy
 *
 * Purpose:	Copies the mpio-specific file access properties.
 *
 * Return:	Success:	Ptr to a new property list
 *		Failure:	NULL
 *
 * Programmer:	Albert Cheng
 *              Jan  8, 2003
 *
 *-------------------------------------------------------------------------
 */
static void *
H5FD__mpio_fapl_copy(const void *_old_fa)
{
    const H5FD_mpio_fapl_t *old_fa = (const H5FD_mpio_fapl_t*)_old_fa;
    H5FD_mpio_fapl_t	*new_fa = NULL;
    void		*ret_value = NULL;

    FUNC_ENTER_STATIC

#ifdef H5FDmpio_DEBUG
if(H5FD_mpio_Debug[(int)'t'])
    HDfprintf(stderr, "%s: entering\n", FUNC);
#endif

    if(NULL == (new_fa = (H5FD_mpio_fapl_t *)H5MM_malloc(sizeof(H5FD_mpio_fapl_t))))
        HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed")

    /* Copy the general information */
    HDmemcpy(new_fa, old_fa, sizeof(H5FD_mpio_fapl_t));

    /* Duplicate communicator and Info object. */
    if(H5FD_mpi_comm_info_dup(old_fa->comm, old_fa->info, &new_fa->comm, &new_fa->info) < 0)
	HGOTO_ERROR(H5E_INTERNAL, H5E_CANTCOPY, NULL, "Communicator/Info duplicate failed")

    /* Set return value */
    ret_value = new_fa;

done:
    if(NULL == ret_value)
	/* cleanup */
	if(new_fa)
	    H5MM_xfree(new_fa);

#ifdef H5FDmpio_DEBUG
if(H5FD_mpio_Debug[(int)'t'])
    HDfprintf(stderr, "%s: leaving\n", FUNC);
#endif

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5FD__mpio_fapl_copy() */


/*-------------------------------------------------------------------------
 * Function:	H5FD__mpio_fapl_free
 *
 * Purpose:	Frees the mpio-specific file access properties.
 *
 * Return:	Success:	0
 *		Failure:	-1
 *
 * Programmer:	Albert Cheng
 *              Jan  8, 2003
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD__mpio_fapl_free(void *_fa)
{
    H5FD_mpio_fapl_t	*fa = (H5FD_mpio_fapl_t*)_fa;
    herr_t		ret_value = SUCCEED;

    FUNC_ENTER_STATIC

#ifdef H5FDmpio_DEBUG
if(H5FD_mpio_Debug[(int)'t'])
    HDfprintf(stderr, "%s: Entering\n", FUNC);
#endif

    /* Sanity checks */
    HDassert(fa);

    /* Free the internal communicator and INFO object */
    HDassert(MPI_COMM_NULL != fa->comm);
    H5FD_mpi_comm_info_free(&fa->comm, &fa->info);
    H5MM_xfree(fa);

#ifdef H5FDmpio_DEBUG
if(H5FD_mpio_Debug[(int)'t'])
    HDfprintf(stderr, "%s: leaving\n", FUNC);
#endif

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5FD__mpio_fapl_free() */


/*-------------------------------------------------------------------------
 * Function:	H5FD_set_mpio_atomicity
 *
 * Purpose:	Sets the atomicity mode
 *
 * Return:	Success:	Non-negative
 *
 * 		Failure:	Negative
 *
 * Programmer:	Mohamad Chaarawi
 *		Feb 14, 2012
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5FD_set_mpio_atomicity(H5FD_t *_file, hbool_t flag)
{
    H5FD_mpio_t *file = (H5FD_mpio_t*)_file;
    int          mpi_code;               /* MPI return code */
    int          temp_flag;
    herr_t       ret_value = SUCCEED;

    FUNC_ENTER_NOAPI_NOINIT

#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'t'])
    	HDfprintf(stdout, "%s: Entering\n", FUNC);
#endif

    if(FALSE == flag)
        temp_flag = 0;
    else
        temp_flag = 1;

    /* set atomicity value */
    if(MPI_SUCCESS != (mpi_code = MPI_File_set_atomicity(file->f, temp_flag)))
        HMPI_GOTO_ERROR(FAIL, "MPI_File_set_atomicity", mpi_code)

done:
#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'t'])
    	HDfprintf(stdout, "%s: Leaving\n", FUNC);
#endif

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5FD_set_mpio_atomicity() */


/*-------------------------------------------------------------------------
 * Function:	H5FD_get_mpio_atomicity
 *
 * Purpose:	Returns the atomicity mode
 *
 * Return:	Success:	Non-negative
 *
 * 		Failure:	Negative
 *
 * Programmer:	Mohamad Chaarawi
 *		Feb 14, 2012
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5FD_get_mpio_atomicity(H5FD_t *_file, hbool_t *flag)
{
    H5FD_mpio_t *file = (H5FD_mpio_t*)_file;
    int          mpi_code;               /* MPI return code */
    int          temp_flag;
    herr_t       ret_value = SUCCEED;

    FUNC_ENTER_NOAPI_NOINIT

#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'t'])
    	HDfprintf(stdout, "%s: Entering\n", FUNC);
#endif

    /* Get atomicity value */
    if(MPI_SUCCESS != (mpi_code = MPI_File_get_atomicity(file->f, &temp_flag)))
        HMPI_GOTO_ERROR(FAIL, "MPI_File_get_atomicity", mpi_code)

    if(0 != temp_flag)
        *flag = TRUE;
    else
        *flag = FALSE;

done:
#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'t'])
    	HDfprintf(stdout, "%s: Leaving\n", FUNC);
#endif

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5FD_get_mpio_atomicity() */

/*
 * Function:    H5FD_mpio_fileinfo_get
 *
 * Purpose:     Implements a normal (posix) file open for MPI rank 0.
 *              Replicates the functionality of H5FD_sec2_open. We
 *              open the file and cache a few key structures before
 *              closing.  These cached structures are those which
 *              are eventually utilized for MPIO file comparisons.
 *
 *              N.B. The file handles returned by the collective MPI
 *              File open function are not guaranteed to have a relation to
 *              an actual posix file handle.  This then, provides the
 *              requirement that we do a "normal" file open to provide
 *              an actual file handle with which we can gather more
 *              detailed information to eventually implement file
 *              comparisons (see: H5FD__mpio_cmp)
 *
 * Return:      Success:        Non-negative
 *
 *              Failure:        Negative
 *              Indicates too, that the information used for
 *              MPIO file comparisons will most likely not
 *              be initialized and this in turn can lead to
 *              runtime issues, e.g. File comparison failures.
 */
static herr_t
H5FD_mpio_fileinfo_get(const char *name, unsigned flags, H5FD_mpio_t *file)
{
    int             status;
    int             fd          = -1;       /* File descriptor          */
    int             o_flags;                /* Flags for open() call    */
    h5_stat_t       sb;
    herr_t          ret_value = SUCCEED;

    FUNC_ENTER_NOAPI_NOINIT

#ifdef H5FDmpio_DEBUG
    if (H5FD_mpio_Debug[(int)'t'])
        fprintf(stdout, "Entering H5FD_mpio_fileinfo_get\n");
#endif

#ifdef H5_HAVE_WIN32_API
    struct _BY_HANDLE_FILE_INFORMATION fileinfo;
    HANDLE          hFile;                  /* Native windows file handle */
#endif
    o_flags = (H5F_ACC_RDWR & flags) ? O_RDWR : O_RDONLY;

    /* Open the file */
    if((fd = HDopen(name, o_flags, H5_POSIX_CREATE_MODE_RW)) < 0)
        HMPI_GOTO_ERROR(FAIL, "HDopen failed", fd)

    if((status = HDfstat(fd, &sb)) < 0)
        HMPI_GOTO_ERROR(FAIL, "HDfstat failed", status)

#ifdef H5_HAVE_WIN32_API
    hFile = (HANDLE)_get_osfhandle(fd);
    if(INVALID_HANDLE_VALUE == hFile)
        HMPI_GOTO_ERROR(FAIL, "_get_osfhandle failed", -1)

    if(!GetFileInformationByHandle((HANDLE)hFile, &fileinfo))
        HMPI_GOTO_ERROR(FAIL, "GetFileInformationByHandle failed", 0)

    file->nFileIndexHigh = fileinfo.nFileIndexHigh;
    file->nFileIndexLow = fileinfo.nFileIndexLow;
    file->dwVolumeSerialNumber = fileinfo.dwVolumeSerialNumber;
#else /* H5_HAVE_WIN32_API */
    file->device = sb.st_dev;
    file->inode = sb.st_ino;
#endif /* H5_HAVE_WIN32_API */

done:
    if(fd >= 0)
        HDclose(fd);

#ifdef H5FDmpio_DEBUG
    if (H5FD_mpio_Debug[(int)'t'])
        fprintf(stdout, "Leaving H5FD_mpio_fileinfo_get\n");
#endif
    FUNC_LEAVE_NOAPI(ret_value)
}


/*-------------------------------------------------------------------------
 * Function:    H5FD__mpio_open
 *
 * Purpose:     Opens a file with name NAME.  The FLAGS are a bit field with
 *		purpose similar to the second argument of open(2) and which
 *		are defined in H5Fpublic.h. The file access property list
 *		FAPL_ID contains the properties driver properties and MAXADDR
 *		is the largest address which this file will be expected to
 *		access.  This is collective.
 *
 * Return:      Success:        A new file pointer.
 *              Failure:        NULL
 *
 * Programmer:  Robert Kim Yates
 *              January 30, 1998
 *
 *-------------------------------------------------------------------------
 */
static H5FD_t *
H5FD__mpio_open(const char *name, unsigned flags, hid_t fapl_id,
    haddr_t H5_ATTR_UNUSED maxaddr)
{
    H5FD_mpio_t			*file = NULL;
    MPI_File			fh;
    unsigned                    file_opened = 0; /* Flag to indicate that the file was successfully opened */
    int				mpi_amode;
    int				mpi_rank;       /* MPI rank of this process */
    int				mpi_size;       /* Total number of MPI processes */
    int				mpi_code;	/* mpi return code */
    MPI_Offset			size;
    const H5FD_mpio_fapl_t	*fa = NULL;
    H5FD_mpio_fapl_t		_fa;
    H5P_genplist_t *plist;      /* Property list pointer */
    MPI_Comm                    comm_dup = MPI_COMM_NULL;
    MPI_Info                    info_dup = MPI_INFO_NULL;
    H5FD_t			*ret_value;     /* Return value */

    FUNC_ENTER_STATIC

#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'t'])
    	HDfprintf(stdout, "%s: Entering - name = \"%s\", flags = 0x%x, fapl_id = %d, maxaddr = %lu\n", FUNC, name, flags, (int)fapl_id, (unsigned long)maxaddr);
#endif

    /* Obtain a pointer to mpio-specific file access properties */
    if(NULL == (plist = H5P_object_verify(fapl_id, H5P_FILE_ACCESS)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, NULL, "not a file access property list")
    if(H5P_FILE_ACCESS_DEFAULT == fapl_id || H5FD_MPIO != H5P_peek_driver(plist)) {
	_fa.comm = MPI_COMM_SELF; /*default*/
	_fa.info = MPI_INFO_NULL; /*default*/
	fa = &_fa;
    } /* end if */
    else
        if(NULL == (fa = (const H5FD_mpio_fapl_t *)H5P_peek_driver_info(plist)))
	    HGOTO_ERROR(H5E_PLIST, H5E_BADVALUE, NULL, "bad VFL driver info")

    /* Duplicate communicator and Info object for use by this file. */
    if(FAIL == H5FD_mpi_comm_info_dup(fa->comm, fa->info, &comm_dup, &info_dup))
	HGOTO_ERROR(H5E_INTERNAL, H5E_CANTCOPY, NULL, "Communicator/Info duplicate failed")

    /* convert HDF5 flags to MPI-IO flags */
    /* some combinations are illegal; let MPI-IO figure it out */
    mpi_amode  = (flags & H5F_ACC_RDWR) ? MPI_MODE_RDWR : MPI_MODE_RDONLY;
    if(flags & H5F_ACC_CREAT)
        mpi_amode |= MPI_MODE_CREATE;
    if(flags & H5F_ACC_EXCL)
        mpi_amode |= MPI_MODE_EXCL;

#ifdef H5FDmpio_DEBUG
    /* Check for debug commands in the info parameter */
    if(MPI_INFO_NULL != info_dup) {
        char debug_str[128];
        int flag;

        MPI_Info_get(fa->info, H5F_MPIO_DEBUG_KEY, sizeof(debug_str) - 1, debug_str, &flag);
        if(flag) {
            int i;

            HDfprintf(stdout, "H5FD_mpio debug flags = '%s'\n", debug_str);
            for(i = 0; debug_str[i]/*end of string*/ && i < 128/*just in case*/; ++i)
                H5FD_mpio_Debug[(int)debug_str[i]] = 1;
        } /* end if */
    } /* end if */
#endif

    if(MPI_SUCCESS != (mpi_code = MPI_File_open(comm_dup, name, mpi_amode, info_dup, &fh)))
        HMPI_GOTO_ERROR(NULL, "MPI_File_open failed", mpi_code)
    file_opened=1;

    /* Get the MPI rank of this process and the total number of processes */
    if (MPI_SUCCESS != (mpi_code=MPI_Comm_rank (comm_dup, &mpi_rank)))
        HMPI_GOTO_ERROR(NULL, "MPI_Comm_rank failed", mpi_code)
    if (MPI_SUCCESS != (mpi_code=MPI_Comm_size (comm_dup, &mpi_size)))
        HMPI_GOTO_ERROR(NULL, "MPI_Comm_size failed", mpi_code)

    /* Build the return value and initialize it */
    if(NULL == (file = (H5FD_mpio_t *)H5MM_calloc(sizeof(H5FD_mpio_t))))
        HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed")
    file->f = fh;
    file->comm = comm_dup;
    file->info = info_dup;
    file->mpi_rank = mpi_rank;
    file->mpi_size = mpi_size;

    /* Only processor p0 will get the filesize and broadcast it. */
    if (mpi_rank == 0) {
        if (MPI_SUCCESS != (mpi_code=MPI_File_get_size(fh, &size)))
            HMPI_GOTO_ERROR(NULL, "MPI_File_get_size failed", mpi_code)
    } /* end if */

    /* Broadcast file size */
    if(MPI_SUCCESS != (mpi_code = MPI_Bcast(&size, (int)sizeof(MPI_Offset), MPI_BYTE, 0, comm_dup)))
        HMPI_GOTO_ERROR(NULL, "MPI_Bcast failed", mpi_code)

    /* Determine if the file should be truncated */
    if(size && (flags & H5F_ACC_TRUNC)) {
        if (MPI_SUCCESS != (mpi_code=MPI_File_set_size(fh, (MPI_Offset)0)))
            HMPI_GOTO_ERROR(NULL, "MPI_File_set_size failed", mpi_code)

        /* Don't let any proc return until all have truncated the file. */
        if (MPI_SUCCESS!= (mpi_code=MPI_Barrier(comm_dup)))
            HMPI_GOTO_ERROR(NULL, "MPI_Barrier failed", mpi_code)

        /* File is zero size now */
        size = 0;
    } /* end if */

    /* Set the size of the file (from library's perspective) */
    file->eof = H5FD_mpi_MPIOff_to_haddr(size);
    file->local_eof = file->eof;

    if (mpi_rank == 0) {
        /* Gather some file info for future comparisons */
        if (H5FD_mpio_fileinfo_get( name, flags, file ) < 0)
            HMPI_GOTO_ERROR(NULL, "H5FD_mpio_fileinfo_get failed", -1)
    }
    /* Set return value */
    ret_value = (H5FD_t*)file;

done:
    if(ret_value == NULL) {
        if(file_opened)
            MPI_File_close(&fh);
	if(MPI_COMM_NULL != comm_dup)
	    MPI_Comm_free(&comm_dup);
	if(MPI_INFO_NULL != info_dup)
	    MPI_Info_free(&info_dup);
	if(file)
	    H5MM_xfree(file);
    } /* end if */

#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'t'])
        HDfprintf(stdout, "%s: Leaving\n", FUNC);
#endif

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5FD__mpio_open() */


/*-------------------------------------------------------------------------
 * Function:    H5FD__mpio_cmp
 *
 * Purpose:     This version of the 'cmp' function is used to compare two
 *              files which have been created and opened using the MPI-IO
 *              driver.
 *              The peculiarity of this is that unlike POSIX io, the
 *              handle returned from an MPI_File_open operation may be
 *              an abstract value and not have any relation to an actual
 *              filestem handle.  The net result is that additional
 *              filesystem information needs to be gathered to subsequently
 *              utilize the stronger filesystem based methodology used in
 *              other HDF5 drivers, e.g. H5FD_sec2_cmp()
 *              The approach is two fold:
 *              1.  The MPI communicators used to access parallel files
 *                  will be compared.
 *              2.  MPI rank 0 is tasked with collecting the additional
 *                  POSIX or Windows NTFS information that is subsequently
 *                  used here for comparison purposes.  The result is
 *                  then broadcast to the participating MPI ranks to effect
 *                  a global result.
 *
 * Return:      An integer value similar to that returned by strcmp()

 *              NOTE: This function can't FAIL.  In those cases where
 *              where we would normally return FAILURE, e.g. when MPI
 *              returns an error, we treat these as unequal comparisons.
 *
 * Programmer:  Richard Warren
 *              Originally borrowed from H5FD_sec2_cmp (Robb Matzke) and
 *              modified as described above.
 *
 *-------------------------------------------------------------------------
 */
static int
H5FD__mpio_cmp(const H5FD_t *_f1, const H5FD_t *_f2)
{
    const H5FD_mpio_t   *f1 = (const H5FD_mpio_t *)_f1;
    const H5FD_mpio_t   *f2 = (const H5FD_mpio_t *)_f2;
    int ret_value = 0;
    int cmp_value = 0;
    int mpi_result;
    MPI_Group f1_grp;
    MPI_Group f2_grp;

    FUNC_ENTER_NOAPI_NOINIT

    if ((mpi_result = MPI_Comm_group(f1->comm, &f1_grp)) != MPI_SUCCESS)
        HMPI_GOTO_ERROR(-1, "MPI_Comm_group(comm1) failed", mpi_result)

    if ((mpi_result = MPI_Comm_group(f2->comm, &f2_grp)) != MPI_SUCCESS)
        HMPI_GOTO_ERROR(-1, "MPI_Comm_group(comm2) failed", mpi_result)

    if ((mpi_result = MPI_Group_compare(f1_grp, f2_grp, &cmp_value)) != MPI_SUCCESS)
        HMPI_GOTO_ERROR(-1, "MPI_Group_compare failed", mpi_result)

    /*  The group compare return values can be one of the following:
     *   MPI_IDENT(0)     == two groups/communicators are identical
     *  ----------------  Those below can lead to unexpected
     *  ----------------  results, so we will return unequal
     *                    for the file comparison.
     *   MPI_CONGRUENT(1) == two groups/communicators are equal but
     *                       are distinct communication domains
     *   MPI_SIMILAR(2)   == two groups have the same members but
     *                       ordering may be different
     *   MPI_UNEQUAL(3)   == self descriptive (unequal)
     *
     *  Note: Congruent groups would seem to satisfy the equality
     *  condition from the file perspective, but there may be conditions
     *  in which collective operations would cause an application to
     *  hang if two different communicators are in use, e.g. any
     *  sort of synchronization (Barrier, Bcast).
     */

    if (cmp_value >= MPI_CONGRUENT)
        HGOTO_DONE(-1)

    if (f1->mpi_rank == 0) {
        /* Because MPI file handles may NOT have any relation to
         * to actual file handle, we utilize a "regular" file open
         * on MPI rank 0 prior to opening with the MPI-IO routines.
         * The H5FD_mpio_t structure is utilized to cache the
         * relevant comparison values which we use for comparisons
         * below.
         */
#ifdef H5_HAVE_WIN32_API
        if ((f1->dwVolumeSerialNumber < f2->dwVolumeSerialNumber) ||
            (f1->nFileIndexHigh       < f2->nFileIndexHigh)       ||
            (f1->nFileIndexLow        < f2->nFileIndexLow))
            cmp_value = -1;
        else
        if ((f1->dwVolumeSerialNumber > f2->dwVolumeSerialNumber) ||
	    (f1->nFileIndexHigh       > f2->nFileIndexHigh)       ||
            (f1->nFileIndexLow        > f2->nFileIndexLow))
            cmp_value = 1;
#else  /* Not WIN32 */
#ifdef H5_DEV_T_IS_SCALAR
        if (f1->device < f2->device)     cmp_value = -1;
        else if(f1->device > f2->device) cmp_value = 1;
#else /* H5_DEV_T_IS_SCALAR */
        /* If dev_t isn't a scalar value on this system, just use memcmp to
         * determine if the values are the same or not.  The actual return value
         * shouldn't really matter...
         */
        cmp_value = HDmemcmp(&(f1->device), &(f2->device), sizeof(dev_t));
#endif /* H5_DEV_T_IS_SCALAR */
        /* Continue checking ONLY if we haven't failed yet */
        if (!cmp_value) {
           if(f1->inode < f2->inode)     cmp_value = -1;
           else if(f1->inode > f2->inode) cmp_value = 1;
#endif /* H5_HAVE_WIN32_API */
        }
    }
    if (MPI_SUCCESS != (mpi_result = MPI_Bcast(&cmp_value, 1, MPI_INT, 0, f1->comm)))
        HMPI_GOTO_ERROR(-1, "MPI_Bcast failed", mpi_result)

    ret_value = cmp_value;
done:

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5FD__mpio_cmp() */


/*-------------------------------------------------------------------------
 * Function:    H5FD__mpio_close
 *
 * Purpose:     Closes a file.  This is collective.
 *
 * Return:      Success:	Non-negative
 * 		Failure:	Negative
 *
 * Programmer:  Unknown
 *              January 30, 1998
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD__mpio_close(H5FD_t *_file)
{
    H5FD_mpio_t	*file = (H5FD_mpio_t*)_file;
    int		mpi_code;	        /* MPI return code */
    herr_t      ret_value = SUCCEED;    /* Return value */

    FUNC_ENTER_STATIC

#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'t'])
    	HDfprintf(stdout, "%s: Entering\n", FUNC);
#endif

    /* Sanity checks */
    HDassert(file);
    HDassert(H5FD_MPIO == file->pub.driver_id);

    /* MPI_File_close sets argument to MPI_FILE_NULL */
    if(MPI_SUCCESS != (mpi_code = MPI_File_close(&(file->f)/*in,out*/)))
        HMPI_GOTO_ERROR(FAIL, "MPI_File_close failed", mpi_code)

    /* Clean up other stuff */
    H5FD_mpi_comm_info_free(&file->comm, &file->info);
    H5MM_xfree(file);

done:
#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'t'])
    	HDfprintf(stdout, "%s: Leaving\n", FUNC);
#endif

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5FD__mpio_close() */


/*-------------------------------------------------------------------------
 * Function:	H5FD__mpio_query
 *
 * Purpose:	Set the flags that this VFL driver is capable of supporting.
 *              (listed in H5FDpublic.h)
 *
 * Return:	Success:	non-negative
 *		Failure:	negative
 *
 * Programmer:	Quincey Koziol
 *              Friday, August 25, 2000
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD__mpio_query(const H5FD_t H5_ATTR_UNUSED *_file, unsigned long *flags /* out */)
{
    FUNC_ENTER_STATIC

    /* Set the VFL feature flags that this driver supports */
    if(flags) {
        *flags = 0;
        *flags |= H5FD_FEAT_AGGREGATE_METADATA;     /* OK to aggregate metadata allocations                             */
        *flags |= H5FD_FEAT_AGGREGATE_SMALLDATA;    /* OK to aggregate "small" raw data allocations                     */
        *flags |= H5FD_FEAT_HAS_MPI;                /* This driver uses MPI                                             */
        *flags |= H5FD_FEAT_ALLOCATE_EARLY;         /* Allocate space early instead of late                             */
        *flags |= H5FD_FEAT_DEFAULT_VFD_COMPATIBLE; /* VFD creates a file which can be opened with the default VFD      */
    } /* end if */

    FUNC_LEAVE_NOAPI(SUCCEED)
} /* end H5FD__mpio_query() */


/*-------------------------------------------------------------------------
 * Function:	H5FD__mpio_get_eoa
 *
 * Purpose:	Gets the end-of-address marker for the file. The EOA marker
 *		is the first address past the last byte allocated in the
 *		format address space.
 *
 * Return:	Success:	The end-of-address marker.
 *		Failure:	HADDR_UNDEF
 *
 * Programmer:	Robb Matzke
 *              Friday, August  6, 1999
 *
 *-------------------------------------------------------------------------
 */
static haddr_t
H5FD__mpio_get_eoa(const H5FD_t *_file, H5FD_mem_t H5_ATTR_UNUSED type)
{
    const H5FD_mpio_t	*file = (const H5FD_mpio_t*)_file;

    FUNC_ENTER_STATIC

    /* Sanity checks */
    HDassert(file);
    HDassert(H5FD_MPIO == file->pub.driver_id);

    FUNC_LEAVE_NOAPI(file->eoa)
} /* end H5FD__mpio_get_eoa() */


/*-------------------------------------------------------------------------
 * Function:	H5FD__mpio_set_eoa
 *
 * Purpose:	Set the end-of-address marker for the file. This function is
 *		called shortly after an existing HDF5 file is opened in order
 *		to tell the driver where the end of the HDF5 data is located.
 *
 * Return:	Success:	0
 *		Failure:	-1
 *
 * Programmer:	Robb Matzke
 *              Friday, August 6, 1999
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD__mpio_set_eoa(H5FD_t *_file, H5FD_mem_t H5_ATTR_UNUSED type, haddr_t addr)
{
    H5FD_mpio_t	*file = (H5FD_mpio_t*)_file;

    FUNC_ENTER_STATIC

    /* Sanity checks */
    HDassert(file);
    HDassert(H5FD_MPIO == file->pub.driver_id);

    file->eoa = addr;

    FUNC_LEAVE_NOAPI(SUCCEED)
} /* end H5FD__mpio_set_eoa() */


/*-------------------------------------------------------------------------
 * Function:	H5FD__mpio_get_eof
 *
 * Purpose:	Gets the end-of-file marker for the file. The EOF marker
 *		is the real size of the file.
 *
 *		The MPIO driver doesn't bother keeping this field updated
 *		since that's a relatively expensive operation. Fortunately
 *		the library only needs the EOF just after the file is opened
 *		in order to determine whether the file is empty, truncated,
 *		or okay.  Therefore, any MPIO I/O function will set its value
 *		to HADDR_UNDEF which is the error return value of this
 *		function.
 *
 *              Keeping the EOF updated (during write calls) is expensive
 *              because any process may extend the physical end of the
 *              file. -QAK
 *
 * Return:	Success:	The end-of-address marker.
 *		Failure:	HADDR_UNDEF
 *
 * Programmer:	Robb Matzke
 *              Friday, August  6, 1999
 *
 *-------------------------------------------------------------------------
 */
static haddr_t
H5FD__mpio_get_eof(const H5FD_t *_file, H5FD_mem_t H5_ATTR_UNUSED type)
{
    const H5FD_mpio_t	*file = (const H5FD_mpio_t*)_file;

    FUNC_ENTER_STATIC

    /* Sanity checks */
    HDassert(file);
    HDassert(H5FD_MPIO == file->pub.driver_id);

    FUNC_LEAVE_NOAPI(file->eof)
} /* end H5FD__mpio_get_eof() */


/*-------------------------------------------------------------------------
 * Function:       H5FD__mpio_get_handle
 *
 * Purpose:        Returns the file handle of MPIO file driver.
 *
 * Returns:        Non-negative if succeed or negative if fails.
 *
 * Programmer:     Raymond Lu
 *                 Sept. 16, 2002
 *
 *-------------------------------------------------------------------------
*/
static herr_t
H5FD__mpio_get_handle(H5FD_t *_file, hid_t H5_ATTR_UNUSED fapl, void** file_handle)
{
    H5FD_mpio_t *file = (H5FD_mpio_t *)_file;
    herr_t      ret_value = SUCCEED;

    FUNC_ENTER_STATIC

    if(!file_handle)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "file handle not valid")

    *file_handle = &(file->f);

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5FD__mpio_get_handle() */


/*-------------------------------------------------------------------------
 * Function:       H5FD_mpio_get_info
 *
 * Purpose:        Returns the file info of MPIO file driver.
 *
 * Returns:        Non-negative if succeed or negative if fails.
 *
 * Programmer:     John Mainzer
 *                 April 4, 2017
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
*/
static herr_t
H5FD_mpio_get_info(H5FD_t *_file, void** mpi_info)
{
    H5FD_mpio_t         *file = (H5FD_mpio_t *)_file;
    herr_t              ret_value = SUCCEED;

    FUNC_ENTER_NOAPI_NOINIT

    if(!mpi_info)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "mpi info not valid")

    *mpi_info = &(file->info);

done:
    FUNC_LEAVE_NOAPI(ret_value)

} /* H5FD_mpio_get_info() */


/*-------------------------------------------------------------------------
 * Function:	H5FD__mpio_read
 *
 * Purpose:	Reads SIZE bytes of data from FILE beginning at address ADDR
 *		into buffer BUF according to data transfer properties in
 *		DXPL_ID using potentially complex file and buffer types to
 *		effect the transfer.
 *
 *		Reading past the end of the MPI file returns zeros instead of
 *		failing.  MPI is able to coalesce requests from different
 *		processes (collective or independent).
 *
 * Return:	Success:	Zero. Result is stored in caller-supplied
 *				buffer BUF.
 *
 *		Failure:	-1, Contents of buffer BUF are undefined.
 *
 * Programmer:	rky, 1998-01-30
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD__mpio_read(H5FD_t *_file, H5FD_mem_t H5_ATTR_UNUSED type,
    hid_t H5_ATTR_UNUSED dxpl_id, haddr_t addr, size_t size, void *buf/*out*/)
{
    H5FD_mpio_t *file = (H5FD_mpio_t*)_file;
    MPI_Offset  mpi_off;
    MPI_Status  mpi_stat;       /* Status from I/O operation */
    int         mpi_code;       /* mpi return code */
    MPI_Datatype buf_type = MPI_BYTE;      /* MPI description of the selection in memory */
    int         size_i;         /* Integer copy of 'size' to read */
#if MPI_VERSION >= 3
    MPI_Count   bytes_read = 0; /* Number of bytes read in */
    MPI_Count   type_size;      /* MPI datatype used for I/O's size */
    MPI_Count   io_size;        /* Actual number of bytes requested */
    MPI_Count   n;
#else
    int         bytes_read = 0; /* Number of bytes read in */
    int         type_size;      /* MPI datatype used for I/O's size */
    int         io_size;        /* Actual number of bytes requested */
    int         n;
#endif
    hbool_t     use_view_this_time = FALSE;
    hbool_t     rank0_bcast = FALSE; /* If read-with-rank0-and-bcast flag was used */
    herr_t      ret_value = SUCCEED;

    FUNC_ENTER_STATIC

#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'t'])
        HDfprintf(stdout, "%s: Entering\n", FUNC);
#endif

    /* Sanity checks */
    HDassert(file);
    HDassert(H5FD_MPIO==file->pub.driver_id);
    HDassert(buf);

    /* Portably initialize MPI status variable */
    HDmemset(&mpi_stat,0,sizeof(MPI_Status));

    /* some numeric conversions */
    if(H5FD_mpi_haddr_to_MPIOff(addr, &mpi_off/*out*/) < 0)
        HGOTO_ERROR(H5E_INTERNAL, H5E_BADRANGE, FAIL, "can't convert from haddr to MPI off")
    size_i = (int)size;
    if((hsize_t)size_i != size)
        HGOTO_ERROR(H5E_INTERNAL, H5E_BADRANGE, FAIL, "can't convert from size to size_i")

#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'r'])
        HDfprintf(stdout, "%s: mpi_off = %ld  size_i = %d\n", FUNC, (long)mpi_off, size_i);
#endif

    /* Only look for MPI views for raw data transfers */
    if(type == H5FD_MEM_DRAW) {
        H5FD_mpio_xfer_t xfer_mode;   /* I/O transfer mode */

        /* Get the transfer mode from the API context */
        if(H5CX_get_io_xfer_mode(&xfer_mode) < 0)
            HGOTO_ERROR(H5E_VFL, H5E_CANTGET, FAIL, "can't get MPI-I/O transfer mode")

        /*
         * Set up for a fancy xfer using complex types, or single byte block. We
         * wouldn't need to rely on the use_view field if MPI semantics allowed
         * us to test that btype=ftype=MPI_BYTE (or even MPI_TYPE_NULL, which
         * could mean "use MPI_BYTE" by convention).
         */
        if(xfer_mode == H5FD_MPIO_COLLECTIVE) {
            MPI_Datatype file_type;

            /* Remember that views are used */
            use_view_this_time = TRUE;

            /* Prepare for a full-blown xfer using btype, ftype, and disp */
            if(H5CX_get_mpi_coll_datatypes(&buf_type, &file_type) < 0)
                HGOTO_ERROR(H5E_VFL, H5E_CANTGET, FAIL, "can't get MPI-I/O datatypes")

            /*
             * Set the file view when we are using MPI derived types
             */
            if(MPI_SUCCESS != (mpi_code = MPI_File_set_view(file->f, mpi_off, MPI_BYTE, file_type, H5FD_mpi_native_g, file->info)))
                HMPI_GOTO_ERROR(FAIL, "MPI_File_set_view failed", mpi_code)

            /* When using types, use the address as the displacement for
             * MPI_File_set_view and reset the address for the read to zero
             */
            mpi_off = 0;
        } /* end if */
    } /* end if */

    /* Read the data. */
    if(use_view_this_time) {
        H5FD_mpio_collective_opt_t coll_opt_mode;

#ifdef H5FDmpio_DEBUG
        if(H5FD_mpio_Debug[(int)'r'])
            HDfprintf(stdout, "%s: using MPIO collective mode\n", FUNC);
#endif
        /* Get the collective_opt property to check whether the application wants to do IO individually. */
        if(H5CX_get_mpio_coll_opt(&coll_opt_mode) < 0)
            HGOTO_ERROR(H5E_VFL, H5E_CANTGET, FAIL, "can't get MPI-I/O collective_op property")

        if(coll_opt_mode == H5FD_MPIO_COLLECTIVE_IO) {
#ifdef H5FDmpio_DEBUG
            if(H5FD_mpio_Debug[(int)'r'])
                HDfprintf(stdout, "%s: doing MPI collective IO\n", FUNC);
#endif
            /* Check whether we should read from rank 0 and broadcast to other ranks */
            if(H5CX_get_mpio_rank0_bcast()) {
#ifdef H5FDmpio_DEBUG
                if(H5FD_mpio_Debug[(int)'r'])
                    HDfprintf(stdout, "%s: doing read-rank0-and-MPI_Bcast\n", FUNC);
#endif
                /* Indicate path we've taken */
                rank0_bcast = TRUE;

                /* Read on rank 0 Bcast to other ranks */
                if(file->mpi_rank == 0)
                    if(MPI_SUCCESS != (mpi_code = MPI_File_read_at(file->f, mpi_off, buf, size_i, buf_type, &mpi_stat)))
                        HMPI_GOTO_ERROR(FAIL, "MPI_File_read_at failed", mpi_code)
                if(MPI_SUCCESS != (mpi_code = MPI_Bcast(buf, size_i, buf_type, 0, file->comm)))
                    HMPI_GOTO_ERROR(FAIL, "MPI_Bcast failed", mpi_code)
            } /* end if */
            else
                if(MPI_SUCCESS != (mpi_code = MPI_File_read_at_all(file->f, mpi_off, buf, size_i, buf_type, &mpi_stat)))
                    HMPI_GOTO_ERROR(FAIL, "MPI_File_read_at_all failed", mpi_code)
        } /* end if */
        else {
#ifdef H5FDmpio_DEBUG
            if(H5FD_mpio_Debug[(int)'r'])
                HDfprintf(stdout, "%s: doing MPI independent IO\n", FUNC);
#endif

            if(MPI_SUCCESS != (mpi_code = MPI_File_read_at(file->f, mpi_off, buf, size_i, buf_type, &mpi_stat)))
                HMPI_GOTO_ERROR(FAIL, "MPI_File_read_at failed", mpi_code)
        } /* end else */

        /*
         * Reset the file view when we used MPI derived types
         */
        if(MPI_SUCCESS != (mpi_code = MPI_File_set_view(file->f, (MPI_Offset)0, MPI_BYTE, MPI_BYTE, H5FD_mpi_native_g, file->info)))
            HMPI_GOTO_ERROR(FAIL, "MPI_File_set_view failed", mpi_code)
    } /* end if */
    else
        if(MPI_SUCCESS != (mpi_code = MPI_File_read_at(file->f, mpi_off, buf, size_i, buf_type, &mpi_stat)))
            HMPI_GOTO_ERROR(FAIL, "MPI_File_read_at failed", mpi_code)

    /* Only retrieve bytes read if this rank _actually_ participated in I/O */
    if(!rank0_bcast || (rank0_bcast && file->mpi_rank == 0) ) {
        /* How many bytes were actually read? */
#if MPI_VERSION >= 3
        if(MPI_SUCCESS != (mpi_code = MPI_Get_elements_x(&mpi_stat, buf_type, &bytes_read)))
#else
        if(MPI_SUCCESS != (mpi_code = MPI_Get_elements(&mpi_stat, MPI_BYTE, &bytes_read)))
#endif
            HMPI_GOTO_ERROR(FAIL, "MPI_Get_elements failed", mpi_code)
    } /* end if */
 
    /* If the rank0-bcast feature was used, broadcast the # of bytes read to
     * other ranks, which didn't perform any I/O.
     */
    /* NOTE: This could be optimized further to be combined with the broadcast
     *          of the data.  (QAK - 2019/1/2)
     */
    if(rank0_bcast)
        if(MPI_SUCCESS != MPI_Bcast(&bytes_read, 1, MPI_LONG_LONG, 0, file->comm))
            HMPI_GOTO_ERROR(FAIL, "MPI_Bcast failed", 0)

    /* Get the type's size */
#if MPI_VERSION >= 3
    if(MPI_SUCCESS != (mpi_code = MPI_Type_size_x(buf_type, &type_size)))
#else
    if(MPI_SUCCESS != (mpi_code = MPI_Type_size(buf_type, &type_size)))
#endif
        HMPI_GOTO_ERROR(FAIL, "MPI_Type_size failed", mpi_code)

    /* Compute the actual number of bytes requested */
    io_size = type_size * size_i;

    /* Check for read failure */
    if(bytes_read < 0 || bytes_read > io_size)
        HGOTO_ERROR(H5E_IO, H5E_READERROR, FAIL, "file read failed")

    /*
     * This gives us zeroes beyond end of physical MPI file.
     */
    if((n = (io_size - bytes_read)) > 0)
        HDmemset((char*)buf+bytes_read, 0, (size_t)n);

done:
#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'t'])
       HDfprintf(stdout, "%s: Leaving\n", FUNC);
#endif

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5FD__mpio_read() */


/*-------------------------------------------------------------------------
 * Function:	H5FD__mpio_write
 *
 * Purpose:	Writes SIZE bytes of data to FILE beginning at address ADDR
 *		from buffer BUF according to data transfer properties in
 *		DXPL_ID using potentially complex file and buffer types to
 *		effect the transfer.
 *
 *		MPI is able to coalesce requests from different processes
 *		(collective and independent).
 *
 * Return:	Success:	Zero. USE_TYPES and OLD_USE_TYPES in the
 *				access params are altered.
 *		Failure:	-1, USE_TYPES and OLD_USE_TYPES in the
 *				access params may be altered.
 *
 * Programmer:	Robert Kim Yates
 *              January 30, 1998
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD__mpio_write(H5FD_t *_file, H5FD_mem_t type, hid_t H5_ATTR_UNUSED dxpl_id,
    haddr_t addr, size_t size, const void *buf)
{
    H5FD_mpio_t			*file = (H5FD_mpio_t*)_file;
    MPI_Offset 		 	mpi_off;
    MPI_Status  		mpi_stat;       /* Status from I/O operation */
    MPI_Datatype		buf_type = MPI_BYTE;      /* MPI description of the selection in memory */
    int			        mpi_code;	/* MPI return code */
#if MPI_VERSION >= 3
    MPI_Count         		bytes_written;
    MPI_Count                   type_size;      /* MPI datatype used for I/O's size */
    MPI_Count                   io_size;        /* Actual number of bytes requested */
#else
    int                         bytes_written;
    int                         type_size;      /* MPI datatype used for I/O's size */
    int                         io_size;        /* Actual number of bytes requested */
#endif
    int                         size_i;
    hbool_t			use_view_this_time = FALSE;
    H5FD_mpio_xfer_t            xfer_mode;   /* I/O transfer mode */
    herr_t              	ret_value = SUCCEED;

    FUNC_ENTER_STATIC

#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'t'])
    	HDfprintf(stdout, "%s: Entering\n", FUNC);
#endif

    /* Sanity checks */
    HDassert(file);
    HDassert(H5FD_MPIO==file->pub.driver_id);
    HDassert(buf);

    /* Verify that no data is written when between MPI_Barrier()s during file flush */
    HDassert(!H5CX_get_mpi_file_flushing());

    /* Portably initialize MPI status variable */
    HDmemset(&mpi_stat, 0, sizeof(MPI_Status));

    /* some numeric conversions */
    if(H5FD_mpi_haddr_to_MPIOff(addr, &mpi_off) < 0)
        HGOTO_ERROR(H5E_INTERNAL, H5E_BADRANGE, FAIL, "can't convert from haddr to MPI off")
    size_i = (int)size;
    if((hsize_t)size_i != size)
        HGOTO_ERROR(H5E_INTERNAL, H5E_BADRANGE, FAIL, "can't convert from size to size_i")

#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'w'])
        HDfprintf(stdout, "%s: mpi_off = %ld  size_i = %d\n", FUNC, (long)mpi_off, size_i);
#endif

    /* Get the transfer mode from the API context */
    if(H5CX_get_io_xfer_mode(&xfer_mode) < 0)
        HGOTO_ERROR(H5E_VFL, H5E_CANTGET, FAIL, "can't get MPI-I/O transfer mode")

    /*
     * Set up for a fancy xfer using complex types, or single byte block. We
     * wouldn't need to rely on the use_view field if MPI semantics allowed
     * us to test that btype=ftype=MPI_BYTE (or even MPI_TYPE_NULL, which
     * could mean "use MPI_BYTE" by convention).
     */
    if(xfer_mode == H5FD_MPIO_COLLECTIVE) {
        MPI_Datatype		file_type;

        /* Remember that views are used */
        use_view_this_time = TRUE;

        /* Prepare for a full-blown xfer using btype, ftype, and disp */
        if(H5CX_get_mpi_coll_datatypes(&buf_type, &file_type) < 0)
            HGOTO_ERROR(H5E_VFL, H5E_CANTGET, FAIL, "can't get MPI-I/O datatypes")

        /*
         * Set the file view when we are using MPI derived types
         */
        if(MPI_SUCCESS != (mpi_code = MPI_File_set_view(file->f, mpi_off, MPI_BYTE, file_type, H5FD_mpi_native_g, file->info)))
            HMPI_GOTO_ERROR(FAIL, "MPI_File_set_view failed", mpi_code)

        /* When using types, use the address as the displacement for
         * MPI_File_set_view and reset the address for the read to zero
         */
        mpi_off = 0;
    } /* end if */

    /* Write the data. */
    if(use_view_this_time) {
        H5FD_mpio_collective_opt_t coll_opt_mode;

#ifdef H5FDmpio_DEBUG
        if(H5FD_mpio_Debug[(int)'w'])
            HDfprintf(stdout, "%s: using MPIO collective mode\n", FUNC);
#endif

        /* Get the collective_opt property to check whether the application wants to do IO individually. */
        if(H5CX_get_mpio_coll_opt(&coll_opt_mode) < 0)
            HGOTO_ERROR(H5E_VFL, H5E_CANTGET, FAIL, "can't get MPI-I/O collective_op property")

        if(coll_opt_mode == H5FD_MPIO_COLLECTIVE_IO) {
#ifdef H5FDmpio_DEBUG
            if(H5FD_mpio_Debug[(int)'w'])
                HDfprintf(stdout, "%s: doing MPI collective IO\n", FUNC);
#endif
            if(MPI_SUCCESS != (mpi_code = MPI_File_write_at_all(file->f, mpi_off, buf, size_i, buf_type, &mpi_stat)))
                HMPI_GOTO_ERROR(FAIL, "MPI_File_write_at_all failed", mpi_code)
        } /* end if */
        else {
            if(type != H5FD_MEM_DRAW)
                HGOTO_ERROR(H5E_PLIST, H5E_BADTYPE, FAIL, "Metadata Coll opt property should be collective at this point")
#ifdef H5FDmpio_DEBUG
            if(H5FD_mpio_Debug[(int)'w'])
                HDfprintf(stdout, "%s: doing MPI independent IO\n", FUNC);
#endif
            if(MPI_SUCCESS != (mpi_code = MPI_File_write_at(file->f, mpi_off, buf, size_i, buf_type, &mpi_stat)))
                HMPI_GOTO_ERROR(FAIL, "MPI_File_write_at failed", mpi_code)
        } /* end else */

        /* Reset the file view when we used MPI derived types */
        if(MPI_SUCCESS != (mpi_code = MPI_File_set_view(file->f, (MPI_Offset)0, MPI_BYTE, MPI_BYTE, H5FD_mpi_native_g,  file->info)))
            HMPI_GOTO_ERROR(FAIL, "MPI_File_set_view failed", mpi_code)
    } /* end if */
    else
        if(MPI_SUCCESS != (mpi_code = MPI_File_write_at(file->f, mpi_off, buf, size_i, buf_type, &mpi_stat)))
            HMPI_GOTO_ERROR(FAIL, "MPI_File_write_at failed", mpi_code)

    /* How many bytes were actually written? */
#if MPI_VERSION >= 3
    if(MPI_SUCCESS != (mpi_code = MPI_Get_elements_x(&mpi_stat, buf_type, &bytes_written)))
#else
    if(MPI_SUCCESS != (mpi_code = MPI_Get_elements(&mpi_stat, MPI_BYTE, &bytes_written)))
#endif
        HMPI_GOTO_ERROR(FAIL, "MPI_Get_elements failed", mpi_code)

    /* Get the type's size */
#if MPI_VERSION >= 3
    if(MPI_SUCCESS != (mpi_code = MPI_Type_size_x(buf_type, &type_size)))
#else
    if(MPI_SUCCESS != (mpi_code = MPI_Type_size(buf_type, &type_size)))
#endif
        HMPI_GOTO_ERROR(FAIL, "MPI_Type_size failed", mpi_code)

    /* Compute the actual number of bytes requested */
    io_size = type_size * size_i;

    /* Check for write failure */
    if(bytes_written != io_size)
        HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "file write failed")

    /* Each process will keep track of its perceived EOF value locally, and
     * ultimately we will reduce this value to the maximum amongst all
     * processes, but until then keep the actual eof at HADDR_UNDEF just in
     * case something bad happens before that point. (rather have a value
     * we know is wrong sitting around rather than one that could only
     * potentially be wrong.) */
    file->eof = HADDR_UNDEF;

    if(bytes_written && ((bytes_written + addr) > file->local_eof))
        file->local_eof = addr + bytes_written;

done:
#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'t'])
    	HDfprintf(stdout, "%s: Leaving, proc %d: ret_value = %d\n", FUNC, file->mpi_rank, ret_value );
#endif

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5FD__mpio_write() */


/*-------------------------------------------------------------------------
 * Function:    H5FD__mpio_flush
 *
 * Purpose:     Makes sure that all data is on disk.  This is collective.
 *
 * Return:      Success:	Non-negative
 *
 * 		Failure:	Negative
 *
 * Programmer:  Robb Matzke
 *              January 30, 1998
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD__mpio_flush(H5FD_t *_file, hid_t H5_ATTR_UNUSED dxpl_id, hbool_t closing)
{
    H5FD_mpio_t		*file = (H5FD_mpio_t*)_file;
    int			mpi_code;	/* mpi return code */
    herr_t              ret_value = SUCCEED;

    FUNC_ENTER_STATIC

#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'t'])
    	HDfprintf(stdout, "%s: Entering\n", FUNC);
#endif

    /* Sanity checks */
    HDassert(file);
    HDassert(H5FD_MPIO == file->pub.driver_id);

    /* Only sync the file if we are not going to immediately close it */
    if(!closing)
        if(MPI_SUCCESS != (mpi_code = MPI_File_sync(file->f)))
            HMPI_GOTO_ERROR(FAIL, "MPI_File_sync failed", mpi_code)

done:
#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'t'])
    	HDfprintf(stdout, "%s: Leaving\n", FUNC);
#endif

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5FD__mpio_flush() */


/*-------------------------------------------------------------------------
 * Function:    H5FD__mpio_truncate
 *
 * Purpose:     Make certain the file's size matches it's allocated size
 *
 *              This is a little sticky in the mpio case, as it is not 
 *              easy for us to track the current EOF by extracting it from
 *              write calls. 
 *
 *              Instead, we first check to see if the eoa has changed since
 *              the last call to this function.  If it has, we call 
 *              MPI_File_get_size() to determine the current EOF, and 
 *              only call MPI_File_set_size() if this value disagrees 
 *              with the current eoa.
 *
 * Return:      Success:	Non-negative
 * 		Failure:	Negative
 *
 * Programmer:  Quincey Koziol
 *              January 31, 2008
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD__mpio_truncate(H5FD_t *_file, hid_t H5_ATTR_UNUSED dxpl_id, hbool_t H5_ATTR_UNUSED closing)
{
    H5FD_mpio_t	*file = (H5FD_mpio_t*)_file;
    herr_t      ret_value = SUCCEED;

    FUNC_ENTER_STATIC

#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'t'])
    	HDfprintf(stdout, "%s: Entering\n", FUNC);
#endif

    /* Sanity checks */
    HDassert(file);
    HDassert(H5FD_MPIO == file->pub.driver_id);

    if(!H5F_addr_eq(file->eoa, file->last_eoa)) {
        int             mpi_code;       /* mpi return code */
        MPI_Offset      size;
        MPI_Offset      needed_eof;

        /* In principle, it is possible for the size returned by the 
         * call to MPI_File_get_size() to depend on whether writes from 
         * all proceeses have completed at the time process 0 makes the
         * call.  
         *
         * In practice, most (all?) truncate calls will come after a barrier
         * and with no interviening writes to the file (with the possible 
         * exception of sueprblock / superblock extension message updates).
         *
         * Check the "MPI file closing" flag in the API context to determine
         * if we can skip the barrier.
         */
        if(!H5CX_get_mpi_file_flushing())
            if(MPI_SUCCESS != (mpi_code = MPI_Barrier(file->comm)))
                HMPI_GOTO_ERROR(FAIL, "MPI_Barrier failed", mpi_code)

        /* Only processor p0 will get the filesize and broadcast it. */
        /* (Note that throwing an error here will cause non-rank 0 processes
         *      to hang in following Bcast.  -QAK, 3/17/2018)
         */
        if(0 == file->mpi_rank)
            if(MPI_SUCCESS != (mpi_code = MPI_File_get_size(file->f, &size)))
                HMPI_GOTO_ERROR(FAIL, "MPI_File_get_size failed", mpi_code)

        /* Broadcast file size */
        if(MPI_SUCCESS != (mpi_code = MPI_Bcast(&size, (int)sizeof(MPI_Offset), MPI_BYTE, 0, file->comm)))
            HMPI_GOTO_ERROR(FAIL, "MPI_Bcast failed", mpi_code)

        if(H5FD_mpi_haddr_to_MPIOff(file->eoa, &needed_eof) < 0)
            HGOTO_ERROR(H5E_INTERNAL, H5E_BADRANGE, FAIL, "cannot convert from haddr_t to MPI_Offset")

        /* eoa != eof.  Set eof to eoa */
        if(size != needed_eof) {
            /* Extend the file's size */
            if(MPI_SUCCESS != (mpi_code = MPI_File_set_size(file->f, needed_eof)))
                HMPI_GOTO_ERROR(FAIL, "MPI_File_set_size failed", mpi_code)

            /* In general, we must wait until all processes have finished 
             * the truncate before any process can continue, since it is 
             * possible that a process would write at the end of the 
             * file, and this write would be discarded by the truncate.
             *
             * While this is an issue for a user initiated flush, it may 
             * not be an issue at file close.  If so, we may be able to 
             * optimize out the following barrier in that case.
             */
            if(MPI_SUCCESS != (mpi_code = MPI_Barrier(file->comm)))
                HMPI_GOTO_ERROR(FAIL, "MPI_Barrier failed", mpi_code)
        } /* end if */

        /* Update the 'last' eoa value */
        file->last_eoa = file->eoa;
    } /* end if */

done:
#ifdef H5FDmpio_DEBUG
    if(H5FD_mpio_Debug[(int)'t'])
    	HDfprintf(stdout, "%s: Leaving\n", FUNC);
#endif

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5FD__mpio_truncate() */


/*-------------------------------------------------------------------------
 * Function:	H5FD__mpio_mpi_rank
 *
 * Purpose:	Returns the MPI rank for a process
 *
 * Return:	Success: non-negative
 *		Failure: negative
 *
 * Programmer:	Quincey Koziol
 *              Thursday, May 16, 2002
 *
 *-------------------------------------------------------------------------
 */
static int
H5FD__mpio_mpi_rank(const H5FD_t *_file)
{
    const H5FD_mpio_t	*file = (const H5FD_mpio_t*)_file;

    FUNC_ENTER_STATIC

    /* Sanity checks */
    HDassert(file);
    HDassert(H5FD_MPIO == file->pub.driver_id);

    FUNC_LEAVE_NOAPI(file->mpi_rank)
} /* end H5FD__mpio_mpi_rank() */


/*-------------------------------------------------------------------------
 * Function:	H5FD__mpio_mpi_size
 *
 * Purpose:	Returns the number of MPI processes
 *
 * Return:	Success: non-negative
 *		Failure: negative
 *
 * Programmer:	Quincey Koziol
 *              Thursday, May 16, 2002
 *
 *-------------------------------------------------------------------------
 */
static int
H5FD__mpio_mpi_size(const H5FD_t *_file)
{
    const H5FD_mpio_t	*file = (const H5FD_mpio_t*)_file;

    FUNC_ENTER_STATIC

    /* Sanity checks */
    HDassert(file);
    HDassert(H5FD_MPIO == file->pub.driver_id);

    FUNC_LEAVE_NOAPI(file->mpi_size)
} /* end H5FD__mpio_mpi_size() */


/*-------------------------------------------------------------------------
 * Function:	H5FD__mpio_communicator
 *
 * Purpose:	Returns the MPI communicator for the file.
 *
 * Return:	Success:	The communicator
 *
 *		Failure:	NULL
 *
 * Programmer:	Robb Matzke
 *              Monday, August  9, 1999
 *
 *-------------------------------------------------------------------------
 */
static MPI_Comm
H5FD__mpio_communicator(const H5FD_t *_file)
{
    const H5FD_mpio_t	*file = (const H5FD_mpio_t*)_file;

    FUNC_ENTER_STATIC

    /* Sanity checks */
    HDassert(file);
    HDassert(H5FD_MPIO == file->pub.driver_id);

    FUNC_LEAVE_NOAPI(file->comm)
} /* end H5FD__mpio_communicator() */


/*-------------------------------------------------------------------------
 * Function:       H5FD__mpio_get_info
 *
 * Purpose:        Returns the file info of MPIO file driver.
 *
 * Returns:        Non-negative if succeed or negative if fails.
 *
 * Programmer:     John Mainzer
 *                 April 4, 2017
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
*/
static herr_t
H5FD__mpio_get_info(H5FD_t *_file, void **mpi_info)
{
    H5FD_mpio_t *file = (H5FD_mpio_t *)_file;
    herr_t      ret_value = SUCCEED;

    FUNC_ENTER_STATIC

    if(!mpi_info)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "mpi info not valid")

    *mpi_info = &(file->info);

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* H5FD__mpio_get_info() */

#endif /* H5_HAVE_PARALLEL */

