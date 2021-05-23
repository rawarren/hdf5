/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
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
 * Purpose: Public, shared definitions for Mirror VFD & remote Writer.
 */

#ifndef H5FDsubfiling_priv_H
#define H5FDsubfiling_priv_H

/********************/
/* Standard Headers */
/********************/

#include <assert.h>
#include <libgen.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/**************/
/* H5 Headers */
/**************/
#include "H5CXprivate.h" /* API Contexts                             */
#include "H5Dprivate.h"  /* Datasets                                 */
#include "H5Eprivate.h"  /* Error handling                           */
#include "H5Iprivate.h"  /* IDs                                      */
#include "H5Ipublic.h"
#include "H5MMprivate.h" /* Memory management                        */
#include "H5Pprivate.h"  /* Property lists                           */
#include "H5private.h"   /* Generic Functions                        */
// #include "H5FDioc.h"

#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    SELECT_IOC_ONE_PER_NODE = 0, /* Default */
    SELECT_IOC_EVERY_NTH_RANK,
    SELECT_IOC_WITH_CONFIG,
    SELECT_IOC_TOTAL,
    ioc_selection_options
} ioc_selection_t;

/****************************************************************************
 *
 * Structure: H5FD_subfiling_fapl_t
 *
 * Purpose:
 *
 *     H5FD_subfiling_fapl_t is a public structure that is used to pass 
 *     subfiling configuration data to the appropriate subfiling VFD via 
 *     the FAPL.  A pointer to an instance of this structure is a parameter 
 *     to H5Pset_fapl_subfiling() and H5Pget_fapl_subfiling().
 *
 * `magic`   (uint32_t)
 *
 *     Magic is a somewhat unique number which identifies this VFD from
 *     other VFDs.  Used in combination with a version number, we can
 *     validate a user generated file access property list (fapl).
 *     This field should be set to H5FD_SUBFILING_FAPL_T_MAGIC.
 *
 * `version` (uint32_t)
 *
 *     Version number of the H5FD_subfiling_fapl_t structure.  Any instance 
 *     passed to the above calls must have a recognized version number, or 
 *     an error will be flagged.
 *
 *     This field should be set to H5FD_CURR_SUBFILING_FAPL_T_VERSION.
 *
 ***   IO Concentrator Info ***
 ***   These fields will be replicated in the stacked IOC VFD which
 ***   provides the extended support for aggregating reads and writes
 ***   and allows global file access to node-local storage containers.
 *
 * `stripe_count` (int32_t)
 *
 *     The integer value which identifies the total number of 
 *     subfiles that have been algorthmically been selected to 
 *     to contain the segments of raw data which make up an HDF5
 *     file.  This value is used to implement the RAID-0 functionality
 *     when reading or writing datasets.
 *
 * `stripe_depth` (int64_t)
 *
 *     The stripe depth defines a limit on the maximum number of contiguous
 *     bytes that can be read or written in a single operation on any
 *     selected subfile.  Larger IO operations can exceed this limit
 *     by utilizing MPI derived types to construct an IO request which
 *     gathers additional data segments from memory for the IO request.
 *
 * `ioc_selection` (enum io_selection datatype)
 *
 *     The io_selection_t defines a specific algorithm by which IO
 *     concentrators (IOCs) and sub-files are identified.  The available
 *     algorthms are: SELECT_IOC_ONE_PER_NODE, SELECT_IOC_EVERY_NTH_RANK,
 *     SELECT_IOC_WITH_CONFIG, and SELECT_IOC_TOTAL.
 *
 ***   STACKING and other VFD support
 ***   i.e. FAPL caching
 ***
 *
 * `ioc_fapl_id` (hid_t)
 *
 *     A valid file access property list (fapl) is cached on each
 *     process and thus enables selection of an alternative provider
 *     for subsequent file operations.
 *     By defalt, Sub-filing employs an additional support VFD that
 *     provides file IO proxy capabilities to all MPI ranks in a
 *     distributed parallel application.  This IO indirection
 *     thus allows application access all sub-files even while
 *     these may actually be node-local and thus not directly
 *     accessable to remote ranks. 
 *
 ***   Subfiling file Info
 *
 * `subfile_dir`  char[]
 *
 *     A file directory name where subfiling files should be
 *     placed. Under normal circumstances, this directory name
 *     should match the directory path of the user defined HDF5
 *     file.
 *
 * `subfile_path` char[]
 *
 *     The full pathname of the user HDF5 file.
 *
 ****************************************************************************/

#ifndef H5FD_SUBFILING_FAPL_T_MAGIC
#define H5FD_CURR_SUBFILING_FAPL_T_VERSION     1
#define H5FD_SUBFILING_FAPL_T_MAGIC 0xFED01331
#endif

#ifndef H5FD_IOC_FAPL_T_MAGIC
#define H5FD_CURR_IOC_FAPL_T_VERSION     1
#define H5FD_IOC_FAPL_T_MAGIC 0xFED21331
#endif

#define H5FD_SUBFILING_PATH_MAX 4096

typedef struct config_common_t {
    uint32_t        magic;
    uint32_t        version;
    int32_t         stripe_count;
    int64_t         stripe_depth;
    ioc_selection_t ioc_selection;
	hid_t           ioc_fapl_id;
    int64_t         context_id;
    char            file_dir[H5FD_SUBFILING_PATH_MAX +1];
	char            file_path[H5FD_SUBFILING_PATH_MAX +1];
} config_common_t;

#define DRIVER_INFO_MESSAGE_MAX_INFO 65536
#define DRIVER_INFO_MESSAGE_MAX_LENGTH 65552 /* MAX_INFO + sizeof(info_header_t) */

#define K(n) ((n) *1024)
#define M(n) ((n) * (1024 * 1024))
#define H5FD_DEFAULT_STRIPE_DEPTH M(32)

typedef struct stat_record {
    int64_t op_count;
    double  min;
	double  max;
	double  total;
} stat_record_t;

typedef enum stat_category {
	WRITE_STAT = 0,
    WRITE_WAIT,
	READ_STAT,
    READ_WAIT,
	FOPEN_STAT,
	FCLOSE_STAT,
    QUEUE_STAT,
	TOTAL_STAT_COUNT
} stat_category_t;

typedef struct _info_header {
	uint8_t   version;
	uint8_t   unused_1;
	uint8_t   unused_2;
	uint8_t   unused_3;
	int32_t   info_length;
	char      vfd_key[8];
} info_header_t;
	

/* THE following definitions are used between H5FDsubfile_mpi.c
 * and H5FDsubfile_threads.c
 *
 * MPI Tags are 32 bits, we treat them as unsigned
 * to allow the use of the available bits for RPC
 * selections:
 *    0000
 *    0001 READ_OP  (Independent)
 *    0010 WRITE_OP (Independent)
 *    0011 /////////
 *    0100 CLOSE_OP (Independent)
 *    -----
 *    1000
 *    1001 COLLECTIVE_READ
 *    1010 COLLECTIVE_WRITE
 *    1011 /////////
 *    1100 COLLECTIVE_CLOSE
 *
 *   31    28      24      20      16      12       8       4       0|
 *   +-------+-------+-------+-------+-------+-------+-------+-------+
 *   |       |       |              ACKS             |      OP       |
 *   +-------+-------+-------+-------+-------+-------+-------+-------+
 *
 */

/* Bit 3 SET indicates collectives */
#    define COLL_FUNC (0x1 << 3)

#    define ACK_PART (0x0acc << 8)
#    define DATA_PART (0xd8da << 8)
#    define READY (0xfeed << 8)
#    define COMPLETED (0xfed1 << 8)

#    define READ_INDEP (READ_OP)
#    define READ_COLL (COLL_FUNC | READ_OP)
#    define WRITE_INDEP (WRITE_OP)
#    define WRITE_COLL (COLL_FUNC | WRITE_OP)

#    define WRITE_INDEP_ACK (ACK_PART | WRITE_OP)
#    define WRITE_INDEP_DATA (DATA_PART | WRITE_OP)

#    define READ_INDEP_DATA (DATA_PART | READ_OP)
#    define SET_LOGGING (LOGGING_OP)

#    define INT32_MASK 0x07FFFFFFFFFFFFFFF


typedef struct _info_message {
	info_header_t header;
	char msg_data[8];
} driver_info_t;

typedef enum io_ops {
    READ_OP = 1,
    WRITE_OP = 2,
    OPEN_OP = 3,
    CLOSE_OP = 4,
    FINI_OP = 8,
    LOGGING_OP = 16
} io_op_t;

typedef enum {
    SF_BADID = (-1),
    SF_TOPOLOGY = 1,
    SF_CONTEXT = 2,
    SF_NTYPES /* number of subfiling object types, MUST BE LAST */
} sf_obj_type_t;

typedef struct {
    long rank;
    long hostid;
} layout_t;


/* This typedef defines a fixed process layout which
 * can be reused for any number of file open operations
 */
typedef struct app_layout_t {
    long             hostid;
    layout_t       * layout;
    int            * node_ranks;
    int              node_count;
    int              node_index;
    int              local_peers;
    int              world_rank;
    int              world_size;
} app_layout_t;
		
/*  This typedef defines things related to IOC selections */
typedef struct topology {
    app_layout_t   * app_layout;    
    bool             rank_is_ioc;
    int              subfile_rank;
    int              n_io_concentrators;
    int            * io_concentrator;
    int            * subfile_fd;
    ioc_selection_t  selection_type;
} sf_topology_t;

typedef struct {
    hid_t            sf_context_id;
    hid_t            h5_file_id;
    int              sf_fid;
    size_t           sf_write_count;
    size_t           sf_read_count;
    size_t           sf_eof;
    /* Copy of the HDF5 File 'serial' number */
    unsigned long    fileno;
    int64_t          sf_stripe_size;
    int64_t          sf_blocksize_per_stripe;
    MPI_Comm         sf_msg_comm;
    MPI_Comm         sf_data_comm;
    MPI_Comm         sf_group_comm;
    MPI_Comm         sf_intercomm;
    int              sf_group_size;
    int              sf_group_rank;
    int              sf_intercomm_root;
    char *           subfile_prefix;
	char *           filename;
    sf_topology_t   *topology;
} subfiling_context_t;

typedef struct {
    /* {Datasize, Offset, FileID} */
    int64_t header[3];
    int     tag;
    int     source;
    int     subfile_rank;
    hid_t   context_id;
    double  start_time;
} sf_work_request_t;

typedef struct {
    hid_t h5_file_id;
    hid_t sf_context_id;
} file_map_to_context_t;

/* 
 * CAUTION::
 * Do we want or need this? 
 * Unfortunately, this structure is ONLY defined
 * in the H5FDsec2.c source file...
 */
typedef struct H5FD_sec2_t {
    H5FD_t          pub;        /* public stuff, must be first      */
    int             fd;         /* the filesystem file descriptor   */
    haddr_t         eoa;        /* end of allocated region          */
    haddr_t         eof;        /* end of file; current file size   */
    haddr_t         pos;        /* current file I/O position        */
    H5FD_file_op_t  op;         /* last operation                   */
    hbool_t         ignore_disabled_file_locks;
    char            filename[H5FD_MAX_FILENAME_LEN]; /* Copy of file name from open operation */
#ifndef H5_HAVE_WIN32_API
    /* On most systems the combination of device and i-node number uniquely
     * identify a file.  Note that Cygwin, MinGW and other Windows POSIX
     * environments have the stat function (which fakes inodes)
     * and will use the 'device + inodes' scheme as opposed to the
     * Windows code further below.
     */
    dev_t            device;    /* file device number   */
    ino_t            inode;     /* file i-node number   */
#else
    /* Files in windows are uniquely identified by the volume serial
     * number and the file index (both low and high parts).
     *
     * There are caveats where these numbers can change, especially
     * on FAT file systems.  On NTFS, however, a file should keep
     * those numbers the same until renamed or deleted (though you
     * can use ReplaceFile() on NTFS to keep the numbers the same
     * while renaming).
     *
     * See the MSDN "BY_HANDLE_FILE_INFORMATION Structure" entry for
     * more information.
     *
     * http://msdn.microsoft.com/en-us/library/aa363788(v=VS.85).aspx
     */
    DWORD nFileIndexLow;
    DWORD nFileIndexHigh;
    DWORD dwVolumeSerialNumber;

    HANDLE hFile; /* Native windows file handle */
#endif /* H5_HAVE_WIN32_API */

    /* Information from properties set by 'h5repart' tool
     *
     * Whether to eliminate the family driver info and convert this file to
     * a single file.
     */
    hbool_t fam_to_single;
} H5FD_sec2_t;

extern int sf_verbose_flag;
extern atomic_int sf_work_pending;

#ifdef __cplusplus
}
#endif

#endif	/* H5FDsubfiling_priv_H */
