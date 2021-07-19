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

#include "H5FDsubfiling.h"

static int sf_close_file_count = 0;
static int sf_ops_after_first_close = 0;
static int sf_enable_directIO = 0;

static int sf_write_ops = 0;
static double sf_pwrite_time = 0.0;
static double sf_write_wait_time = 0.0;

static int sf_read_ops = 0;
static double sf_pread_time = 0.0;
static double sf_read_wait_time = 0.0;
static double sf_queue_delay_time = 0.0;

/* The following is our basic template for a subfile filename.
 * Note that eventually we shouldn't use 0_of_N since we
 * intend to use the user defined HDF5 filename for a
 * zeroth subfile as well as for all metadata.
 */
#define SF_FILENAME_TEMPLATE "%ld_node_local_temp_%d_of_%d"
static const char *sf_filename_template = SF_FILENAME_TEMPLATE;

static int *request_count_per_rank = NULL;

atomic_int sf_workinprogress = 0;
atomic_int sf_work_pending = 0;
atomic_int sf_file_open_count = 0;
atomic_int sf_file_close_count = 0;
atomic_int sf_file_refcount = 0;
atomic_int sf_ioc_fini_refcount = 0;
atomic_int sf_ioc_ready = 0;
volatile int sf_shutdown_flag = 0;

/* 
 * Structure definitions to enable async io completions
 * We first define a structure which contains the basic
 * input arguments for the functions which were originally
 * invoked.  See below.
 */
typedef struct _client_io_args {
    int        ioc;            /* ID of the IO Concentrator handling this IO.   */
    hid_t      context_id;     /* The context id provided for the read or write */
    int64_t    offset;         /* The file offset for the IO operation          */
    int64_t    elements;       /* How many bytes                                */
    void      *data;           /* A pointer to the (contiguous) data segment    */
    MPI_Request io_req;        /* An MPI request to allow the code to loop while */
                               /* making progress on multiple IOs               */
} io_args_t;

/* pre-define */
typedef struct _client_io_func io_func_t;

struct _client_io_func {
    int (*io_function)(void *this_io); /* pointer to a completion function */
    io_args_t io_args;         /* arguments passed to the completion function   */
    int pending;               /* The function is complete (0) or pending (1)?  */
};

typedef struct _io_req {
    struct _io_req   *prev;    /* A simple list structure containing completion */
    struct _io_req   *next;    /* functions. These should get removed as IO ops */
    io_func_t completion_func; /* are completed */
} io_req_t;


int      n_io_pending = 0;
io_req_t pending_io_requests;



typedef struct _client_xfer_info {
    int64_t    offset;
    int64_t    length;
    int        ioc_targets;
    io_op_t    op;
} client_xfer_info_t;


typedef struct _xfer_info {
    int64_t    offset;
    int64_t    length;
} xfer_info_t;

#define STAT_BLOCKSIZE 1024
typedef struct _ioc_stats {
    int         read_index;
    int         read_size;
    xfer_info_t *read_info;
    int         write_index;
    int         write_size;
    xfer_info_t *write_info;
} ioc_stats_t;

static ioc_stats_t ioc_xfer_records;


int client_op_index = 0;
int client_op_size = 0;
client_xfer_info_t *client_ops = NULL;

/* const char *sf_subfile_prefix = "."; */

#define MAX_WORK_PER_RANK 2
#    define K(n) ((n) *1024)
#    define M(n) ((n) * (1024 * 1024))
#    define DEFAULT_STRIPE_SIZE M(32)
#    define MAX_DEPTH 1024


/*
=========================================
Private functions
=========================================
*/

static char * get_ioc_subfile_path(int ioc, int ioc_count, subfiling_context_t *sf_context);
static int async_completion(void *arg);

/* ===================================================================== */
/* MPI_Datatype Creation functions.
 * These are catagorized by usage paterns, i.e. when data is sent to or
 * received from and IOC, the initial data offset provided by the user
 * may or may NOT start on a stripe boundary.  Because this, the initial
 * data segment to the selected IOC will often be less than 'stripe_size'
 * in length.  The purpose of these Datatype creation functions is to
 * enable the gathering of all data from this client to the IOC target
 * into a single MPI message.  The MPI datatype will the be utilized by
 * the sending function to pack data into a contiguous block of memory
 * which enables the IOC to write to disk in an effective manner.
 * ===================================================================== */

/*-------------------------------------------------------------------------
 * Function:    H5FD__create_first_mpi_type
 *
 * Purpose:     Return an appropriate MPI datatype to represent the initial
 *              IO operation when reading or writing data to or from an IO
 *              Concentrator (IOC).
 *
 *              If the 'first_io' is sufficient to complete the IO to the
 *              IOC, then the returned MPI datatype will simply be MPI_BYTE.
 *              For all other non-zero length IO operations, we create a
 *              derived MPI datatype using MPI_Type_indexed. The 'ioc_depth'
 *              input will define the number of blocks/disps pairs that are
 *              required to represent the desired IO operation.
 *
 * Return:      The MPI_Datatype that will be used to send or receive data.
 * Errors:      MPI_Type_NULL if for any reason, the MPI_Datatype creation
 *              fails.
 *
 * Programmer:  Richard Warren
 *              7/17/2020
 *
 * Changes:     Initial Version/None.
 *
 *-------------------------------------------------------------------------
 */
#ifdef H5_ALLOW_VFD_DERIVED_TYPES
static MPI_Datatype
H5FD__create_first_mpi_type(subfiling_context_t *context, int ioc_depth,
    int64_t offset, int64_t target_write_bytes, int64_t first_io)
{
    MPI_Datatype newType = MPI_DATATYPE_NULL;
    int64_t      stripe_size = context->sf_stripe_size;
    int64_t      blocksize_per_stripe = context->sf_blocksize_per_stripe;
    int64_t      offset_in_stripe = offset % stripe_size;
    int64_t      next_offset = blocksize_per_stripe - offset_in_stripe;
    int64_t      total_bytes = first_io;

    if (first_io == target_write_bytes) {
        if (first_io > 0) {
            return MPI_BYTE;
        }
    }
    if (first_io) {
        int  k;
        int  temp_blocks[64];
        int  temp_disps[64];
        int *blocks = temp_blocks;
        int *disps = temp_disps;
        if (ioc_depth > 64) {
            blocks = (int *) calloc((size_t) ioc_depth, sizeof(int));
            if (blocks == NULL) {
                perror("calloc");
                return newType;
            }
            disps = (int *) calloc((size_t) ioc_depth, sizeof(int));
            if (disps == NULL) {
                perror("calloc");
                return newType;
            }
        }
        blocks[0] = (int) first_io;
        disps[0] = (int) 0;
        for (k = 1; k <= ioc_depth; k++) {
            disps[k] = (int) next_offset;
            blocks[k] = (int) stripe_size;
            total_bytes += stripe_size;
            next_offset += context->sf_blocksize_per_stripe;
        }
        if (total_bytes != target_write_bytes) {
            printf("Warning (%s): total_SUM(%ld) != target_bytes(%ld)\n",
                __func__, total_bytes, target_write_bytes);
        }

        if (MPI_Type_indexed(k, blocks, disps, MPI_BYTE, &newType) !=
            MPI_SUCCESS) {
            perror("MPI_Type_indexed failed!");
            return newType;
        }
        MPI_Type_commit(&newType);
        if (1) {
            int type_size;
            MPI_Type_size(newType, &type_size);
            if (type_size != target_write_bytes) {
                printf("%s: type_size=%d should be: %ld\n", __func__, type_size,
                    target_write_bytes);
            }
        }
        if (ioc_depth > 64) {
            if (blocks != temp_blocks) {
                free(blocks);
                blocks = NULL;
            }
            if (disps != temp_disps) {
                free(disps);
                disps = NULL;
            }
        }
    }
    return newType;
} /* end H5FD__create_first_mpi_type() */
#else

/* Fill the output vectors 'io_offset', 'io_datasize' and 'io_f_offset'
 * All calculations are in terms of bytes.
*/
static void
H5FD__create_first_mpi_type(subfiling_context_t *context, int ioc_depth,
                            int64_t src_offset, int64_t target_datasize, int64_t f_offset,
                            int64_t *io_offset, int64_t *io_datasize, int64_t *io_f_offset,
                            int64_t first_io)
{
    int64_t      stripe_size = context->sf_stripe_size;
    int64_t      blocksize_per_stripe = context->sf_blocksize_per_stripe;
    int64_t      offset_in_stripe = f_offset % stripe_size;
    int64_t      next_offset = blocksize_per_stripe - offset_in_stripe;
    int64_t      total_bytes = first_io;

    io_offset[0] = src_offset;
    io_datasize[0] = first_io;
    io_f_offset[0] = f_offset;

    if (first_io == target_datasize) {
        return;
    }
    if (first_io) {
        int  k;
        f_offset += first_io;
        for (k = 1; k <= ioc_depth; k++) {
            io_offset[k] = next_offset;
            io_datasize[k] = stripe_size;
            io_f_offset[k] = f_offset;
            f_offset += stripe_size;
            total_bytes += stripe_size;
            next_offset += context->sf_blocksize_per_stripe;
        }
        if (total_bytes != target_datasize) {
            printf("Warning (%s): total_SUM(%ld) != target_bytes(%ld)\n",
                __func__, total_bytes, target_datasize);
        }
    }
    return;
} /* end H5FD__create_first_mpi_type() */

#endif

#ifdef H5_ALLOW_VFD_DERIVED_TYPES
/*-------------------------------------------------------------------------
 * Function:    H5FD__create_final_mpi_type
 *
 * Purpose:     Return an appropriate MPI datatype to represent the final
 *              IO operation when reading or writing data to or from an IO
 *              Concentrator (IOC).
 *
 *              The data that we're sending to an IO concentrator (IOC)
 *              contains the final collection of bytes. Other than that detail,
 *              this is pretty much like the typical' IO case, i.e. all block
 *              sizes are identical (execpt for the very last block).
 *Furthermore, they all start at relative stripe offset of 0, in other words on
 *a 'stripe_size' boundary.
 *
 * Return:      The MPI_Datatype that will be used to send or receive data.
 * Errors:      MPI_Type_NULL if for any reason, the MPI_Datatype creation
 *              fails.
 *
 * Programmer:  Richard Warren
 *              7/17/2020
 *
 * Changes:     Initial Version/None.
 *
 *-------------------------------------------------------------------------
 */
static MPI_Datatype
H5FD__create_final_mpi_type(subfiling_context_t *context, int ioc_depth,
    int64_t target_write_bytes, int64_t last_write)
{
    MPI_Datatype newType = MPI_DATATYPE_NULL;
    int64_t      stripe_size = context->sf_stripe_size;
    int64_t      depth_in_bytes = (stripe_size * ioc_depth) + last_write;
    int64_t      total_bytes = last_write;

    if (depth_in_bytes == target_write_bytes) {
        if (depth_in_bytes > 0) {
            return MPI_BYTE;
        }
    }

    if (depth_in_bytes) {
        int  k;
        int  temp_blocks[64];
        int  temp_disps[64];
        int *blocks = temp_blocks;
        int *disps = temp_disps;
        if (ioc_depth > 64) {
            blocks = (int *) calloc((size_t) ioc_depth, sizeof(int));
            if (blocks == NULL) {
                perror("calloc");
                return newType;
            }
            disps = (int *) calloc((size_t) ioc_depth, sizeof(int));
            if (disps == NULL) {
                perror("calloc");
                return newType;
            }
        }

        for (k = 0; k < ioc_depth; k++) {
            disps[k] = (int) (k * context->sf_blocksize_per_stripe);
            blocks[k] = (int) stripe_size;
            total_bytes += stripe_size;
        }
        blocks[k - 1] = (int) last_write;
        if (total_bytes != target_write_bytes) {
            printf("Warning (%s): total_SUM(%ld) != target_bytes(%ld)\n",
                __func__, total_bytes, target_write_bytes);
        }

        if (MPI_Type_indexed(ioc_depth, blocks, disps, MPI_BYTE, &newType) !=
            MPI_SUCCESS) {
            return MPI_DATATYPE_NULL;
        }
        MPI_Type_commit(&newType);
        if (ioc_depth > 64) {
            if (blocks != temp_blocks) {
                free(blocks);
                blocks = NULL;
            }
            if (disps != temp_disps) {
                free(disps);
                disps = NULL;
            }
        }
    }
    return newType;
}
#else  /* #ifdef H5_ALLOW_VFD_DERIVED_TYPES */

/* Fill the output vectors 'io_offset', 'io_datasize' and 'io_f_offset'
 * All calculations are in terms of bytes.
*/
static void
H5FD__create_final_mpi_type(subfiling_context_t *context, int ioc_depth,
                            int64_t src_offset, int64_t target_datasize, int64_t f_offset,
                            int64_t *io_offset, int64_t *io_datasize, int64_t *io_f_offset,
                            int64_t last_io)
{
    int64_t      stripe_size = context->sf_stripe_size;
    int64_t      blocksize_per_stripe = context->sf_blocksize_per_stripe;
    int64_t      offset_in_stripe = f_offset % stripe_size;
    int64_t      next_offset = src_offset;
    int64_t      total_bytes = 0;
    int64_t      overcount = stripe_size - last_io;

    if (last_io == target_datasize) {
        io_offset[0] = src_offset;
        io_f_offset[0] = f_offset;
        io_datasize[0] = last_io;
        return;
    }

    if (last_io) {
        int  k;
        for (k = 0; k < ioc_depth; k++) {
            io_offset[k] = next_offset;
            io_datasize[k] = stripe_size;
            io_f_offset[k] = f_offset;
            f_offset += stripe_size;
            total_bytes += stripe_size;
            next_offset += context->sf_blocksize_per_stripe;
        }
        io_datasize[ioc_depth -1] = last_io;
        total_bytes -= overcount;

        if (total_bytes != target_datasize) {
            printf("Warning (%s): total_SUM(%ld) != target_bytes(%ld)\n",
                __func__, total_bytes, target_datasize);
        }
    }
    return;
} /* end H5FD__create_final_mpi_type() */

#endif /* #ifdef H5_ALLOW_VFD_DERIVED_TYPES */


#ifdef H5_ALLOW_VFD_DERIVED_TYPES
/*-------------------------------------------------------------------------
 * Function:    H5FD__create_f_l_mpi_type
 *
 * Purpose:     Return an appropriate MPI datatype which includes both the
 *              first and final IO data segments.
 *
 *              A special case where the current IOC has both the first and
 *              final write blocks. This function is basically a merge of
 *              the first_mpi_type and final_mpi_type functions.
 *
 * Return:      The MPI_Datatype that will be used to send or receive data.
 * Errors:      MPI_Type_NULL if for any reason, the MPI_Datatype creation
 *              fails.
 *
 * Programmer:  Richard Warren
 *              7/17/2020
 *
 * Changes:     Initial Version/None.
 *
 *-------------------------------------------------------------------------
 */
static MPI_Datatype
H5FD__create_f_l_mpi_type(subfiling_context_t *context, int ioc_depth,
    int64_t offset, int64_t target_write_bytes, int64_t first_write,
    int64_t last_write)
{
    MPI_Datatype newType = MPI_DATATYPE_NULL;
    int64_t      stripe_size = context->sf_stripe_size;
    int64_t      blocksize_per_stripe = context->sf_blocksize_per_stripe;
    int64_t      offset_in_stripe = offset % stripe_size;
    int64_t      next_offset = blocksize_per_stripe - offset_in_stripe;
    int64_t      total_bytes = first_write + last_write;
    int          sf_world_rank = context->topology->app_layout->world_rank;
    int          status = 0;

    /* We might actaully check that the 'target_write_bytes'
     * input variable exceeds 2Gb.  If so, then we should
     * always create a derived type.
     */
    if ((total_bytes == target_write_bytes) &&
        (context->topology->n_io_concentrators == 1)) {
        return MPI_BYTE;
    } else if (first_write) {
        int  k;
        int  temp_blocks[64];
        int  temp_disps[64];
        int *blocks = temp_blocks;
        int *disps = temp_disps;
        if (ioc_depth > 64) {
            blocks = (int *) calloc((size_t) ioc_depth, sizeof(int));
            if (blocks == NULL) {
                perror("calloc");
                return newType;
            }
            disps = (int *) calloc((size_t) ioc_depth, sizeof(int));
            if (disps == NULL) {
                perror("calloc");
                return newType;
            }
        }

        blocks[0] = (int) first_write;
        disps[0] = 0;

        for (k = 1; k < ioc_depth; k++) {
            blocks[k] = (int) stripe_size;
            disps[k] = (int) next_offset;
            next_offset += context->sf_blocksize_per_stripe;
            total_bytes += stripe_size;
        }
        if (k == 1) {
            disps[k] = (int) next_offset;
        }
        blocks[k] = (int) last_write;

        if (total_bytes != target_write_bytes) {
            printf("[%d] Warning (%s): total_SUM(%ld) != target_bytes(%ld)\n",
                sf_world_rank, __func__, total_bytes, target_write_bytes);
        }

        status = MPI_Type_indexed(k + 1, blocks, disps, MPI_BYTE, &newType);
        if (status != MPI_SUCCESS) {
            puts("MPI_Type_indexed failed!");
            return MPI_DATATYPE_NULL;
        }

        status = MPI_Type_commit(&newType);
        if (ioc_depth > 64) {
            if (blocks != temp_blocks) {
                free(blocks);
                blocks = NULL;
            }
            if (disps != temp_disps) {
                free(disps);
                disps = NULL;
            }
        }

        if (status != MPI_SUCCESS) {
            MPI_Type_free(&newType);
            puts("MPI_Type_commit failed!");
            return MPI_DATATYPE_NULL;
        }
    }
    return newType;
} /* end H5FD__create_f_l_mpi_type() */

#else  /* #ifdef H5_ALLOW_VFD_DERIVED_TYPES */
static void
H5FD__create_f_l_mpi_type(subfiling_context_t *context, int ioc_depth,
                          int64_t src_offset, int64_t target_datasize, int64_t f_offset,
                          int64_t *io_offset, int64_t *io_datasize, int64_t *io_f_offset,
                          int64_t first_io, int64_t last_io )
{
    int64_t      stripe_size = context->sf_stripe_size;
    int64_t      blocksize_per_stripe = context->sf_blocksize_per_stripe;
    int64_t      offset_in_stripe = f_offset % stripe_size;
    int64_t      next_offset = src_offset + (blocksize_per_stripe - offset_in_stripe);
    int64_t      total_bytes = first_io;
    int          sf_world_rank = context->topology->app_layout->world_rank;

    io_offset[0] = src_offset;
    io_datasize[0] = first_io;
    io_f_offset[0] = f_offset;

    if (total_bytes == target_datasize) {
        return;
    }

    if (total_bytes) {
        int  k;
        f_offset += (blocksize_per_stripe - f_offset);
        for (k = 1; k < ioc_depth; k++) {
            io_offset[k] = next_offset;
            io_datasize[k] = stripe_size;
            io_f_offset[k] = f_offset;
            total_bytes += stripe_size;
            f_offset += context->sf_blocksize_per_stripe;
            next_offset += context->sf_blocksize_per_stripe;
        }
        io_datasize[ioc_depth] = last_io;
        io_f_offset[ioc_depth] = f_offset;
        io_offset[ioc_depth] = next_offset;

        total_bytes += last_io;

        if (total_bytes != target_datasize) {
            printf("[%d] Warning (%s): total_SUM(%ld) != target_bytes(%ld)\n",
                sf_world_rank, __func__, total_bytes, target_datasize);
        }
    }
    return;
} /* end H5FD__create_f_l_mpi_type() */

#endif /* #ifdef H5_ALLOW_VFD_DERIVED_TYPES */

#ifdef H5_ALLOW_VFD_DERIVED_TYPES
/*-------------------------------------------------------------------------
 * Function:    H5FD__create_mpi_uniform_type
 *
 * Purpose:     Return an appropriate MPI datatype to represent the typical
 *              IO operation when reading or writing data to or from an IO
 *              Concentrator (IOC).
 *
 *              Each data segment is of 'stripe_size' length and will be
 *              seperated from a previous or following segment by
 *              'sf_blocksize_per_stripe' bytes of data.
 *
 * Return:      The MPI_Datatype that will be used to send or receive data.
 * Errors:      MPI_Type_NULL if for any reason, the MPI_Datatype creation
 *              fails.
 *
 * Programmer:  Richard Warren
 *              7/17/2020
 *
 * Changes:     Initial Version/None.
 *
 *-------------------------------------------------------------------------
 */
static MPI_Datatype
H5FD__create_mpi_uniform_type(
    subfiling_context_t *context, int ioc_depth, int64_t target_write_bytes)
{
    MPI_Datatype newType = MPI_DATATYPE_NULL;
    int64_t      stripe_size = context->sf_stripe_size;
    int64_t      check_depth = stripe_size * ioc_depth;
    int64_t      total_bytes = 0;

    if (check_depth == stripe_size) {
        if (target_write_bytes > 0) {
            return MPI_BYTE;
        }
    }

    if (target_write_bytes) {
        int  k;
        int  temp_blocks[64];
        int  temp_disps[64];
        int *blocks = temp_blocks;
        int *disps = temp_disps;
        if (ioc_depth > 64) {
            blocks = (int *) calloc((size_t) ioc_depth, sizeof(int));
            if (blocks == NULL) {
                perror("calloc");
                return newType;
            }
            disps = (int *) calloc((size_t) ioc_depth, sizeof(int));
            if (disps == NULL) {
                perror("calloc");
                return newType;
            }
        }
        for (k = 0; k < ioc_depth; k++) {
            disps[k] = (int) (k * context->sf_blocksize_per_stripe);
            blocks[k] = (int) (stripe_size);
            total_bytes += stripe_size;
        }

        if (total_bytes != target_write_bytes) {
            printf("Warning (%s): total_SUM(%ld) != target_bytes(%ld)\n",
                __func__, total_bytes, target_write_bytes);
        }

        if (MPI_Type_indexed(ioc_depth, blocks, disps, MPI_BYTE, &newType) !=
            MPI_SUCCESS) {
            perror("MPI_Type_indexed failed!");
            return MPI_DATATYPE_NULL;
        }
        MPI_Type_commit(&newType);
        if (1) {
            int type_size;
            MPI_Type_size(newType, &type_size);
            if (type_size != target_write_bytes) {
                printf("%s: type_size=%d should be: %ld\n", __func__, type_size,
                    target_write_bytes);
            }
        }

        if (ioc_depth > 64) {
            if (blocks != temp_blocks) {
                free(blocks);
                blocks = NULL;
            }
            if (disps != temp_disps) {
                free(disps);
                disps = NULL;
            }
        }
    }
    return newType;
} /* end H5FD__create_mpi_uniform_type() */

#else  /* #ifdef H5_ALLOW_VFD_DERIVED_TYPES */

static void
H5FD__create_mpi_uniform_type(subfiling_context_t *context, int ioc_depth,
                              int64_t src_offset, int64_t target_datasize, int64_t f_offset,
                              int64_t *io_offset, int64_t *io_datasize, int64_t *io_f_offset)
{
    int64_t      stripe_size = context->sf_stripe_size;
    int64_t      blocksize_per_stripe = context->sf_blocksize_per_stripe;
    int64_t      next_offset = src_offset + blocksize_per_stripe;
    int64_t      total_bytes = 0;

    io_offset[0] = src_offset;
	
    io_datasize[0] = stripe_size;
    io_f_offset[0] = f_offset;
    if (target_datasize == 0) {
		io_datasize[0] = 0;
        return;
    }
	else io_datasize[0] = stripe_size;
    f_offset += blocksize_per_stripe;
    total_bytes = stripe_size;

    if (target_datasize > stripe_size) {
        int  k;

        for (k = 1; k < ioc_depth; k++) {
            io_offset[k] = next_offset;
            io_datasize[k] = stripe_size;
            io_f_offset[k] = f_offset;
            total_bytes += stripe_size;
            f_offset += blocksize_per_stripe;
            next_offset += blocksize_per_stripe;
        }

        if (total_bytes != target_datasize) {
            printf("Warning (%s): total_SUM(%ld) != target_bytes(%ld)\n",
                __func__, total_bytes, target_datasize);
        }
    }
    return;
} /* end H5FD__create_mpi_uniform_type() */

#endif /* #ifdef H5_ALLOW_VFD_DERIVED_TYPES */

#ifdef H5_ALLOW_VFD_DERIVED_TYPES
/*-------------------------------------------------------------------------
 * Function:    init__indep_io
 *
 * Purpose:     Utility function to initialize the set of IO transactions
 *              used to communicate with IO concentrators for read and write
 *              IO operations.
 *
 * Return:      A filled set of vectors (1 entry per IO concentrator) which
 *              fully describe the IO transactions for read and writes.
 *              At most, every IO concentrator will have a descriptor which
 *              identifies the local memory offset, the virtual FILE offset,
 *              and the total length of the IO which will be sent to or
 *              received from the individual IOCs.
 *
 *              For IO operations which involve a subset of IO concentrators,
 *              the vector entries for the unused IOCs will have lengths of
 *              zero and MPI NULL datatypes.
 *
 * Errors:      Cannot fail.
 *
 * Programmer:  Richard Warren
 *              7/17/2020
 *
 * Changes:     Initial Version/None.
 *
 *-------------------------------------------------------------------------
 */
int
init__indep_io(void *_sf_context, int64_t *sf_source_data_offset,
    int64_t *sf_datasize, int64_t *sf_offset, MPI_Datatype *sf_dtype,
    int64_t offset, int64_t elements, int dtype_extent)
{
    subfiling_context_t *sf_context = _sf_context;
    int     container_count = sf_context->topology->n_io_concentrators;
    int64_t stripe_size = sf_context->sf_stripe_size;
    int64_t data_size = elements * dtype_extent;

    int64_t start_id = offset / stripe_size;
    int64_t offset_in_stripe = offset % stripe_size;
    int64_t start_length = MIN(data_size, (stripe_size - offset_in_stripe));
    int64_t start_row = start_id / container_count;
    int64_t ioc_start = start_id % container_count;

    int64_t final_offset = offset + data_size;
    int64_t final_id = final_offset / stripe_size;
    int64_t final_length =
        (start_length == data_size ? 0 : final_offset % stripe_size);
    int64_t ioc_final = final_id % container_count;
    int64_t container_bytes = 0, total_bytes = 0;
    int64_t source_offset = 0;

    int     row_id_start = (int) (start_id - ioc_start);
    int     row_id_final = (int) (final_id - ioc_final);
    int     i, k, depth = ((row_id_final - row_id_start) / container_count) + 1;
    int     container_id = (int) start_id;
    int64_t row_offset = (int64_t)(start_row * stripe_size);

    for (i = 0, k = (int) ioc_start; i < container_count; i++) {
        int container_depth = depth;
        hbool_t is_first = false, is_last = false;
        container_bytes = 0;
        sf_datasize[k] = container_bytes;
        if (total_bytes < data_size) {
            if (k == ioc_start) {
                is_first = true;
                container_bytes = start_length;
                container_depth--; /* Account for the start_length */
                if (ioc_final < ioc_start) {
                    container_depth--;
                    depth--;
                }
            }
            if (k == ioc_final) {
                is_last = true;
                container_bytes += final_length;
                if (container_depth)
                    container_depth--; /* Account for the final_length */
                if (depth)
                    depth--;
            }
            container_bytes += container_depth * stripe_size;
            total_bytes += container_bytes;
        }

        sf_source_data_offset[k] = source_offset;
        sf_datasize[k] = container_bytes;
        sf_offset[k] = row_offset + offset_in_stripe;

        if (container_count == 1) {
            sf_dtype[k] = MPI_BYTE;
        } else {
            /* Fill the IO datatypes */
            if (is_first) {
                if (is_last) { /* First + Last */
                    sf_dtype[k] = H5FD__create_f_l_mpi_type(sf_context,
                        container_depth + 1, sf_offset[k], container_bytes,
                        start_length, final_length);
                } else { /* First ONLY */
                    sf_dtype[k] =
                        H5FD__create_first_mpi_type(sf_context, container_depth,
                            sf_offset[k], container_bytes, start_length);
                }
                source_offset += start_length;
                offset_in_stripe = 0;
            } else if (is_last) { /* Last ONLY */
                source_offset += stripe_size;
                sf_dtype[k] = H5FD__create_final_mpi_type(
                    sf_context, container_depth, container_bytes, final_length);
            } else { /* Everything else (uniform) */
                source_offset += stripe_size;
                sf_dtype[k] = H5FD__create_mpi_uniform_type(
                    sf_context, container_depth, container_bytes);
            }
        }
        k++;
        container_id++;

        if (k == container_count) {
            k = 0;
            depth = ((row_id_final - container_id) / container_count) + 1;
            row_offset += stripe_size;
        }
    }
    if (total_bytes != data_size) {
        printf("Error: total_bytes != data_size\n");
    }

    return container_count;
} /* end init__indep_io() */

#else  /* not H5_ALLOW_VFD_DERIVED_TYPES */
/*-------------------------------------------------------------------------
 * Function:    init__indep_io
 *
 * Purpose:     Utility function to initialize the set of IO transactions
 *              used to communicate with IO concentrators for read and write
 *              IO operations.
 *
 * Return:      A filled set of vectors.  As a consequence of not allowing
 *              use of MPI derived datatypes in the VFD layer, we need to
 *              accomodate the possiblity that large IO transactions will
 *              be required to use multiple IOs per IOC.  
 *
 *              Example: Using 4 IOCs, each with 1M stripe-depth; when
 *              presented an IO request for 8MB then at a minimum each IOC
 *              will require 2 IOs of 1MB each.  Depending on the starting
 *              file offset, the 2 IOs can instead be 3...
 *              
 *              To fully describe the IO transactions for read and writes we
 *              we thus use a return type where each IOC vector element is
 *              instead a vector itself and has a vector length of which
 *              corresponds to the max number of IO transactions per IOC.
 *              In the example above, these vector lengths can be 2 or 3.
 *              The actual length is determined by the 'container_depth'
 *              variable.
 *
 *              For IO operations which involve a subset of IO concentrators,
 *              the vector entries for the unused IOCs will have lengths of
 *              zero and MPI NULL datatypes.  The 'container_depth' in this
 *              case will always be 1.
 *
 * Return value: The vector "depth" or max number of IOs per IOC.
 *
 * Errors:      Cannot fail.
 *
 * Programmer:  Richard Warren
 *              7/17/2020
 *
 * Changes:     Initial Version/None.
 *
 *-------------------------------------------------------------------------
 */
#if 1
int
init__indep_io(void *_sf_context, size_t maxdepth, int ioc_total,
			   int64_t *sf_source_data_offset,
               int64_t *sf_datasize,
			   int64_t *sf_offset,
               int64_t offset, int64_t elements, int dtype_extent)
#else
int
init__indep_io(void *_sf_context, size_t maxdepth, int ioc_total,
			   int64_t sf_source_data_offset[maxdepth][ioc_total],
               int64_t sf_datasize[maxdepth][ioc_total],
			   int64_t sf_offset[maxdepth][ioc_total],
               int64_t offset, int64_t elements, int dtype_extent)
#endif	
{
    subfiling_context_t *sf_context = _sf_context;
    int     container_count = sf_context->topology->n_io_concentrators;
    int64_t stripe_size = sf_context->sf_stripe_size;
    int64_t data_size = elements * dtype_extent;

    int64_t start_id = offset / stripe_size;
    int64_t offset_in_stripe = offset % stripe_size;
    int64_t start_length = MIN(data_size, (stripe_size - offset_in_stripe));
    int64_t start_row = start_id / container_count;
    int64_t ioc_start = start_id % container_count;

    int64_t final_offset = offset + data_size;
    int64_t final_id = final_offset / stripe_size;
    int64_t final_length =
        (start_length == data_size ? 0 : final_offset % stripe_size);
    int64_t ioc_final = final_id % container_count;
    int64_t container_bytes = 0, total_bytes = 0;
    int64_t source_offset = 0;

    int     row_id_start = (int) (start_id - ioc_start);
    int     row_id_final = (int) (final_id - ioc_final);
    int     i, k, depth = ((row_id_final - row_id_start) / container_count) + 1;
    int     container_id = (int) start_id;
    int64_t row_offset = (int64_t)(start_row * stripe_size);

    /* Given the IO parameters, we loop thru the set of IOCs
     * to determine the various vector components for each.
     * Those IOCs whose datasize is zero (0), will not have
     * IO requests passed to them.
     */
    for (i = 0, k = (int) ioc_start; i < container_count; i++) {
        int index;
		/* We use 'output_offset' as an index into a linear 
         * version of a 2D array. In 'C' the last subscript
         * is the one that varies most rapidly.
         * In our case, the 2D array is represented as
         * array[ maxdepth ][ container_count ]
         */
        int output_offset = k;	
        int container_depth = depth;

        hbool_t is_first = false, is_last = false;
        int64_t __sf_source_data_offset[maxdepth];
        int64_t __sf_datasize[maxdepth];
        int64_t __sf_offset[maxdepth];

        memset(__sf_source_data_offset, 0, sizeof(__sf_offset));
        memset(__sf_datasize, 0, sizeof(__sf_offset));
        memset(__sf_offset, 0, sizeof(__sf_offset));

        container_bytes = 0;
        __sf_datasize[k] = container_bytes;
        if (total_bytes < data_size) {
            if (k == ioc_start) {
                is_first = true;
                container_bytes = start_length;
                container_depth--; /* Account for the start_length */
                if (ioc_final < ioc_start) {
                    container_depth--;
                    depth--;
                }
            }
            if (k == ioc_final) {
                is_last = true;
                container_bytes += final_length;
                if (container_depth)
                    container_depth--; /* Account for the final_length */
                if (depth)
                    depth--;
            }
            container_bytes += container_depth * stripe_size;
            total_bytes += container_bytes;
        }

        __sf_source_data_offset[0] = source_offset;
        __sf_datasize[0] = container_bytes;
        __sf_offset[0] = row_offset + offset_in_stripe;

        if (container_count == 1) {
            // sf_dtype[k] = MPI_BYTE;
        } else {
            /* Fill the IO datatypes */
            if (is_first) {
                if (is_last) { /* First + Last */
                    H5FD__create_f_l_mpi_type(sf_context,
                                              container_depth + 1,
											  source_offset,
											  container_bytes,
											  row_offset + offset_in_stripe,
                                              __sf_source_data_offset,
                                              __sf_datasize,
                                              __sf_offset,
                                              start_length,
                                              final_length );
                } else { /* First ONLY */
                    H5FD__create_first_mpi_type(sf_context,
                                                container_depth,
                                                source_offset,
                                                container_bytes,
                                                row_offset + offset_in_stripe,
                                                __sf_source_data_offset,
                                                __sf_datasize,
                                                __sf_offset,
                                                start_length);

                }
                /* Move the memory pointer to the starting location
                 * for next IOC request.
                 */
                source_offset += start_length;
                // offset_in_stripe = 0;
            } else if (is_last) { /* Last ONLY */
                H5FD__create_final_mpi_type( sf_context,
                                             container_depth,
                                             source_offset,
                                             container_bytes,
                                             row_offset + offset_in_stripe,
                                             __sf_source_data_offset,
                                             __sf_datasize,
                                             __sf_offset,
                                             final_length);
                /* Probably not needed... */
                source_offset += stripe_size;
            } else { /* Everything else (uniform) */
                H5FD__create_mpi_uniform_type(sf_context,
                                              container_depth,
											  source_offset,
                                              container_bytes,
											  row_offset + offset_in_stripe,
											  __sf_source_data_offset,
											  __sf_datasize,
											  __sf_offset );
                source_offset += stripe_size;
            }
        }

        /*  Copy the per ioc subvector values into the outputs */
        for (index = 0; index <= depth; index++) {
            /* The output 2D arrays are indexed by [depth][ioc]
             * to index by IOC (1 of 'container_count' containers)
             * we simply add the container_count after each loop
             * to arrive at the array info for the next container.
             */
#if 1
            sf_source_data_offset[output_offset] = __sf_source_data_offset[index];
            sf_datasize[output_offset]           = __sf_datasize[index];
            sf_offset[output_offset]             = __sf_offset[index];
            output_offset += container_count;
#else
            sf_source_data_offset[index][k] = __sf_source_data_offset[index];
            sf_datasize[index][k]           = __sf_datasize[index];
            sf_offset[index][k]             = __sf_offset[index];
#endif
        }

        k++;
		offset_in_stripe += __sf_datasize[0];
        container_id++;

        if (k == container_count) {
            k = 0;
            depth = ((row_id_final - container_id) / container_count) + 1;
            row_offset += stripe_size;
        }
    }
    if (total_bytes != data_size) {
        printf("Error: total_bytes != data_size\n");
    }

    return depth +1;
} /* end init__indep_io() */
#endif

/*-------------------------------------------------------------------------
 * Function:    Internal read__independent_async
 *
 * Purpose:     The IO operations can be striped across a selection of
 *              IO concentrators.  The read and write independent calls
 *              compute the group of 1 or more IOCs and further create
 *              derived MPI datatypes when required by the size of the
 *              contiguous read or write requests.
 *
 *              IOC(0) contains the logical data storage for file offset
 *              zero and all offsets that reside within modulo range of
 *              the subfiling stripe_size.
 *
 *              We cycle through all 'n_io_conentrators' and send a
 *              descriptor to each IOC that has a non-zero sized IO
 *              request to fullfill.
 *
 *              Sending descriptors to an IOC usually gets an ACK or
 *              NACK in response.  For the read operations, we post
 *              asynch READs to receive the file data and wait until
 *              all pending operations have completed.
 *
 * Return:      Success (0) or Faiure (non-zero)
 * Errors:      If MPI operations fail for some reason.
 *
 * Programmer:  Richard Warren
 *              7/17/2020
 *
 * Changes:     Initial Version/None.
 *-------------------------------------------------------------------------
 */
static int
read__independent_async(int n_io_concentrators, hid_t context_id, int64_t offset,
						int64_t elements, int dtype_extent, void *data, io_req_t **io_req)
{
    int          sf_world_size, sf_world_rank;
    int          status = 0, errors = 0;
    int64_t      stripe_size, ioc_row, start_id, ioc_start, ioc_offset, offset_in_stripe;
    int *        io_concentrator = NULL;
    int *        subfile_fd = NULL;
    io_req_t    *sf_io_request = NULL;
    useconds_t   delay = 10;
    int64_t      msg[3] = {0,};

    subfiling_context_t *sf_context = get__subfiling_object(context_id);
    assert(sf_context != NULL);

    /* Calculate the IOC that we'll send the IO request to */
    stripe_size = sf_context->sf_stripe_size;

    start_id         = offset / stripe_size;
    offset_in_stripe = offset % stripe_size;
    ioc_row          = offset / sf_context->sf_blocksize_per_stripe;

    ioc_start   = start_id % n_io_concentrators;
    ioc_offset  = (ioc_row * stripe_size) + offset_in_stripe;

    io_concentrator = sf_context->topology->io_concentrator;
    assert(io_concentrator != NULL);

    /* Make sure that we can return a request structure
     * if everything is working correctly
     */
    assert(io_req);

    sf_world_size = sf_context->topology->app_layout->world_size;
    sf_world_rank = sf_context->topology->app_layout->world_rank;

    /* Prepare an IO request.
     * This gets sent to the ioc identified by the file offset
     */
    msg[0] = elements;
    msg[1] = ioc_offset;
    msg[2] = context_id;

#if 0
    printf("[%d %s] offset=%ld, ioc=%d, elements=%ld, ioc_offset=%ld\n",
		   sf_world_rank, __func__, offset, (int)ioc_start, elements, ioc_offset );
	fflush(stdout);
#endif
    status = MPI_Send(msg, 3, MPI_INT64_T, io_concentrator[ioc_start],
            READ_INDEP, sf_context->sf_msg_comm);

    if (status != MPI_SUCCESS) {
        int  len;
        char estring[MPI_MAX_ERROR_STRING];
        MPI_Error_string(status, estring, &len);
        printf("[%d] ERROR! MPI_Send request header (%ld) "
                   "bytes to %d returned an error(%s)\n",
               sf_world_rank, sizeof(msg), io_concentrator[ioc_start], estring);
        fflush(stdout);
        return -1;
    }

    /* At this point in the new implementation, we should queue
     * the async recv so that when the top level VFD tells us
     * to complete all pending IO requests, we have all the info
     * we need to accomplish that.
     */
    sf_io_request = (io_req_t *)malloc(sizeof(io_req_t));
    assert(sf_io_request);

    sf_io_request->completion_func.io_args.ioc        = (int)ioc_start;
    sf_io_request->completion_func.io_args.context_id = context_id;
    sf_io_request->completion_func.io_args.offset     = offset;
    sf_io_request->completion_func.io_args.elements   = elements;
    sf_io_request->completion_func.io_args.data       = data;
    sf_io_request->completion_func.io_args.io_req     = MPI_REQUEST_NULL;
    sf_io_request->completion_func.io_function        = async_completion;
    sf_io_request->completion_func.pending            = 0;

    sf_io_request->prev = sf_io_request->next = NULL;
    /* Start the actual data transfer */

    status = MPI_Irecv(data, (int)elements,
                       MPI_BYTE, io_concentrator[ioc_start],
                       READ_INDEP_DATA,
                       sf_context->sf_data_comm,
                       &sf_io_request->completion_func.io_args.io_req);

    if (status == MPI_SUCCESS) {
        sf_io_request->completion_func.pending = 1;
		*io_req = sf_io_request;
    }
	else {
        puts("MPI_Irecv must have failed!");
        free(sf_io_request);
		*io_req = NULL;
	}

    return status;
} /* end read__independent_async() */

#ifdef USE_ORIGINAL_CODE
static int
read__independent(int n_io_concentrators, hid_t context_id, int64_t offset,
    int64_t elements, int dtype_extent, void *data)
{
    int          i, ioc, active_sends=0 ,n_waiting = 0, errors = 0, status = 0;
    int          sf_world_rank;
    int *        io_concentrator = NULL;
    int *        subfile_fd = NULL;

    int          indices[n_io_concentrators];
    MPI_Request  reqs[n_io_concentrators];
    MPI_Status   stats[n_io_concentrators];
    int64_t      sourceOffset;
    int64_t      source_data_offset[n_io_concentrators];
    int64_t      ioc_read_datasize[n_io_concentrators];
    int64_t      ioc_read_offset[n_io_concentrators];
    MPI_Datatype ioc_read_type[n_io_concentrators];
    char *       sourceData = (char *) data;

    useconds_t   delay = 10;

    subfiling_context_t *sf_context = get__subfiling_object(context_id);
    assert(sf_context != NULL);
    sf_world_rank = sf_context->topology->app_layout->world_rank;
    subfile_fd = sf_context->topology->subfile_fd;
    assert(subfile_fd != NULL);

    io_concentrator = sf_context->topology->io_concentrator;

    /* Prepare the IOCs with a message which indicates the length
     * and file offset for the actual data to be provided.
     */
    for (ioc = 0; ioc < n_io_concentrators; ioc++) {
        int64_t msg[3] = {ioc_read_datasize[ioc], ioc_read_offset[ioc], sf_context->sf_context_id};
        int     packsize = 0;

        /* We may not require data from this IOC...
         * or we may read the data directly from the file!
         * Check the size to verify!
         */

        reqs[ioc] = MPI_REQUEST_NULL;

        if (ioc_read_datasize[ioc] == 0) {
            continue;
        }

        sourceOffset = source_data_offset[ioc];

        if (sf_enable_directIO && (subfile_fd[ioc] == 0)) {
            char *subfile_path = get_ioc_subfile_path(ioc, n_io_concentrators, sf_context);
            if (subfile_path != NULL) {
                subfile_fd[ioc] = open(subfile_path, O_RDONLY);
                if (subfile_fd[ioc] < 0) {
                    perror("opening subfile failed\n");
                }
                else {
                    if (ioc_read_type[ioc] == MPI_BYTE) {
                        errors += sf_read_data(subfile_fd[ioc], ioc_read_offset[ioc], &sourceData[sourceOffset], ioc_read_datasize[ioc], ioc);
                    }
                    else {
                        char *contig_buffer = NULL;
                        int position = 0;
                        int datasize = (int)ioc_read_datasize[ioc];
                        contig_buffer = (char *)HDmalloc((size_t)datasize);
                        assert(contig_buffer != NULL);
                        errors += sf_read_data(subfile_fd[ioc], ioc_read_offset[ioc], contig_buffer, ioc_read_datasize[ioc], ioc);
                        if (MPI_Unpack(contig_buffer, datasize, &position, &sourceData[sourceOffset],
                                       1, ioc_read_type[ioc], MPI_COMM_SELF) != MPI_SUCCESS) {
                            puts("MPI_Unpack error!");
                        }
                        if (contig_buffer)
                            HDfree(contig_buffer);
                    }

                    if (close(subfile_fd[ioc]) < 0) {
                        perror("closing subfile failed");
                    }

                    subfile_fd[ioc] = 0;

                    if (errors) {
                        puts("sf_read_data returned an error!");
                    }
                    else continue;
                }
            }
        }

        active_sends++;
#ifndef NDEBUG
        if (sf_verbose_flag && (client_log != NULL)) {
            fprintf(client_log,
                    "[%d %s]: read_src[ioc(%d), "
                    "sourceOffset=%ld, datasize=%ld, foffset=%ld]\n",
                    sf_world_rank, __func__, ioc, sourceOffset,
                    ioc_read_datasize[ioc], ioc_read_offset[ioc]);
            fflush(client_log);
        }
#endif

        status = MPI_Send(msg, 3, MPI_INT64_T, io_concentrator[ioc],
                           READ_INDEP, sf_context->sf_msg_comm);
        if (status != MPI_SUCCESS) {
            printf("[%d] MPI_Send failure!", sf_world_rank);
            return status;
        } else {
            /* FIXME: We should be checking that we don't exceed
             * the 2GB size limit. This could probably be part
             * of the vector construction rather than validating here.
             */
            if (ioc_read_type[ioc] == MPI_BYTE) {
                int bytes = (int) ioc_read_datasize[ioc];
                status = MPI_Irecv(&sourceData[sourceOffset], bytes,
                    ioc_read_type[ioc], io_concentrator[ioc], READ_INDEP_DATA,
                    sf_context->sf_data_comm, &reqs[ioc]);
            } else {
                MPI_Pack_size(1, ioc_read_type[ioc], MPI_COMM_WORLD, &packsize);
                status = MPI_Irecv(&sourceData[sourceOffset], 1,
                    ioc_read_type[ioc], io_concentrator[ioc], READ_INDEP_DATA,
                    sf_context->sf_data_comm, &reqs[ioc]);
            }
            if (status != MPI_SUCCESS) {
                int  length = 256;
                char error_string[length];
                MPI_Error_string(status, error_string, &length);
                printf("(%s) MPI_Irecv error: %s\n", __func__, error_string);
                return status;
            }
            n_waiting++;
        }
    }

    if (active_sends > 1) {
#ifndef NDEBUG
       if (sf_verbose_flag && (client_log != NULL)) {
          fprintf(client_log,
                  "[%d %s] ioc spillovers = %d\n",
                  sf_world_rank, __func__, active_sends-1);
       }
#endif
    }
    /* We've queued all of the Async READs, now we just need to
     * complete them in any order...
     */
    while (n_waiting) {
        int ready = 0;
        status = MPI_Testsome(n_io_concentrators, reqs, &ready, indices, stats);
        if (status != MPI_SUCCESS) {
            int  length = 256;
            char error_string[length];
            MPI_Error_string(status, error_string, &length);
            printf("(%s) MPI_Waitsome error: %s\n", __func__, error_string);
            for (i = 0; i < n_waiting; i++) {
                printf(
                    "stats[%d].SOURCE=%d, stats.TAG=%d, stats.MPI_ERROR=%d\n",
                    i, stats[i].MPI_SOURCE, stats[i].MPI_TAG,
                    stats[i].MPI_ERROR);
                fflush(stdout);
            }
            return status;
        }

        for (i = 0; i < ready; i++) {
#ifndef NDEBUG
            if (sf_verbose_flag && (client_log != NULL)) {
                fprintf(client_log,
                        "[%d] READ bytes(%ld) of data from ioc_concentrator %d "
                        "complete\n",
                        sf_world_rank, ioc_read_datasize[indices[i]],
                        indices[i]);
                fflush(client_log);
            }
#endif
            if (ioc_read_type[indices[i]] != MPI_BYTE) {
                MPI_Type_free(&ioc_read_type[indices[i]]);
            }
            n_waiting--;
        }

        if (n_waiting) {
            usleep(delay);
        }
    }
    return status;
}
#endif

/*-------------------------------------------------------------------------
 * Function:    get_ioc_subfile_path
 *
 * Purpose:     We provide a utility function to generate a subfiling
 *              filename from a template.  While the user provides a
 *              name which will serve as the HDF5 file name, sub-filing
 *              files are related to the user filename via the filesytem 
 *              inode identifier.  The inode id can be utilized as a
 *              global unique identifier (GUID) which provides a
 *              grouping ID to easily distinguish subfiles.
 *
 *              The inode_id is contained in the 'sf_context' structure.
 *
 * Return:      A full filepath which should be copied, e.g. using strdup
 *-------------------------------------------------------------------------
 */
static char *
get_ioc_subfile_path(int ioc, int ioc_count, subfiling_context_t *sf_context)
{
    static char filepath[PATH_MAX];
    char *subfile_dir = NULL;
    char *prefix = sf_context->subfile_prefix;

    if (prefix != NULL) {
        sprintf(filepath, "%s/" SF_FILENAME_TEMPLATE, prefix,
                sf_context->h5_file_id, ioc, ioc_count);
    }
    else {
        strcpy(filepath,sf_context->filename);
        subfile_dir = strrchr(filepath, '/');
        assert(subfile_dir);
        sprintf(subfile_dir+1, SF_FILENAME_TEMPLATE, sf_context->h5_file_id,
                ioc, ioc_count);
    }
    return filepath;
} /* end get_ioc_subfile_path() */



/*-------------------------------------------------------------------------
 * Utility functions in support of a first pass attempt at handling async
 * IO.  The working assumption is that reads and writes to a collection
 * of IO Concentrators (IOCs) will proceed by stages.  In the first stage,
 * each MPI rank will get their individual IOs started by preping the IOC
 * with a message which indicates (via the MPI tag) what operation is
 * starting, along with the file offset, data size, and a context_id.
 * The latter will be used to access the actual open file descriptor.
 * 
 *-------------------------------------------------------------------------
 * Function:    progress_this_pending_io
 *
 * Purpose:     In this initial example, we can progress an individual
 *              IO request which is described by the io_req_t input arg.
 *
 * Return:      an integer status.  Zero(0) indicates success. Negative
 *              values (-1) indicates an error.
 *-------------------------------------------------------------------------
 */
static int
progress_this_pending_io(io_req_t *this_req)
{
    assert(this_req);
    assert(this_req->completion_func.io_function);
    return (*this_req->completion_func.io_function)(&this_req->completion_func);
}

/*-------------------------------------------------------------------------
 * Function:    write_data
 *
 * Purpose:     Given a io_func_t structure containing the function pointer
 *              and it's input arguments, we write the supplied data out
 *              asynchronous using MPI_Isend, to the appropriate IOC.
 *
 * Return:      an integer status.  Zero(0) indicates success. Negative
 *              values (-1) indicates an error.
 *-------------------------------------------------------------------------
 */
static int
write_data(io_func_t *this_func)
{
    int ioc, status;
    hid_t context_id;
    int64_t offset, elements;
    void *data;
    MPI_Request *io_req = NULL;
    int *io_concentrator = NULL;
    subfiling_context_t *sf_context = NULL;
    assert(this_func);

    sf_context = get__subfiling_object(this_func->io_args.context_id);
    // printf("%s: context_id = %ld\n", __func__, this_func->io_args.context_id);

    assert(sf_context);

    io_concentrator = sf_context->topology->io_concentrator;
    ioc = this_func->io_args.ioc;

    status = MPI_Isend(data, (int)elements,
                       MPI_BYTE, io_concentrator[ioc], WRITE_INDEP_DATA,
                       sf_context->sf_data_comm, &this_func->io_args.io_req);
    // printf("%s - Isend: io_req = %u\n",__func__, this_func->io_args.io_req);
	// fflush(stdout);
    return status;
}

/*-------------------------------------------------------------------------
 * Function:    async_completion
 *
 * Purpose:     Given a single io_func_t structure containing the function
 *              pointer and it's input arguments and a single MPI_Request
 *              argument which needs to be completed, we make progress
 *              by calling MPI_Test.  In this initial example, we loop
 *              until the request is completed as indicated by a non-zero
 *              flag variable.
 *
 *              As we go further with the implemention, we anticipate that
 *              rather than testing a single request variable, we will
 *              deal with a collection of all pending IO requests (on
 *              this rank).
 *
 * Return:      an integer status.  Zero(0) indicates success. Negative
 *              values (-1) indicates an error.
 *-------------------------------------------------------------------------
 */
static int
async_completion(void *arg)
{
    struct async_arg {
		int n_reqs;
		MPI_Request *sf_reqs;
	} *in_progress = (struct async_arg *)arg;

    assert(arg);
    int status, errors = 0;
    int count = in_progress->n_reqs;
    int n_waiting = count;
    int indices[count];
	MPI_Status stats[count];
    useconds_t delay = 5;

    while (n_waiting) {
        int i, ready = 0;
        status = MPI_Testsome(count, in_progress->sf_reqs, &ready, indices, stats);
        if (status != MPI_SUCCESS) {
            int  len;
            char estring[MPI_MAX_ERROR_STRING];
            MPI_Error_string(status, estring, &len);
            printf("[%s] MPI_ERROR! MPI_Testsome returned an error(%s)\n",
                __func__, estring);
            fflush(stdout);
            errors++;
			return -1;
        }

        if (ready == 0) {
            usleep(delay);
        }

        for (i = 0; i < ready; i++) {
            n_waiting--;
        }
    }
    return errors;	
}


/*-------------------------------------------------------------------------
 * Function:    Internal write__independent_async.
 *
 * Purpose:     The IO operations can be striped across a selection of
 *              IO concentrators.  The read and write independent calls
 *              compute the group of 1 or more IOCs and further create
 *              derived MPI datatypes when required by the size of the
 *              contiguous read or write requests.
 *
 *              IOC(0) contains the logical data storage for file offset
 *              zero and all offsets that reside within modulo range of
 *              the subfiling stripe_size.
 *
 *              We cycle through all 'n_io_conentrators' and send a
 *              descriptor to each IOC that has a non-zero sized IO
 *              request to fullfill.
 *
 *              Sending descriptors to an IOC usually gets an ACK or
 *              NACK in response.  For the write operations, we post
 *              asynch READs to receive ACKs from IOC ranks that have
 *              allocated memory receive the data to write to the
 *              subfile.  Upon receiving an ACK, we send the actual
 *              user data to the IOC.
 *
 * Return:      Success (0) or Faiure (non-zero)
 * Errors:      If MPI operations fail for some reason.
 *
 * Programmer:  Richard Warren
 *              7/17/2020
 *
 * Changes:     Initial Version/None.
 *-------------------------------------------------------------------------
 */
static int
write__independent_async(int n_io_concentrators, hid_t context_id, int64_t offset,
						 int64_t elements, int dtype_extent, void *data, io_req_t **io_req)
{
    int          sf_world_size, sf_world_rank;
    int          ack = 0, active_sends = 0, n_waiting = 0, status = 0, errors = 0;
    int64_t      stripe_size, ioc_row, start_id, ioc_start, ioc_offset, offset_in_stripe;
    int *        io_concentrator = NULL;
    int *        subfile_fd = NULL;
    io_req_t    *sf_io_request = NULL;
    MPI_Request  ackrequest;
    int64_t      msg[3] = {0,};

    subfiling_context_t *sf_context = get__subfiling_object(context_id);
    assert(sf_context != NULL);

    /* Calculate the IOC that we'll send the IO request to */
    stripe_size = sf_context->sf_stripe_size;

    start_id         = offset / stripe_size;
    offset_in_stripe = offset % stripe_size;
    ioc_row          = offset / sf_context->sf_blocksize_per_stripe;

    ioc_start   = start_id % n_io_concentrators;
    ioc_offset  = (ioc_row * stripe_size) + offset_in_stripe;

    io_concentrator = sf_context->topology->io_concentrator;
    assert(io_concentrator != NULL);

    /* Make sure that we can return a request structure
     * if everything is working correctly
     */
    assert(io_req);

    sf_world_size = sf_context->topology->app_layout->world_size;
    sf_world_rank = sf_context->topology->app_layout->world_rank;

    /* Prepare an IO request.
     * This gets sent to the ioc identified by the file offset.
     * (see above: Calculate the IOC))
     */
    msg[0] = elements;
    msg[1] = ioc_offset;
    msg[2] = context_id;

#if 0
    printf("[%d %s] offset=%ld, ioc=%d, elements=%ld, ioc_offset=%ld\n",
		   sf_world_rank, __func__, offset, (int)ioc_start, elements, ioc_offset );
	fflush(stdout);
#endif
    status = MPI_Send(msg, 3, MPI_INT64_T, io_concentrator[ioc_start],
            WRITE_INDEP, sf_context->sf_msg_comm);

    if (status != MPI_SUCCESS) {
        int  len;
        char estring[MPI_MAX_ERROR_STRING];
        MPI_Error_string(status, estring, &len);
        printf("[%d] ERROR! MPI_Send of %ld bytes to %d returned an "
                   "error(%s)\n",
               sf_world_rank, sizeof(msg), io_concentrator[ioc_start], estring);
        fflush(stdout);
        return -1;
    }
    else active_sends++;
    /* 
     * We wait for memory to be allocated on the target IOC so that we can
     * start sending user data. Once memory is allocated, we will receive
     * an ACK (or NACK) message from the IOC to allow us to proceed.
     */
    // printf("Posting Irec for an ACK\n");
    // fflush(stdout);
    status = MPI_Irecv(&ack, 1, MPI_INT, io_concentrator[ioc_start],
                       WRITE_INDEP_ACK, sf_context->sf_data_comm, &ackrequest);

    if (status != MPI_SUCCESS) {
        printf("[%d %s] MPI_Irecv failed\n", sf_world_rank, __func__);
        fflush(stdout);
        return -1;
    }
    
    n_waiting = active_sends;
    // printf("Waiting for an ACK\n");
    // fflush(stdout);

    while( n_waiting ) {
        int flag = 0;
        status = MPI_Test(&ackrequest, &flag, MPI_STATUS_IGNORE);
        if (status == MPI_SUCCESS) {
            if (flag == 0)
                usleep(0);
            else { 
                n_waiting--;
                if (ack == 0) {	/* NACK */
                    printf("%s - Received NACK!\n", __func__);
                }
#if 0
                else {
                    puts("Received ACK");
                }
#endif
            }
        }
    }

    /* At this point in the new implementation, we should queue
     * the async write so that when the top level VFD tells us
     * to complete all pending IO requests, we have all the info
     * we need to accomplish that.
     */
    sf_io_request = (io_req_t *)malloc(sizeof(io_req_t));
    assert(sf_io_request);

    sf_io_request->completion_func.io_args.ioc        = (int)ioc_start;
    sf_io_request->completion_func.io_args.context_id = context_id;
    sf_io_request->completion_func.io_args.offset     = offset;
    sf_io_request->completion_func.io_args.elements   = elements;
    sf_io_request->completion_func.io_args.data       = data;
    sf_io_request->completion_func.io_args.io_req     = MPI_REQUEST_NULL;
    sf_io_request->completion_func.io_function        = async_completion;
    sf_io_request->completion_func.pending            = 0;

    sf_io_request->prev = sf_io_request->next = NULL;
    /* Start the actual data transfer */

    status = MPI_Isend(data, (int)elements,
                       MPI_BYTE, io_concentrator[ioc_start],
                       WRITE_INDEP_DATA,
                       sf_context->sf_data_comm,
                       &sf_io_request->completion_func.io_args.io_req);

	// printf("isend of %ld bytes to ioc(%d) returned %d, io_req=0x%x\n", elements,
    //        (int)ioc_start, status, sf_io_request->completion_func.io_args.io_req);
    // fflush(stdout);

    /* When we actually have the async IO support,
     * the request should be queued before we 
     * return to the caller. 
     * Having queued the IO operation, we might want to
     * get additional work started before allowing the
     * queued IO requests to make further progress and/or
     * to complete, so we just return to the caller.
     */

    if (status == MPI_SUCCESS) {
        sf_io_request->completion_func.pending = 1;
		*io_req = sf_io_request;
    }
	else {
        puts("MPI_Isend must have failed!");
        free(sf_io_request);
		*io_req = NULL;
	}
    return status;
} /* end write__independent_async() */


#ifdef USE_ORIGINAL_CODE
static int
write__independent(int n_io_concentrators, hid_t context_id, int64_t offset,
    int64_t elements, int dtype_extent, void *data)
{
    int          sf_world_size, sf_world_rank;
    int *        io_concentrator = NULL;
    int *        subfile_fd = NULL;

    int          acks[n_io_concentrators];
    int          indices[n_io_concentrators];
    MPI_Request  reqs[n_io_concentrators];
    MPI_Request  ackreq[n_io_concentrators];
    MPI_Status   stats[n_io_concentrators];
    int64_t      sourceOffset;
    int64_t      source_data_offset[n_io_concentrators];
    int64_t      ioc_write_datasize[n_io_concentrators];
    int64_t      ioc_write_offset[n_io_concentrators];
    MPI_Datatype ioc_write_type[n_io_concentrators];
    int          active_sends = 0, n_waiting = 0, status = 0, errors = 0;
    int          direct_writes = 0, skipped_writes=0;
    int          i, target, ioc;
    double       t0, t1, t2;
    char *       sourceData = (char *) data;
    useconds_t   delay = 20;


    subfiling_context_t *sf_context = get__subfiling_object(context_id);
    assert(sf_context != NULL);

    io_concentrator = sf_context->topology->io_concentrator;
    assert(io_concentrator != NULL);

    subfile_fd = sf_context->topology->subfile_fd;
    assert(subfile_fd != NULL);

    sf_world_size = sf_context->topology->app_layout->world_size;
    sf_world_rank = sf_context->topology->app_layout->world_rank;

    if (sf_context->topology->rank_is_ioc) {
        sf_context->sf_write_count++;
    }

    /* The following function will initialize the collection of IO transfer
     * parameters, i.e. local memory (source) offsets, target file offsets,
     * target data sizes (in bytes), and a MPI Datatype for each of the
     * IO concentrator transactions.
     *
     * For small transfers, at least 1 IOC instance will have valid info.
     * For larger transfers, it is likely that the full set of
     * n_io_concentrators will be utilized.  If the total transaction size is
     * less than n_io_concentrators X stripe_size, then the MPI datatype should
     * probably be MPI_BYTE.  Larger tranactions will create MPI derived
     * datatypes to span the entire logical collection of stripes.  Said
     * differently, the largest IO requests will require a stripe depth greater
     * than one.
     */
#if 0
    if (init__indep_io(sf_context, source_data_offset, ioc_write_datasize,
            ioc_write_offset, ioc_write_type, offset, elements,
            dtype_extent) < 0) {
        return -1;
    }
#endif
    /* 
     * Prepare the IOCs with a message which indicates the length
     * of the actual data to be written.  We also provide the file
     * offset so that when the IOC recieves the data (in whatever order)
     * they can lseek to the correct offset and write the data.
     * These messages are all posted with MPI_Isend, so that we can
     * have overlapping operations running in parallel...
     */
    t0 = MPI_Wtime();           /* Clock snapshot */

    for (target = 0; target < n_io_concentrators; target++) {
        int64_t msg[3] = { 0, };
        ioc = (sf_world_rank + target) % n_io_concentrators;
        sourceOffset = source_data_offset[ioc];
        // printf("preparing for ioc(%d) write: datasize = %ld\n", ioc, ioc_write_datasize[ioc]);
        msg[0] = ioc_write_datasize[ioc];
        msg[1] = ioc_write_offset[ioc];
        msg[2] = sf_context->sf_context_id;
        acks[ioc] = 0;
        reqs[ioc] = MPI_REQUEST_NULL;
        ackreq[ioc] = MPI_REQUEST_NULL;
        /* Having refactor'ed things, we aren't currently
         * handling anything other than MPI_BYTE datatypes.
         * The refactoring changes avoids using derived types
         * by breaking larger data transfers into multiple
         * segments.  As a result, for larger IO operations
         * each IOC may require multiple IO requests to satisfy
         * a single user IO.
         * 
         * Here, we presume that each transaction with an IOC
         * thread will consist of a data IO which does not
         * require a derived datatype, so simply assign
         * the operational type as MPI_BYTE.
         */
        ioc_write_type[ioc] = MPI_BYTE;

        if (ioc_write_datasize[ioc] == 0) {
            skipped_writes++;
            continue;
        }

        sourceOffset = source_data_offset[ioc];

        if (sf_enable_directIO && (subfile_fd[ioc] == 0)) {
            char *subfile_path = get_ioc_subfile_path(ioc, n_io_concentrators, sf_context);
            if (subfile_path != NULL) {
                subfile_fd[ioc] = open(subfile_path, O_RDWR /* | O_DSYNC */ );
                if (subfile_fd[ioc] < 0) {
                    perror("opening subfile failed\n");
                }
                else {
                    if (ioc_write_type[ioc] == MPI_BYTE) {
                        // printf("[%d] Write ioc(%d): MPI_BYTE elements = %lld, sourceOffset = %lld, foffset = %lld\n",
                        // sf_world_rank, ioc, elements, sourceOffset, ioc_write_offset[ioc]);
                        errors += sf_write_data(subfile_fd[ioc], ioc_write_offset[ioc],
                                                &sourceData[sourceOffset], ioc_write_datasize[ioc], ioc);
                    }
                    else {      /* This data is larger than a single stripe */
                        char *contig_buffer = NULL;
                        int type_size = 0;
                        int position = 0;
                        int outsize = (int)ioc_write_datasize[ioc];
                        contig_buffer = (char *)HDmalloc((size_t)outsize);
                        assert(contig_buffer != NULL);
                        MPI_Type_size(ioc_write_type[ioc], &type_size);
                        // printf("[%d] Write ioc(%d): type_size = %d, outsize = %d elements = %lld\n",
                        // sf_world_rank, ioc, type_size, outsize, elements);
                        if (MPI_Pack(data, 1, ioc_write_type[ioc], contig_buffer,
                            outsize, &position, MPI_COMM_SELF) == MPI_SUCCESS) {
                            ioc_write_type[ioc] = (int64_t) outsize;
                            errors += sf_write_data(subfile_fd[ioc], ioc_write_offset[ioc],
                                      contig_buffer, ioc_write_datasize[ioc], ioc);
                        }
                        else {
                            puts("MPI_Pack failure!");
                            errors++;
                        }
                        if (contig_buffer)
                            HDfree(contig_buffer);
                    }

                    if (close(subfile_fd[ioc]) < 0) {
                        perror("closing subfile failed");
                    }
                    subfile_fd[ioc] = 0;

                    if (errors == 0) {
                        direct_writes++;
                        continue;
                    }
                }
            }
        }
#ifndef NDEBUG
        if (sf_verbose_flag && (client_log != NULL)) {
            fprintf(client_log,
                    "[%d %s]: write_dest[ioc(%d), "
                    "sourceOffset=%ld, datasize=%ld, foffset=%ld]\n",
                    sf_world_rank, __func__, ioc, sourceOffset,
                    ioc_write_datasize[ioc], ioc_write_offset[ioc]);
            fflush(client_log);
        }
#endif

        active_sends++;

        /*
         * Send the Message HEADER which indicates the requested IO operation
         * (via the message TAG) along with the data size and file offset.
         */

        status = MPI_Send(msg, 3, MPI_INT64_T, io_concentrator[ioc],
            WRITE_INDEP, sf_context->sf_msg_comm);

        if (status != MPI_SUCCESS) {
            int  len;
            char estring[MPI_MAX_ERROR_STRING];
            MPI_Error_string(status, estring, &len);
            printf("[%d] ERROR! MPI_Send of %ld bytes to %d returned an "
                   "error(%s)\n",
                sf_world_rank, sizeof(msg), io_concentrator[ioc], estring);
            fflush(stdout);
            break; /* If unable to send to an IOC, we can call it quits...  */
        }

        /* 
         * We wait for memory to be allocated on the target IOC so that we can
         * start sending user data. Once memory is allocated, we will receive
         * an ACK (or NACK) message from the IOC to allow us to proceed.
         */
        status = MPI_Irecv(&acks[ioc], 1, MPI_INT, io_concentrator[ioc],
                           WRITE_INDEP_ACK, sf_context->sf_data_comm, &ackreq[ioc]);

        if (status != MPI_SUCCESS) {
            printf("[%d %s] MPI_Irecv failed\n", sf_world_rank, __func__);
            fflush(stdout);
        }
    }

    n_waiting = active_sends;

    while (n_waiting) {
        int ready = 0;
        status = MPI_Testsome(n_io_concentrators, ackreq, &ready, indices, stats);
        if (status != MPI_SUCCESS) {
            int  len = MPI_MAX_ERROR_STRING;
            char estring[MPI_MAX_ERROR_STRING];
            MPI_Error_string(status, estring, &len);
            printf("[%d %s] MPI_ERROR! MPI_Testsome returned an error(%s)\n",
                sf_world_rank, __func__, estring);
            fflush(stdout);
            errors++;
        }

        if (ready == 0) {
            usleep(delay);
        }

        for (i=0; i < ready; i++) {
            int from = indices[i];
            sourceOffset = source_data_offset[from];
            if (ackreq[from] > 0) {
                if (ioc_write_type[from] == MPI_BYTE) {
                    int datasize = (int) (ioc_write_datasize[from] & INT32_MASK);
                    status = MPI_Isend(&sourceData[sourceOffset], datasize,
                        MPI_BYTE, io_concentrator[from], WRITE_INDEP_DATA,
                        sf_context->sf_data_comm, &reqs[from]);
                } else {
                    status = MPI_Isend(&sourceData[sourceOffset], 1,
                        ioc_write_type[from], io_concentrator[from],
                        WRITE_INDEP_DATA, sf_context->sf_data_comm, &reqs[from]);
                }
            }
            else {
                errors++;
                puts("ACK error!");
                fflush(stdout);
            }
#ifndef NDEBUG
            if (sf_verbose_flag && (client_log != NULL)) {
                fprintf(client_log, "[%d] received ack(%d) from ioc(%d)\n",
                        sf_world_rank, acks[from], from);
                fflush(client_log);
             }
#endif
            /* Check the status of our MPI_Isend... */
            if (status != MPI_SUCCESS) {
                int  len = MPI_MAX_ERROR_STRING;
                char estring[MPI_MAX_ERROR_STRING];
                MPI_Error_string(status, estring, &len);
                printf("[%d %s] MPI_ERROR! MPI_Isend returned an error(%s)\n",
                       sf_world_rank, __func__, estring);
                fflush(stdout);
                errors++;
            }
            n_waiting--;
        }
    }

#ifndef NDEBUG
    t1 = MPI_Wtime();
    if (sf_verbose_flag && (client_log != NULL)) {
        fprintf(client_log, "[%d %s] time awaiting ACKs: %lf seconds)\n",
                sf_world_rank, __func__, (t1 - t0));
    }
#endif


    /* Reset n_waiting from waiting on ACKS to waiting on the data Isends */
    n_waiting = active_sends;

    while (n_waiting) {
        int ready = 0;
        status = MPI_Testsome(n_io_concentrators, reqs, &ready, indices, stats);
        if (status != MPI_SUCCESS) {
            int  len;
            char estring[MPI_MAX_ERROR_STRING];
            MPI_Error_string(status, estring, &len);
            printf("[%d %s] MPI_ERROR! MPI_Waitsome returned an error(%s)\n",
                sf_world_rank, __func__, estring);
            fflush(stdout);
            errors++;
        }

        if (ready == 0) {
            usleep(delay);
        }

        for (i = 0; i < ready; i++) {
            /* One of the Issend calls has completed
             * If we used a derived type to send data, then should free
             * that datatype instance.
             */
            if (ioc_write_type[indices[i]] != MPI_BYTE) {
                MPI_Type_free(&ioc_write_type[indices[i]]);
            }
            n_waiting--;
        }
    }


#ifndef NDEBUG
    t2 = MPI_Wtime();
    if (sf_verbose_flag && (client_log != NULL)) {
        if (active_sends > 0) {
            fprintf(client_log, "[%d %s] sending data = %lf seconds of total_time = %lf\n",
                sf_world_rank, __func__, (t2 - t1), (t2 - t0));
        }
        else {
            fprintf(client_log, "[%d %s] direct_writes(%d) of data = %lf seconds of total_time = %lf; skipped_writes(%d)\n",
                    sf_world_rank, __func__, direct_writes, (t2 - t1), (t2 - t0), skipped_writes);
        }
    }
#endif

    if (active_sends > 1) {
#ifndef NDEBUG
       if (sf_verbose_flag && (client_log != NULL)) {
          fprintf(client_log,
                  "[%d %s] ioc spillovers = %d\n",
                  sf_world_rank, __func__, active_sends-1);
          fflush(client_log);
       }
#endif
    }

    if (errors)
        return -1;
    return status;
} /* end write__independent() */

/*-------------------------------------------------------------------------
 * Function:    Public/Client sf_write_independent
 *
 * Purpose:     A public function which wraps the Internal version
 *              and allows the addition of the additional 'n_io_concentrator'
 *              argument.  This is important as it allows me to skip
 *              memory allocation functions since storage for the various
 *              vector variables is on the call stack...
 *
 * Return:      The integer status returned by the Internal read_independent
 *              function.  Successful operations will return 0.
 * Errors:      An MPI related error value.
 *
 * Programmer:  Richard Warren
 *              7/17/2020
 *
 * Changes:     Initial Version/None.
 *
 *-------------------------------------------------------------------------
 */
int
sf_write_independent(hid_t sf_fid, int64_t offset, int64_t elements,
    int dtype_extent, const void *data)
{
    hid_t sf_context_id = fid_map_to_context(sf_fid);
    subfiling_context_t *sf_context = get__subfiling_object(sf_context_id);

    assert(sf_context != NULL);
    return write__independent(sf_context->topology->n_io_concentrators,
        sf_context_id, offset, elements, dtype_extent, data);
}

/*-------------------------------------------------------------------------
 * Function:    Public/Client sf_write_vector
 *
 * Purpose:     Another write__independent wrapper.  As with the
 *              sf_read_vector function, we simply loop over the vector
 *              elements and call the underlying write_independent function.
 *
 * Return:      The integer status returned by the Internal read_independent
 *              function.  Successful operations will return 0.
 * Errors:      An MPI related error value.
 *
 * Programmer:  Richard Warren
 *              7/17/2020
 *
 * Changes:     Initial Version/None.
 *
 *-------------------------------------------------------------------------
 */
herr_t
sf_write_vector(hid_t h5_fid, hssize_t count, haddr_t *addrs, hsize_t sizes[],
    void *bufs[] /* data_in */)
{
    hssize_t             status=0, k = 0;
    herr_t               ret_value = SUCCEED;
    hid_t                sf_context_id = fid_map_to_context(h5_fid);
    subfiling_context_t *sf_context = get__subfiling_object(sf_context_id);

    assert(sf_context != NULL);

    /*
     * Call the underlying write function for each vector element.
     * 'count' represents the number of IOCs that we are dealing with.
     * Each loop should be writing data to a new subfile. 
     * 'addrs' represents the filespace address which needs to be decoded
     * into the appropriate IOC and subfile offset.
     */
    for (k = 0; k < count; k++) {
        status = write__independent(sf_context->topology->n_io_concentrators,
                                    sf_context_id, (int64_t) addrs[k],
                                    (int64_t) sizes[k], 1, bufs[k]);
        if (status < 0) {
            printf("%s - encountered an internal error!\n", __func__);
            goto errors;
        }
    }
    return ret_value;

errors:
    return FAIL;
}   /* end sf_write_vector() */

#endif	/* #ifdef USE_ORIGINAL_CODE */

/* 
 * FIXME::
 * Refactored version of the above sf_write_vector() function.
 * The H5FD__ioc_write_vector VFD call includes additional 'hid_t dxpl'
 * and 'H5FD_mem_t types[]'. Those are included here at present, but
 * need to be addressed at some point.
 */
herr_t
H5FD__write_vector_internal(hid_t h5_fid, hssize_t count, haddr_t addrs[], hsize_t sizes[], void *bufs[] /* data_in */)
{
    herr_t               ret_value = SUCCEED;   
    hssize_t             status=0, k = 0;
    hid_t                sf_context_id = fid_map_to_context(h5_fid);
    subfiling_context_t *sf_context = NULL;
    io_req_t            **sf_async_reqs = NULL;
    MPI_Request         *active_reqs = NULL;
    struct __mpi_req {
        int n_reqs;
        MPI_Request *active_reqs;
    } *mpi_reqs = NULL;

    sf_context = get__subfiling_object(sf_context_id);
    assert(sf_context != NULL);

    active_reqs = (MPI_Request *)calloc(count+2, sizeof(MPI_Request));
    assert(active_reqs);

    sf_async_reqs = (io_req_t *)calloc(count, sizeof(io_req_t *));
    assert(sf_async_reqs);

    /* 
     * Note: We allocated extra space in the active_requests (above).
     * The extra should be enough for an integer plus a pointer.
     */
    mpi_reqs = (struct __mpi_req *)&active_reqs[count];
    mpi_reqs->n_reqs = count;    
    mpi_reqs->active_reqs = active_reqs;

    /* Each pass thru the following should queue an MPI write
     * to a new IOC. Both the IOC selection and offset within the
     * particular subfile are based on the combinatation of striping
     * factors and the virtual file offset (addrs[k]).
     */
    for (k=0; k < count; k++) {
        status = write__independent_async(count,
										  sf_context_id, (int64_t) addrs[k],
										  (int64_t) sizes[k], 1, bufs[k],
                                          &sf_async_reqs[k] );
        if (status < 0) {
            printf("%s - encountered an internal error!\n", __func__);
            goto errors;
        }
		else {
            mpi_reqs->active_reqs[k] = sf_async_reqs[k]->completion_func.io_args.io_req;
		}
    }

    /* Here, we should have queued 'count' async requests.
     * We can can now try to complete those before returning
     * to the caller for the next set of IO operations.
     */
    if (sf_async_reqs[0]->completion_func.io_function)
        ret_value = (*sf_async_reqs[0]->completion_func.io_function)(mpi_reqs);

    if (active_reqs) 
        free(active_reqs);

    if (sf_async_reqs) {
        for (k=0; k < count; k++) {
            if (sf_async_reqs[k]) {
                free(sf_async_reqs[k]);
            }
        }
        free(sf_async_reqs);
    }
    return ret_value;

errors:
    return FAIL;

}

/* 
 * FIXME::
 * Refactored version of the above sf_read_vector() function.
 * The H5FD__ioc_read_vector VFD call includes additional 'hid_t dxpl'
 * and 'H5FD_mem_t types[]'. Those are included here at present, but
 * need to be addressed at some point.
 */
herr_t H5FD__read_vector_internal(hid_t h5_fid, hssize_t count, haddr_t addrs[], hsize_t sizes[], void *bufs[] /* data_out */)
{
    herr_t               ret_value = SUCCEED;   
    hssize_t             status=0, k = 0;
    hid_t                sf_context_id = fid_map_to_context(h5_fid);
    subfiling_context_t *sf_context = NULL;
	io_req_t            **sf_async_reqs = NULL;
	MPI_Request         *active_reqs = NULL;
    struct __mpi_req {
		int n_reqs;
		MPI_Request *active_reqs;
	} *mpi_reqs = NULL;

    sf_context = get__subfiling_object(sf_context_id);
    assert(sf_context != NULL);

    active_reqs = (MPI_Request *)calloc(count+2, sizeof(MPI_Request));
    assert(active_reqs);

    sf_async_reqs = (io_req_t *)calloc(count, sizeof(io_req_t *));
    assert(sf_async_reqs);

    /* 
     * Note: We allocated extra space in the active_requests (above).
     * The extra should be enough for an integer plus a pointer.
     */
    mpi_reqs = (struct __mpi_req *)&active_reqs[count];
    mpi_reqs->n_reqs = count;
    mpi_reqs->active_reqs = active_reqs;

    for (k=0; k < count; k++) {
        status = read__independent_async(sf_context->topology->n_io_concentrators,
										  sf_context_id, (int64_t) addrs[k],
										  (int64_t) sizes[k], 1, bufs[k],
                                          &sf_async_reqs[k] );
        if (status < 0) {
            printf("%s - encountered an internal error!\n", __func__);
            goto errors;
        }
		else {
            mpi_reqs->active_reqs[k] = sf_async_reqs[k]->completion_func.io_args.io_req;
		}
    }
    /* Here, we should have queued 'count' async requests 
     * (one to each required IOC).
     *
     * We can can now try to complete those before returning
     * to the caller for the next set of IO operations.
     */
    if (sf_async_reqs[0]->completion_func.io_function)
        ret_value = (*sf_async_reqs[0]->completion_func.io_function)(mpi_reqs);

    if (active_reqs) 
        free(active_reqs);

    if (sf_async_reqs) {
        for (k=0; k < count; k++) {
            if (sf_async_reqs[k]) {
                free(sf_async_reqs[k]);
            }
        }
        free(sf_async_reqs);
    }
    return ret_value;

errors:
    return FAIL;

}


int
sf_truncate(hid_t h5_fid, haddr_t H5_ATTR_PARALLEL_UNUSED addr)
{
    hid_t                sf_context_id = fid_map_to_context(h5_fid);
    subfiling_context_t *sf_context = get__subfiling_object(sf_context_id);

    assert(sf_context != NULL);
    return 0;
}



int
sf_shutdown_local_ioc(hid_t fid)
{
    hid_t context_id = fid_map_to_context(fid);
    subfiling_context_t *sf_context = get__subfiling_object(context_id);
    assert(sf_context != NULL);
    if (sf_context->topology->rank_is_ioc) {
        sf_shutdown_flag = 1;
    }
    return 0;
}

/*-------------------------------------------------------------------------
 * Function:    Public/IOC ioc_main
 *
 * Purpose:     This is the principal function run by the IO Concentrator
 *              main thread.  It remains within a loop until allowed to
 *              exit by means of setting the 'sf_shutdown_flag'.   This
 *              usually accomplished as part of the file close operation.
 *
 *              The function implements an asynchronous polling approach
 *              for incoming messages. These messages can be thought of
 *              as a primitive RPC which utilizes MPI TAGs to code and
 *              implement the desired subfiling functionality.
 *
 *              As each incoming message is received, it get added to
 *              a queue for processing by a thread_pool thread.
 *              The message handlers are dispatched via the
 *              "handle_work_request" ftn (see H5FDsubfile_thread.c)

 *              Subfiling is effectively a software RAID-0 implementation
 *              where having multiple IO Concentrators and independent
 *              subfiles is equated to the multiple disks and a true
 *              hardware base RAID implementation.
 *
 *              IO Concentrators are ordered according to their MPI rank.
 *              In the simplest interpretation, IOC(0) will always contain
 *              the initial bytes of the logical disk image.  Byte 0 of
 *              IOC(1) will contain the byte written to the logical disk
 *              offset "stripe_size" X IOC(number).
 *
 *              Example: If the stripe size is defined to be 256K, then
 *              byte 0 of subfile(1) is at logical offset 262144 of the
 *              file.   Similarly, byte 0 of subfile(2) represents the
 *              logical file offset = 524288.   For logical files larger
 *              than 'N' X stripe_size, we simply "wrap around" back to
 *              subfile(0).  The following shows the mapping of 30
 *              logical blocks of data over 3 subfiles:
 *              +--------+--------+--------+--------+--------+--------+
 *              | blk(0 )| blk(1) | blk(2 )| blk(3 )| blk(4 )| blk(5 )|
 *              | IOC(0) | IOC(1) | IOC(2) | IOC(0) | IOC(1) | IOC(2) |
 *              +--------+--------+--------+--------+--------+--------+
 *              | blk(6 )| blk(7) | blk(8 )| blk(9 )| blk(10)| blk(11)|
 *              | IOC(0) | IOC(1) | IOC(2) | IOC(0) | IOC(1) | IOC(2) |
 *              +--------+--------+--------+--------+--------+--------+
 *              | blk(12)| blk(13)| blk(14)| blk(15)| blk(16)| blk(17)|
 *              | IOC(0) | IOC(1) | IOC(2) | IOC(0) | IOC(1) | IOC(2) |
 *              +--------+--------+--------+--------+--------+--------+
 *              | blk(18)| blk(19)| blk(20)| blk(21)| blk(22)| blk(23)|
 *              | IOC(0) | IOC(1) | IOC(2) | IOC(0) | IOC(1) | IOC(2) |
 *              +--------+--------+--------+--------+--------+--------+
 *              | blk(24)| blk(25)| blk(26)| blk(27)| blk(28)| blk(29)|
 *              | IOC(0) | IOC(1) | IOC(2) | IOC(0) | IOC(1) | IOC(2) |
 *              +--------+--------+--------+--------+--------+--------+
 *
 * Return:      None
 * Errors:      None
 *
 * Programmer:  Richard Warren
 *              7/17/2020
 *
 * Changes:     Initial Version/None.
 *-------------------------------------------------------------------------
 */
int
ioc_main(int64_t context_id)
{
    int                  sf_world_rank, sf_world_size;
    int                  subfile_rank;
    int                  flag, ret;
    int                  max_work_depth;
    int                  shutdown_requested;
    MPI_Status           status, msg_status;
    sf_work_request_t *  incoming_requests = NULL;
    useconds_t           delay = 20;
    subfiling_context_t *context = get__subfiling_object(context_id);
    double               queue_start_time;

    assert(context != NULL);
    /* We can't have opened any files at this point.. 
     * The file open approach has changed so that the normal
     * application rank (hosting this thread) does the file open.
     * We can simply utilize the file descriptor (which should now
     * represent an open file).
     */
    // context->sf_fid = -1;

    subfile_rank = context->sf_group_rank;
    sf_world_size = context->topology->app_layout->world_size;
    sf_world_rank = context->topology->app_layout->world_rank;

    if (request_count_per_rank == NULL) {
        request_count_per_rank =
            (int *) calloc((size_t) sf_world_size, sizeof(int));
        assert(request_count_per_rank != NULL);
    }

    max_work_depth = MAX(8, sf_world_size * MAX_WORK_PER_RANK);
    incoming_requests = (sf_work_request_t *) calloc(
        (size_t)(max_work_depth + 1), sizeof(sf_work_request_t));

    /* Validate that the allocation succeeded */
    assert(incoming_requests != NULL);

    /* Initialize atomic vars */
    atomic_init(&sf_workinprogress, 0);
    atomic_init(&sf_work_pending, 0);
    atomic_init(&sf_file_close_count, 0);
    atomic_init(&sf_file_refcount, 0);
    atomic_init(&sf_ioc_fini_refcount, 0);
    atomic_init(&sf_ioc_ready, 1);
    shutdown_requested = 0;

    while (!shutdown_requested || sf_work_pending) {
        flag = 0;
        ret = MPI_Iprobe(
            MPI_ANY_SOURCE, MPI_ANY_TAG, context->sf_msg_comm, &flag, &status);
        if ((ret == MPI_SUCCESS) && (flag != 0)) {
            sf_work_request_t *msg = NULL;
            int                count;
            int                index;
            int                request_size = (int) sizeof(sf_work_request_t);
            int                source = status.MPI_SOURCE;
            int                tag = status.MPI_TAG;

            MPI_Get_count(&status, MPI_BYTE, &count);
            if (count > request_size) {
                msg = (sf_work_request_t *) malloc((size_t) count);
                ret = MPI_Recv(msg, count, MPI_BYTE, source, tag,
                    context->sf_msg_comm, &msg_status);
            } else {
                index = atomic_load(&sf_workinprogress);
                ret = MPI_Recv(&incoming_requests[index], count,
                    MPI_BYTE, source, tag, context->sf_msg_comm, &msg_status);
                if (MPI_SUCCESS == ret) {
                    int howmany = 0;
                    MPI_Get_count(&msg_status, MPI_BYTE, &howmany);
                    if (howmany != count) {
                        printf("%s: MPI_Recv completed %d bytes of %d\n", __func__, howmany, count);
                        fflush(stdout);
                    }
                }
            }
            queue_start_time = MPI_Wtime();
            if (ret == MPI_SUCCESS) {
                if (msg) {
                    printf("%s: non-std msg=(%p) from %d\n", __func__, msg, source);
                    fflush(stdout);

                    msg->source = source;
                    msg->subfile_rank = subfile_rank;
                    msg->context_id = context->sf_context_id;
                    msg->start_time = queue_start_time;
                    tpool_add_work(msg);
                } else {
                    incoming_requests[index].tag = tag;
                    incoming_requests[index].source = source;
                    incoming_requests[index].subfile_rank = subfile_rank;
                    incoming_requests[index].start_time = queue_start_time;
					incoming_requests[index].buffer = NULL;
					incoming_requests[index].completed = 0;
                    tpool_add_work(&incoming_requests[index]);
                    if (index == max_work_depth - 1) {
                        atomic_init(&sf_workinprogress, 0);
                    } else {
                        atomic_fetch_add(&sf_workinprogress, 1); // atomic
                    }
                }
            }
        } else {
            usleep(delay);
        }
        shutdown_requested = sf_shutdown_flag;
    }

    if (incoming_requests) {
        free(incoming_requests);
    }

    /* Reset the shutdown flag */
    sf_shutdown_flag = 0;

    return 0;
}

/*
=========================================
Private helper functions
=========================================
*/

static int
send_ack__(int target, int subfile_rank, int tag, MPI_Comm comm)
{
    int ack = 1;
    int ret = MPI_Send(&ack, 1, MPI_INT, target, tag, comm);
#ifndef NDEBUG
    if (sf_verbose_flag) {
        if (sf_logfile) {
            fprintf(sf_logfile, "[ioc(%d): Sending ACK to MPI_rank(%d)\n",
                subfile_rank, target);
        }
    }
#endif
    return ret;
}

static int
send_nack__(int target, int subfile_rank, int tag, MPI_Comm comm)
{
    int nack = 0;
    int ret = MPI_Send(&nack, 1, MPI_INT, target, tag, comm);

#ifndef NDEBUG
    if (sf_verbose_flag) {
        if (sf_logfile) {
            fprintf(sf_logfile, "[ioc(%d): Sending NACK to MPI_rank(%d)\n",
                subfile_rank, target);
        }
    }
#endif
    return ret;
}

/*
=========================================
queue_xxx functions that should be run
from the thread pool threads...
=========================================
*/


/*-------------------------------------------------------------------------
 * Function:    Public/IOC queue_write_indep
 *
 * Purpose:     Implement the IOC independent write function.  The
 *              function is invoked as a result of the IOC receiving the
 *              "header"/RPC.  What remains is to allocate memory for the
 *              data sent by the client and then write the data to our
 *              subfile.  We utilize pwrite for the actual file writing.
 *              File flushing is done at file close.
 *
 * Return:      The integer status returned by the Internal read_independent
 *              function.  Successful operations will return 0.
 * Errors:      An MPI related error value.
 *
 * Programmer:  Richard Warren
 *              7/17/2020
 *
 * Changes:     Initial Version/None.
 *
 *-------------------------------------------------------------------------
 */
int
queue_write_indep(
    sf_work_request_t *msg, int subfile_rank, int source, MPI_Comm comm)
{
    int                  fd;
    char *               recv_buffer = NULL;
    int                  ret = MPI_SUCCESS;
    MPI_Status           msg_status;
    int64_t              data_size = msg->header[0];
    int64_t              file_offset = msg->header[1];
    int64_t              file_context_id = msg->header[2];
    double               t_start, t_end;
    double               t_write, t_wait, t_queue_delay;
    subfiling_context_t *sf_context = get__subfiling_object(file_context_id);
    assert(sf_context != NULL);

    /* flag that we've attempted to write data to the file */
    sf_context->sf_write_count++;
    /* For debugging performance */
    sf_write_ops++;

    t_start = MPI_Wtime();
    t_queue_delay = t_start - msg->start_time;

#ifndef NDEBUG
    if (sf_verbose_flag) {
        if (sf_logfile) {
            fprintf(sf_logfile,
            "[ioc(%d) %s]: msg from %d: datasize=%ld\toffset=%ld, queue_delay = %lf seconds\n",
            subfile_rank, __func__, source, data_size, file_offset, t_queue_delay);
        }
    }
#endif

    if (recv_buffer == NULL) {
        if ((recv_buffer = (char *) malloc((size_t) data_size)) == NULL) {
            perror("malloc");
            send_nack__(source, subfile_rank, WRITE_INDEP_ACK, comm);
            return -1;
        }
    }

    send_ack__(source, subfile_rank, WRITE_INDEP_ACK, comm);
    ret = MPI_Recv(recv_buffer, (int) data_size, MPI_BYTE, source,
        WRITE_INDEP_DATA, comm, &msg_status);

    t_end = MPI_Wtime();
    t_wait = t_end - t_start;
    sf_write_wait_time += t_wait;
    t_start = t_end;
#ifndef NDEBUG
    if (sf_verbose_flag) {
        if (sf_logfile) {
            fprintf(sf_logfile,
                "[ioc(%d) %s] MPI_Recv(%ld bytes, from = %d) status = %d\n",
                subfile_rank, __func__, data_size, source, ret);
        }
    }
#endif

    if (ret != MPI_SUCCESS) {
        int  len;
        char estring[MPI_MAX_ERROR_STRING];
        MPI_Error_string(ret, estring, &len);
        printf("[ioc(%d) %s] MPI_ERROR(%d)! MPI_Recv of %ld bytes from %d "
               "returned an error(%s)\n",
            subfile_rank, __func__, msg_status.MPI_ERROR, data_size, source,
            estring);
        fflush(stdout);
        return ret;
    }

    fd = sf_context->sf_fid;
    if (fd < 0) {
        printf("[ioc(%d)] WARNING: %s called while subfile_fid = %d (closed)\n",
            subfile_rank, __func__, fd);
        fflush(stdout);
    } else {
        if (sf_write_data(fd, file_offset, recv_buffer, data_size, subfile_rank) < 0) {
            free(recv_buffer);
            recv_buffer = NULL;
            printf("[ioc(%d) %s] sf_write_data returned an error!\n", subfile_rank,
                   __func__);
           fflush(stdout);
           return -1;
        }
        t_end = MPI_Wtime();
        t_write = t_end - t_start;
        sf_pwrite_time += t_write;
    }

    sf_queue_delay_time += t_queue_delay;

    /* Done... */
	msg->completed = 1;
    if (recv_buffer) {
        free(recv_buffer);
    }
    return 0;
}

/*-------------------------------------------------------------------------
 * Function:    Public/IOC queue_read_indep
 *
 * Purpose:     Implement the IOC independent read function.  The
 *              function is invoked as a result of the IOC receiving the
 *              "header"/RPC.  What remains is to allocate memory for
 *              reading the data and then to send this to the client.
 *              We utilize pread for the actual file reading.
 *
 * Return:      The integer status returned by the Internal read_independent
 *              function.  Successful operations will return 0.
 * Errors:      An MPI related error value.
 *
 * Programmer:  Richard Warren
 *              7/17/2020
 *
 * Changes:     Initial Version/None.
 *
 *-------------------------------------------------------------------------
 */
int
queue_read_indep(
    sf_work_request_t *msg, int subfile_rank, int source, MPI_Comm comm)
{
    int                  fd;
    char *               send_buffer = NULL;
    int                  ret = MPI_SUCCESS;
    int64_t              data_size = msg->header[0];
    int64_t              file_offset = msg->header[1];
    int64_t              file_context_id = msg->header[2];
    double               t_start, t_end;
    double               t_read, t_queue_delay;

    subfiling_context_t *sf_context = get__subfiling_object(file_context_id);
    assert(sf_context != NULL);

    sf_context->sf_read_count++;
    /* For debugging performance */
    sf_read_ops++;

    t_start = MPI_Wtime();
    t_queue_delay = t_start - msg->start_time;

    fd = sf_context->sf_fid;
    if (fd < 0) {
        printf("[ioc(%d) %s] subfile(%d) file descriptor not valid\n",
            subfile_rank, __func__, fd);
        return -1;
    }

#ifndef NDEBUG
    if (sf_verbose_flag && (sf_logfile != NULL)) {
        fprintf(sf_logfile,
                "[ioc(%d) %s] msg from %d: datasize=%ld\toffset=%ld queue_delay=%lf seconds\n",
                subfile_rank, __func__, source, data_size, file_offset, t_queue_delay);
    }
#endif
    if ((send_buffer = (char *) malloc((size_t) data_size)) == NULL) {
        perror("malloc");
        return -1;
    }

    if (sf_read_data(fd, file_offset, send_buffer, data_size, subfile_rank) < 0) {
        printf("[%d] %s - sf_read_data for source(%d) returned an error! "
               "read_count=%ld\n",
            subfile_rank, __func__, source, sf_context->sf_read_count);
        fflush(stdout);
        return -1;
    }

    ret = MPI_Send(send_buffer, (int) data_size, MPI_BYTE, source, READ_INDEP_DATA, comm);
    if (ret != MPI_SUCCESS) {
        int  len;
        char estring[MPI_MAX_ERROR_STRING];
        MPI_Error_string(ret, estring, &len);
        printf("[ioc(%d)] ERROR! MPI_Send of %ld bytes to %d returned an "
               "error(%s)\n",
            subfile_rank, data_size, source, estring);
        fflush(stdout);
        return ret;
    }
    t_end = MPI_Wtime();
    t_read = t_end - t_start;
    sf_pread_time += t_read;
    sf_queue_delay_time += t_queue_delay;
#if 0
    printf("[ioc(%d)] MPI_Send %ld bytes to source(%d) completed\n",
            subfile_rank, data_size, source);
    fflush(stdout);
#endif

#ifndef NDEBUG
    if (sf_verbose_flag && (sf_logfile != NULL)) {
        fprintf(sf_logfile, "[ioc(%d)] MPI_Send to source(%d) completed\n",
                subfile_rank, source);
    }
#endif

    if (send_buffer) {
        free(send_buffer);
        send_buffer = NULL;
    }

    return 0;
}


/*-------------------------------------------------------------------------
 * Function:    Public/IOC subfiling_open_file
 *
 * Purpose:     This function gets called when a client invokes a OPEN_OP.
 *              The HDF5 file opening protocol actually attempts to open
 *              a file; first without any truncate other flags which would
 *              modify the file state if it already exists.  A file close
 *              and then the second file open using the user supplied open
 *              flags is invoked.   The OPEN_OP provides the user flags as
 *              part of the RPC message.  The file prefix info doesn't
 *              transmited as part of the RPC since it is available as
 *              part of the client context which can be utilized by the
 *              IOC thread.  We access the sf_context by reading the
 *              cache of contexts at the index provided with the RPC msg.
 *
 * Return:      The integer status returned by the Internal read_independent
 *              function.  Successful operations will return 0.
 * Errors:      An MPI related error value.
 *
 * Programmer:  Richard Warren
 *              7/17/2020
 *
 * Changes:     Initial Version/None.
 *
 *-------------------------------------------------------------------------
 */
int
subfiling_open_file( sf_work_request_t *msg, int subfile_rank, int flags)
{
    int errors = 0;
    int subfile_fid = -1;
    char filepath[PATH_MAX];

    char *prefix = NULL;
    char *subfile_dir = NULL;

    double t_start=0.0, t_end=0.0;
    /* Only the real IOCs open the subfiles
     * Once a file is opened, all additional file open requests
     * can return immediately.
     */

    t_start = MPI_Wtime();
    /* Only allow the actual IO concentrator ranks to create sub-files */
    if (subfile_rank >= 0) {
        char                 config[PATH_MAX];
        int64_t              h5_file_id = msg->header[1];
        int64_t              file_context_id = msg->header[2];
        subfiling_context_t *sf_context = get__subfiling_object(file_context_id);
        assert(sf_context != NULL);

        begin_thread_exclusive();
        // printf("%s (before) sf_fid = %d\n", __func__, sf_context->sf_fid );
        /* Check to see whether we need to create the subfile(s) */
        if (sf_context->sf_fid < 0) {
            int  n_io_concentrators = sf_context->topology->n_io_concentrators;
            int *io_concentrator = sf_context->topology->io_concentrator;
            const char *dotconfig = ".subfile_config";
            mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;

            if ((prefix = sf_context->subfile_prefix)) {
                mkdir(prefix, S_IRWXU);
                sprintf(filepath, "%s/" SF_FILENAME_TEMPLATE, prefix,
                    h5_file_id, subfile_rank, n_io_concentrators);
            } else {
              strcpy(filepath,sf_context->filename);
              subfile_dir = strrchr(filepath, '/');
              assert(subfile_dir);
              sprintf(subfile_dir+1, SF_FILENAME_TEMPLATE, h5_file_id,
                      subfile_rank, n_io_concentrators);
            }
            if ((sf_context->sf_fid = open(filepath, flags, mode)) < 0) {
                end_thread_exclusive();
                // printf("[%d %s] file open(%s) failed!\n", subfile_rank, __func__, filepath);
                // fflush(stdout);

#ifndef NDEBUG
                if (sf_verbose_flag) {
                    printf("[%d %s] file open(%s) failed!\n", subfile_rank, __func__, filepath);
                    fflush(stdout);
                }
#endif
                errors++;
                goto done;
            }

            // printf("%s: (after) sf_fid = %d\n", __func__, sf_context->sf_fid); 
            // fflush(stdout);

            strcpy(filepath,sf_context->filename);
            subfile_dir = strrchr(filepath, '/');
            if (subfile_dir) {
                sprintf(subfile_dir,"/%ld%s", h5_file_id, dotconfig);
                        strcpy(config, filepath);
            }

            if ((subfile_rank == 0) && (flags & O_CREAT)) {
                FILE *f = NULL;
                                
                /* If a config file already exists, AND
                 * the user wants to truncate subfiles (if they exist),
                 * then we should also truncate an existing config file.
                 */
                if (access(config, flags) == 0) {
                    truncate(config, 0);
                }
                f = fopen(config, "w+");
                if (f != NULL) {
                    int k;
                    char linebuf[PATH_MAX];
                    sprintf( linebuf, "stripe_size=%ld\n", sf_context->sf_stripe_size);
                    fwrite(linebuf, strlen(linebuf), 1, f);
                    sprintf(linebuf, "aggregator_count=%d\n", n_io_concentrators);
                    fwrite(linebuf, strlen(linebuf), 1, f);
                    sprintf(linebuf,"hdf5_file=%s\n", sf_context->filename);
                    fwrite(linebuf, strlen(linebuf), 1, f);

                    for (k = 0; k < n_io_concentrators; k++) {
                        if (prefix)
                            sprintf(linebuf, "%s/%ld_node_local_temp_%d_of_%d:%d\n", prefix,
                            h5_file_id, subfile_rank, n_io_concentrators, io_concentrator[k]);
                        else
                            sprintf(linebuf, "%ld_node_local_temp_%d_of_%d:%d\n", h5_file_id,
                            subfile_rank, n_io_concentrators, io_concentrator[k]);

                        fwrite(linebuf, strlen(linebuf), 1, f);
                    }

                    fclose(f);
                } else {
                    perror("fopen(config)");
                    errors++;
                    goto done;
                }
            }
#ifndef NDEBUG
            if (sf_verbose_flag) {
                if (sf_logfile) {
                    fprintf(sf_logfile, "[ioc:%d] Opened subfile %s\n",
                        subfile_rank, filepath);
                }
            }
#endif
        }
        end_thread_exclusive();
    }
done:
    t_end = MPI_Wtime();

#ifndef NDEBUG
    if (sf_verbose_flag) {
        printf("[%s %d] open completed in %lf seconds with %d errors\n",
               __func__, subfile_rank, (t_end - t_start), errors);
        fflush(stdout);
    }
#endif
    return errors;
}

/*-------------------------------------------------------------------------
 * Function:    UTILITY FUNCTIONS:
 *
 *              sf_get_mpi_rank  - (not used) retrieves the MPI rank of the
 *                                 calling process.  Was used when pairing
 *                                 the subfiling VFD with the SUBFILING VFD.
 *
 *              sf_get_mpi_size  - (not used) retrieves the MPI size of the
 *                                 communicator associated with the open
 *                                 file.
 *
 *              sf_get_group_com - (not used) retrieves the MPI Comm object
 *                                 associated with the open file/sf_context.
 *
 *              sf_subfile_set_logging - (not used) informs one or all IOC
 *                                 instances to set the verbose/logging flag
 *                                 to the value provided by the user.
 *
 * Return:      none
 * Errors:      none
 *
 * Programmer:  Richard Warren
 *              7/17/2020
 *
 * Changes:     Initial Version/None.
 *
 *-------------------------------------------------------------------------
 */

int
sf_get_mpi_rank(hid_t fid, int *rank)
{
    hid_t                context_id = fid_map_to_context(fid);
    subfiling_context_t *sf_context = get__subfiling_object(context_id);
    assert(sf_context != NULL);
    assert(rank != NULL);
    *rank = sf_context->sf_group_rank;
    return 0;
}

int
sf_get_mpi_size(hid_t fid, int *size)
{
    hid_t                context_id = fid_map_to_context(fid);
    subfiling_context_t *sf_context = get__subfiling_object(context_id);
    assert(sf_context != NULL);
    assert(size != NULL);
    *size = sf_context->sf_group_size;
    return 0;
}

int
sf_get_group_comm(hid_t fid, MPI_Comm *comm)
{
    hid_t                context_id = fid_map_to_context(fid);
    subfiling_context_t *sf_context = get__subfiling_object(context_id);
    assert(sf_context != NULL);
    assert(comm != NULL);
    *comm = sf_context->sf_group_comm;
    return 0;
}

int
sf_subfile_set_logging(hid_t sf_fid, int ioc_rank, int flag)
{
    int                  ioc;
    int                  status = 0;
    hid_t                context_id = fid_map_to_context(sf_fid);
    subfiling_context_t *sf_context = get__subfiling_object(context_id);
    int                  n_io_concentrators;
    int *                io_concentrator = NULL;
    int64_t              lflag = (int64_t)(flag & 0xFF);
    int64_t              msg[3];

    assert(sf_context != NULL);

    msg[0] = lflag;
    msg[1] = 0;
    msg[2] = sf_context->sf_context_id;

    n_io_concentrators = sf_context->topology->n_io_concentrators;
    io_concentrator = sf_context->topology->io_concentrator;

    for (ioc = 0; ioc < n_io_concentrators; ioc++) {
        if ((flag < 0) || (flag == ioc_rank)) {
            status = MPI_Ssend(msg, 3, MPI_INT64_T, io_concentrator[ioc],
                LOGGING_OP, sf_context->sf_msg_comm);
        }
    }
    return status;
}
