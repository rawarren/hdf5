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

#include "H5private.h"
#include "h5diff.h"
#include "ph5diff.h"
#include "h5diff_common.h"
#include "h5tools.h"
#include "h5tools_utils.h"

/* Name of tool */
#define PROGRAMNAME "ph5diff"

/* List of table structures.  There is one table structure for each file */

#ifndef PATH_MAX
#define PATH_MAX 512
#endif
#ifndef MAX_NAME
#define MAX_NAME 256
#endif

/* File_Access_type bits */
#define FACC_DEFAULT 0x0 /* default */
#define FACC_MPIO    0x1 /* MPIO */
#define FACC_SPLIT   0x2 /* Split File */

#define DXFER_COLLECTIVE_IO  0x1 /* Collective IO*/
#define DXFER_INDEPENDENT_IO 0x2 /* Independent IO collectively */

/* Name of tool */
#define PROGRAMNAME "ph5diff"

static herr_t     gather_obj_info(hid_t file, const char *oname, iter_t *iter);
static void       ph5diff_worker(int);
static size_t     dset_compare(hid_t dset_obj1, hid_t dset_obj2, int ndims, hsize_t *dims, size_t n_elements,
                               int mpi_size, int mpi_rank);
static void       diff_single(int mpi_rank, int mpi_size);
static diff_err_t diff_collection(void);

/* What if the ndims counts aren't equal?
 * That doesn't mean the data isn't the same...
 * e.g. [1000,1,1,1] is the same as [1000]!
 * Should we just compare elements?
 */
static int
compare_dims(dset_t *ds1, dset_t *ds2)
{
    if (ds1->ndims == ds2->ndims) {
        int count;
        int k;
        for (k = 0, count = 0; k < ds1->ndims; k++) {
            if (ds1->dims[k] != ds2->dims[k])
                count++;
        }
        return count;
    }
    return 1;
}

static int
get_local_dims(size_t total_rows, int mpi_size, int mpi_rank, int row_index, int ndims, size_t *start_dims,
               size_t *block_dims, size_t *dims)
{
    size_t local_rows = (size_t)(total_rows / (size_t)mpi_size);
    size_t rank_start = (size_t)(local_rows * (size_t)mpi_rank);

    if (ndims > 3) {
        printf("ndims > 3:  what now?");
        return -1;
    }
    else if (ndims > 2) {
        start_dims[0]         = rank_start / dims[row_index];
        block_dims[0]         = local_rows / dims[row_index];
        start_dims[row_index] = rank_start % dims[row_index];
        block_dims[row_index] =
            (local_rows % dims[row_index] != 0 ? local_rows % dims[row_index] : dims[row_index]);
        start_dims[ndims - 1] = 0;
        block_dims[ndims - 1] = dims[ndims - 1];
    }
    else {
        /* 2D */
        start_dims[0]         = rank_start / dims[row_index];
        block_dims[0]         = local_rows;
        block_dims[ndims - 1] = dims[ndims - 1];
    }
    block_dims[ndims - 1] = dims[ndims - 1];
    return 0;
}

static size_t
get_rank_dims(int mpi_size, int mpi_rank /*, int row_index */, int ndims, hsize_t *start_dims,
              hsize_t *block_dims, hsize_t *dims)
{
    int     k;
    hsize_t blocksize = (hsize_t)(dims[0] / (hsize_t)mpi_size);
    hsize_t remaining = (hsize_t)(dims[0] % (hsize_t)mpi_size);

    block_dims[0] = blocksize;
    start_dims[0] = (blocksize * (hsize_t)mpi_rank);
    if (remaining && (mpi_rank == (mpi_size - 1))) {
        block_dims[0] += remaining;
        blocksize += remaining;
    }
    for (k = 1; k < ndims; k++) {
        block_dims[k] = dims[k];
        start_dims[k] = 0;
        blocksize *= dims[k];
    }
    return blocksize;
}

static size_t
dset_compare(hid_t dset_obj1, hid_t dset_obj2, int ndims, hsize_t *dims, size_t n_elements, int mpi_size,
             int mpi_rank)
{
    herr_t ret;
    hid_t  dtype           = -1;
    hid_t  file_dataspace1 = -1, file_dataspace2 = -1; /* File dataspace IDs */
    hid_t  mem_dataspace1 = -1, mem_dataspace2 = -1;
    hid_t  xfer_pl     = -1;
    int    total_diffs = 0;
    int    diffs       = 0;
    double t0, t1, io_time;

    if (mpi_size > 1) {
        int     k, row_index = ndims - 2;
        size_t  dtype_size    = 0;
        size_t  rank_elements = 0;
        size_t  total_rows;
        size_t  slice_size     = 0;
        hsize_t start_dims[32] = {
            0,
        };
        hsize_t block_dims[32] = {
            0,
        };
        hsize_t count_dims[32] = {
            1,
        };
        void *rdata1 = NULL, *rdata2 = NULL;

        t0         = MPI_Wtime();
        total_rows = dims[row_index];
        for (k = 0; k < row_index; k++)
            total_rows *= dims[k];

        for (k = 0; k < ndims; k++)
            count_dims[k] = 1;

        slice_size = get_rank_dims(mpi_size, mpi_rank /*, row_index */, ndims, start_dims, block_dims, dims);
        if ((xfer_pl = H5Pcreate(H5P_DATASET_XFER)) < 0) {
            printf("Unable to create a dset xfer property!\n");
        }
        if (H5Pset_dxpl_mpio(xfer_pl, H5FD_MPIO_INDEPENDENT) < 0) {
            printf("Unable to set the file xfer property!\n");
        }

        file_dataspace1 = H5Dget_space(dset_obj1);
        H5Sselect_hyperslab(file_dataspace1, H5S_SELECT_SET, start_dims, NULL, count_dims, block_dims);
        mem_dataspace1 = H5Screate_simple(3, block_dims, NULL);
        dtype          = H5Dget_type(dset_obj1);
        if (dtype > 0)
            dtype_size = H5Tget_size(dtype);
        if (dtype_size > 0)
            rdata1 = malloc(slice_size * dtype_size);
        if (rdata1 == NULL) {
            printf("Unable to allocate buffer for reading\n");
            goto done;
        }

        ret = H5Dread(dset_obj1, dtype, mem_dataspace1, file_dataspace1, xfer_pl, rdata1);
        if (ret < 0) {
            printf("H5Dread dset_obj1: returned an error!\n");
            goto done;
        }

        file_dataspace2 = H5Dget_space(dset_obj2);
        H5Sselect_hyperslab(file_dataspace2, H5S_SELECT_SET, start_dims, NULL, count_dims, block_dims);
        mem_dataspace2 = H5Screate_simple(3, block_dims, NULL);

        rdata2 = malloc(slice_size * dtype_size);
        if (rdata2 == NULL) {
            printf("Unable to allocate buffer for reading\n");
            goto done;
        }

        ret = H5Dread(dset_obj2, dtype, mem_dataspace2, file_dataspace2, xfer_pl, rdata2);
        if (ret < 0) {
            printf("H5Dread dset_obj2: returned an error!\n");
            goto done;
        }
        t1      = MPI_Wtime();
        io_time = t1 - t0;
        t0      = t1;
        if (dtype_size == 2) {
            unsigned short *i2_dset1 = (unsigned short *)rdata1;
            unsigned short *i2_dset2 = (unsigned short *)rdata2;
            size_t          index;
            for (index = 0; index < rank_elements; index++) {
                if (i2_dset1[index] != i2_dset2[index])
                    diffs++;
            }
        }
        else if (dtype_size == 4) {
            unsigned int *i4_dset1 = (unsigned int *)rdata1;
            unsigned int *i4_dset2 = (unsigned int *)rdata2;
            size_t        index;
            for (index = 0; index < rank_elements; index++) {
                if (i4_dset1[index] != i4_dset2[index])
                    diffs++;
            }
        }
        else if (dtype_size == 8) {
            uint64_t *i8_dset1 = (uint64_t *)rdata1;
            uint64_t *i8_dset2 = (uint64_t *)rdata2;
            size_t    index;
            for (index = 0; index < rank_elements; index++) {
                if (i8_dset1[index] != i8_dset2[index])
                    diffs++;
            }
        }

        if (MPI_Allreduce(&diffs, &total_diffs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
            printf("MPI_Allreduce returned an error\n");
        }
        t1 = MPI_Wtime();
        if (mpi_rank == 0) {
            double comparison_time = t1 - t0;
            double total_time      = io_time + comparison_time;
            printf("total_diffs = %d\n", total_diffs);
            printf("total time = %lf == (io+cmp) %lf + %lf\n", total_time, io_time, comparison_time);
        }

done:
        if (rdata1 != NULL)
            free(rdata1);
        if (rdata2 != NULL)
            free(rdata2);
        if (file_dataspace1 >= 0)
            H5Sclose(file_dataspace1);
        if (file_dataspace2 >= 0)
            H5Sclose(file_dataspace2);
        if (mem_dataspace1 >= 0)
            H5Sclose(mem_dataspace1);
        if (mem_dataspace2 >= 0)
            H5Sclose(mem_dataspace2);
        if (xfer_pl >= 0)
            H5Pclose(xfer_pl);
    }
    return (size_t)diffs;
}

static size_t
dset_compare_seq(hid_t dset_obj1, hid_t dset_obj2, int ndims, hsize_t *dims, size_t n_elements)
{
    printf("Entered dset_compare_seq\n");
    return 0;
}

static void
print_class(H5T_class_t obj_class, hid_t type1, int indent)
{
    if (indent > 0)
        printf("\t");
    switch (obj_class) {
        case H5T_INTEGER:
            puts("H5T_INTEGER");
            break;
        case H5T_FLOAT:
            puts("H5T_FLOAT");
            break;
        case H5T_TIME:
            puts("H5T_TIME");
            break;
        case H5T_STRING:
            puts("H5T_STRING");
            break;
        case H5T_BITFIELD:
            puts("H5T_BITFIELD");
            break;
        case H5T_OPAQUE:
            puts("H5T_OPAQUE");
            break;
        case H5T_COMPOUND:
            puts("H5T_COMPOUND");
            {
                int i, nmembers = H5Tget_nmembers(type1);
                for (i = 0; i < nmembers; i++) {
                    char *mname = H5Tget_member_name(type1, (unsigned)i);
                    hid_t mtype = H5Tget_member_type(type1, (unsigned)i);
                    print_class(H5Tget_class(mtype), mtype, (unsigned)1);
                    printf("member_name = %s\n", mname);
                }
            }
            break;
        case H5T_REFERENCE:
            puts("H5T_REFERENCE");
            break;
        case H5T_ENUM:
            puts("H5T_ENUM");
            break;
        case H5T_VLEN:
            puts("H5T_VLEN");
            break;
        case H5T_ARRAY:
            puts("H5T_ARRAY");
            break;
        default:
            puts("Unknown class\n");
    }
}

static size_t
dset_attr_compare(hid_t obj1, hid_t obj2, H5O_info_t *objInfo1, H5O_info_t *objInfo2)
{
    if (objInfo1->num_attrs == objInfo2->num_attrs) {
        int         i, rank1, rank2;
        ssize_t     len1, len2;
        size_t      type1_extent, type2_extent;
        hid_t       attr1, attr2;
        hid_t       atype1, atype2;
        hid_t       aspace1, aspace2;
        char        attr1_name[64], attr2_name[64];
        H5T_class_t type_class1, type_class2;

        /* The following assumes that the iteration order results
         * in the same attribute: 'name','type', etc. from each dataset...
         */
        for (i = 0; i < (int)objInfo1->num_attrs; i++) {
            attr1 = H5Aopen_by_idx(obj1, ".", H5_INDEX_CRT_ORDER, H5_ITER_INC, (hsize_t)i, H5P_DEFAULT,
                                   H5P_DEFAULT);
            attr2 = H5Aopen_by_idx(obj2, ".", H5_INDEX_CRT_ORDER, H5_ITER_INC, (hsize_t)i, H5P_DEFAULT,
                                   H5P_DEFAULT);
            len1  = H5Aget_name(attr1, sizeof(attr1_name), attr1_name);
            len2  = H5Aget_name(attr2, sizeof(attr2_name), attr2_name);
            if (len1 == len2) {
                printf("comparing: %s , %s\n", attr1_name, attr2_name);
                atype1 = H5Aget_type(attr1);
                atype2 = H5Aget_type(attr2);

                type_class1 = H5Tget_class(atype1);
                type_class2 = H5Tget_class(atype2);
                if (type_class1 == type_class2) {
                    print_class(type_class1, atype1, 0);
                }

                aspace1 = H5Aget_space(attr1);
                aspace2 = H5Aget_space(attr2);

                rank1 = H5Sget_simple_extent_ndims(aspace1);
                rank2 = H5Sget_simple_extent_ndims(aspace2);

                type1_extent = H5Tget_size(atype1);
                type2_extent = H5Tget_size(atype2);
                printf("type_extents = %ld, %ld\n", type1_extent, type2_extent);
            }
            else {
                printf("Attribute ordering may different!  Attribute names are different!\n");
            }
            H5Aclose(attr1);
            H5Aclose(attr2);
        }
    }
    return 0;
}

static size_t
dset_compare_contig(hid_t obj1, hid_t obj2, size_t n_elements, size_t dtype_size, int mpi_size, int mpi_rank)
{
    hid_t xfer_plist;     /* Dataset transfer properties list */
    hid_t file_dataspace; /* File dataspace ID */
    hid_t mem_dataspace;  /* memory dataspace ID */
    if (mpi_size > 1) {
        size_t ranks          = (size_t)mpi_size;
        size_t slice_elements = n_elements / ranks;
        if (n_elements % ranks) {
            slice_elements++;
        }
    }
    return 0;
}

static void
dtype_list1(hid_t tid)
{
    H5T_class_t t_class;
    t_class = H5Tget_class(tid);
    if (t_class < 0) {
        puts(" Invalid datatype.\n");
    }
    else {
        /*
         * Each class has specific properties that can be
         * retrieved, e.g., size, byte order, exponent, etc.
         */
        if (t_class == H5T_INTEGER) {
            puts(" Datatype is 'H5T_INTEGER'.\n");
            /* display size, signed, endianess, etc. */
        }
        else if (t_class == H5T_FLOAT) {
            puts(" Datatype is 'H5T_FLOAT'.\n");
            /* display size, endianess, exponennt, etc. */
        }
        else if (t_class == H5T_STRING) {
            puts(" Datatype is 'H5T_STRING'.\n");
            /* display size, padding, termination, etc. */
        }
        else if (t_class == H5T_BITFIELD) {
            puts(" Datatype is 'H5T_BITFIELD'.\n");
            /* display size, label, etc. */
        }
        else if (t_class == H5T_OPAQUE) {
            puts(" Datatype is 'H5T_OPAQUE'.\n");
            /* display size, etc. */
        }
        else if (t_class == H5T_COMPOUND) {
            puts(" Datatype is 'H5T_COMPOUND'.\n");
            /* recursively display each member: field name, type  */
        }
        else if (t_class == H5T_ARRAY) {
            puts(" Datatype is 'H5T_COMPOUND'.\n");
            /* display  dimensions, base type  */
        }
        else if (t_class == H5T_ENUM) {
            puts(" Datatype is 'H5T_ENUM'.\n");
            /* display elements: name, value   */
        }
        else {
            puts(" Datatype is 'Other'.\n");
            /* eg. Object Reference, ...and so on ... */
        }
    }
}

static herr_t
group_list1(hid_t gid)
{
    int     i;
    ssize_t len;
    hsize_t nobj;
    herr_t  err;
    int     otype;
    hid_t   grpid, typeid, dsid;
    char    group_name[MAX_NAME];
    char    memb_name[MAX_NAME];

    /*
     * Information about the group:
     *  Name and attributes
     *
     *  Other info., not shown here: number of links, object id
     */
    len = H5Iget_name(gid, group_name, MAX_NAME);
    printf("Group Name: %s\n", group_name);

    /*
     *  Get all the members of the groups, one at a time.
     */
    err = H5Gget_num_objs(gid, &nobj);
    for (i = 0; i < nobj; i++) {
        /*
         *  For each object in the group, get the name and
         *   what type of object it is.
         */
        printf("  Member: %d ", i);
        fflush(stdout);
        len = H5Gget_objname_by_idx(gid, (hsize_t)i, memb_name, (size_t)MAX_NAME);
        printf("   %d ", len);
        fflush(stdout);
        printf("  Member: %s ", memb_name);
        fflush(stdout);
        otype = H5Gget_objtype_by_idx(gid, (size_t)i);

        /*
         * process each object according to its type
         */
        switch (otype) {
            case H5G_GROUP:
                printf(" GROUP:\n");
                grpid = H5Gopen(gid, memb_name, H5P_DEFAULT);
                group_list1(grpid);
                H5Gclose(grpid);
                break;
            case H5G_DATASET:
                printf(" DATASET:\n");
                dsid = H5Dopen(gid, memb_name, H5P_DEFAULT);
                dataset_list1(dsid);
                H5Dclose(dsid);
                break;
            case H5G_TYPE:
                printf(" DATA TYPE:\n");
                typeid = H5Topen(gid, memb_name, H5P_DEFAULT);
                dtype_list1(typeid);
                H5Tclose(typeid);
                break;
            default:
                printf(" unknown?\n");
                break;
        }
    }
}

/*-------------------------------------------------------------------------
 * Function:    datatype_list2
 *
 * Purpose:     List information about a datatype which should appear after
 *              information which is general to all objects.
 *
 * Return:      Success: 0
 *              Failure: -1
 *-------------------------------------------------------------------------
 */
static herr_t
datatype_list2(hid_t type, const char H5_ATTR_UNUSED *name)
{
    printf("datatype_list2 entered.  What now?\n");
    return 0;
}

/*-------------------------------------------------------------------------
 * Function:    list_obj
 *
 * Purpose:     Prints information about an object
 *
 * Return:      Success: 0
 *              Failure: -1
 *-------------------------------------------------------------------------
 */
static herr_t
list_obj(const char *name, const H5O_info_t *oinfo, const char *first_seen, void *_iter)
{
    H5O_type_t        obj_type = oinfo->type; /* Type of the object */
    iter_t *          iter     = (iter_t *)_iter;
    h5tools_str_t     buffer; /* string into which to render   */
    h5tools_context_t ctx;    /* print context  */
    dset_info_t *     dset = &dsets1;

    H5TOOLS_DEBUG("enter");

    HDmemset(&ctx, 0, sizeof(ctx));
    HDmemset(&buffer, 0, sizeof(h5tools_str_t));

    h5tools_str_reset(&buffer);
    if (dsets2.iter == iter)
        dset = &dsets2;

    H5TOOLS_DEBUG("Object name:%s", name);
    /* Print the link's name, either full name or base name */
    if (!iter->symlink_target) {
        if (dset->dset_count == dset->max_count) {
            alloc_dset_table(dset);
        }
    }

    /* Check object information */
    if (oinfo->type < 0 || oinfo->type >= H5O_TYPE_NTYPES) {
        h5tools_str_append(&buffer, "Unknown type(%d)", (int)oinfo->type);
        obj_type = H5O_TYPE_UNKNOWN;
    }
    if (iter->symlink_target)
        h5tools_str_append(&buffer, "{");
    if (obj_type >= 0 && dispatch_g[obj_type].name)
        h5tools_str_append(&buffer, "%s", dispatch_g[obj_type].name);

    /* Check if we've seen this object before */
    if (first_seen) {
        h5tools_str_reset(&buffer);
        h5tools_str_append(&buffer, ", same as ");
        printf("%s %s\n", buffer.s, first_seen);
    } /* end if */
    else {
        hid_t obj = H5I_INVALID_HID; /* ID of object opened */

        /* Open the object.  Not all objects can be opened.  If this is the case
         * then return right away.
         */
        H5TOOLS_DEBUG("Open object name=%s", name);
        if (obj_type >= 0 && (obj = H5Oopen(iter->fid, name, H5P_DEFAULT)) < 0) {
            h5tools_str_reset(&buffer);
            h5tools_str_append(&buffer, " *ERROR*\n");
            goto done;
        } /* end if */

        /* List the first line of information for the object. */
        H5TOOLS_DEBUG("Object type:%d", obj_type);

        if (obj_type == H5O_TYPE_GROUP) {
            iter_t *icopy = (iter_t *)malloc(sizeof(iter_t));
            memcpy(icopy, iter, sizeof(iter_t));
            H5Oclose(obj);
            gather_obj_info(iter->fid, name, iter);
        }
        else if (obj_type == H5O_TYPE_DATASET) {
            dset->known_dsets[dset->dset_count].name  = strdup(name);
            dset->known_dsets[dset->dset_count].token = oinfo->token;

            if (obj_type >= 0 && dispatch_g[obj_type].list1) {
                dset->new_dset = obj;
                (dispatch_g[obj_type].list1)(obj);
            }
            if (!iter->symlink_target) {
                h5tools_str_reset(&buffer);
                h5tools_str_append(&buffer, "\n");
            }

            /* Detailed list for object */
            if (obj_type >= 0 && dispatch_g[obj_type].list2)
                (dispatch_g[obj_type].list2)(obj, name);

            /* Close the object. */
            if (obj_type >= 0)
                H5Oclose(obj);
        } /* end else */
    }

done:
    if (iter->symlink_target) {
        h5tools_str_reset(&buffer);
        h5tools_str_append(&buffer, "}\n");
        iter->symlink_target = FALSE;
    }
    h5tools_str_close(&buffer);

    H5TOOLS_ENDDEBUG("exit");

    return 0;
} /* end list_obj() */

/*-------------------------------------------------------------------------
 * Function:    gather_obj_info
 *
 * Purpose:     Begins iteration on an object
 *
 * Return:      Success: 0
 *              Failure: -1
 *-------------------------------------------------------------------------
 */
static herr_t
gather_obj_info(hid_t file, const char *oname, iter_t *iter)
{
    int               retval = 0;
    H5O_info_t        oi;           /* Information for object */
    hsize_t           curr_pos = 0; /* total data element position   */
    h5tools_str_t     buffer;       /* string into which to render   */
    h5tools_context_t ctx;          /* print context  */

    HDmemset(&ctx, 0, sizeof(ctx));
    HDmemset(&buffer, 0, sizeof(h5tools_str_t));

    h5tools_str_reset(&buffer);

    /* Retrieve info for object to list */
    if (H5Oget_info_by_name(file, oname, &oi, H5O_INFO_BASIC | H5O_INFO_TIME, H5P_DEFAULT) < 0) {
        if (iter->symlink_target) {
            h5tools_str_append(&buffer, "{**NOT FOUND**}\n");
            iter->symlink_target = FALSE;
        }
        else
            printf("%s %s **NOT FOUND**\n", buffer.s, oname);
        retval = -1;
        goto done;
    } /* end if */

    /* Check for group iteration */
    if (H5O_TYPE_GROUP == oi.type) {
        /* Get ID for group */
        if (!iter->symlink_target && (iter->gid = H5Gopen2(file, oname, H5P_DEFAULT)) < 0) {
            h5tools_str_append(&buffer, "%s: unable to open '%s' as group\n", iter->fname, oname);
            goto done; /* Previously "continue", when this code was in main().
                        * We don't "continue" here in order to close the file
                        * and free the file name properly. */
        }              /* end if */

        /* Delay specifying the name start point so the original object name is
         * displayed if it is a link or non-group object */
        iter->name_start = iter->base_len;

        /* Specified name is a group. List the complete contents of the group. */
        h5trav_visit(file, oname, iter->symlink_target, FALSE, list_obj, NULL, iter, H5O_INFO_BASIC);

#if 0
        /* Close group */
        if (!iter->symlink_target)
            H5Gclose(iter->gid);
#endif
    } /* end if */
    else {
        /* Use file ID for root group ID */
        iter->gid = file;
        /* Specified name is a non-group object -- list that object */
        list_obj(oname, &oi, NULL, iter);
    } /* end else */

done:
    h5tools_str_close(&buffer);

    return retval;
}

#ifdef __NOT_USED__
#if 1
static hsize_t
create_object_tables(const char *fname1, const char *fname2, const char *objname1, const char *objname2,
                     diff_opt_t *opts)
{
    hid_t          file1_id = -1;
    hid_t          file2_id = -1;
    iter_t *       iter1, *iter2;
    static char    root_name[] = "/";
    symlink_trav_t symlink_list1, symlink_list2;
    diff_err_t     ret_value = H5DIFF_NO_ERR;
    H5O_info2_t    oi;

    char *oname = NULL;

    /*-------------------------------------------------------------------------
     * open the files first; if they are not valid, no point in continuing
     *-------------------------------------------------------------------------*/

    if ((file1_id = h5tools_fopen(fname1, H5F_ACC_RDONLY, H5P_DEFAULT, NULL, NULL, (size_t)0)) < 0) {
        parallel_print("h5diff: <%s>: unable to open file\n", fname1);
        H5TOOLS_GOTO_ERROR(H5DIFF_ERR, "<%s>: unable to open file\n", fname1);
    } /* end if */
    H5TOOLS_DEBUG("file1_id = %s", fname1);

    if ((file2_id = h5tools_fopen(fname2, H5F_ACC_RDONLY, H5P_DEFAULT, NULL, NULL, (size_t)0)) < 0) {
        parallel_print("h5diff: <%s>: unable to open file\n", fname2);
        H5TOOLS_GOTO_ERROR(H5DIFF_ERR, "<%s>: unable to open file\n", fname2);
    } /* end if */
    H5TOOLS_DEBUG("file2_id = %s", fname2);

    /* Build object display table */
    DISPATCH(H5O_TYPE_GROUP, "Group", NULL, NULL);
    DISPATCH(H5O_TYPE_DATASET, "Dataset", dataset_list1, dataset_list2);
    DISPATCH(H5O_TYPE_NAMED_DATATYPE, "Type", NULL, datatype_list2);

    iter1 = (iter_t *)malloc(sizeof(iter_t));
    iter2 = (iter_t *)malloc(sizeof(iter_t));
    initialize_dset_info(&dsets1, iter1);
    initialize_dset_info(&dsets2, iter2);

    iter1->fname = strdup(fname1);
    iter2->fname = strdup(fname2);
    iter1->fid   = file1_id;
    iter2->fid   = file2_id;
    iter1->gid = iter2->gid = -1;
    iter1->symlink_target = iter2->symlink_target = FALSE;
    iter1->symlink_list                           = &symlink_list1;
    iter2->symlink_list                           = &symlink_list2;
    iter1->symlink_list->dangle_link = iter2->symlink_list->dangle_link = FALSE;

    iter1->base_len = iter2->base_len = 0;
    iter1->name_start = iter2->name_start = 0;

    /* Initialize list of visited symbolic links */
    symlink_list1.nused = symlink_list1.nalloc = 0;
    symlink_list2.nused = symlink_list2.nalloc = 0;
    symlink_list1.objs = symlink_list2.objs = NULL;
    oname                                   = root_name;

    if (gather_obj_info(file1_id, oname, iter1) < 0) {
        printf("gather_obj_info returned an error\n");
    }
    if (gather_obj_info(file2_id, oname, iter2) < 0) {
        printf("gather_obj_info returned an error\n");
    }

done:
    if (file1_id > 0) {
        H5Fclose(file1_id);
        iter1->fid = -1;
    }
    if (file2_id > 0) {
        H5Fclose(file2_id);
        iter2->fid = -1;
    }
    return ret_value;
}

#else
static hsize_t
create_object_tables(const char *fname1, const char *fname2, const char *objname1, const char *objname2,
                     diff_opt_t *opts)
{
    hid_t          file1_id = -1;
    hid_t          file2_id = -1;
    iter_t *       iter1, *iter2;
    static char    root_name[] = "/";
    symlink_trav_t symlink_list1, symlink_list2;
    diff_err_t     ret_value = H5DIFF_NO_ERR;
    char *         oname     = NULL;

    /*-------------------------------------------------------------------------
     * open the files first; if they are not valid, no point in continuing
     *-------------------------------------------------------------------------*/

    if ((file1_id = h5tools_fopen(fname1, H5F_ACC_RDONLY, H5P_DEFAULT, NULL, NULL, (size_t)0)) < 0) {
        parallel_print("h5diff: <%s>: unable to open file\n", fname1);
        H5TOOLS_GOTO_ERROR(H5DIFF_ERR, "<%s>: unable to open file\n", fname1);
    } /* end if */
    H5TOOLS_DEBUG("file1_id = %s", fname1);

    if ((file2_id = h5tools_fopen(fname2, H5F_ACC_RDONLY, H5P_DEFAULT, NULL, NULL, (size_t)0)) < 0) {
        parallel_print("h5diff: <%s>: unable to open file\n", fname2);
        H5TOOLS_GOTO_ERROR(H5DIFF_ERR, "<%s>: unable to open file\n", fname2);
    } /* end if */
    H5TOOLS_DEBUG("file2_id = %s", fname2);

    /* Build object display table */
    DISPATCH(H5O_TYPE_GROUP, "Group", NULL, NULL);
    DISPATCH(H5O_TYPE_DATASET, "Dataset", dataset_list1, dataset_list2);
    DISPATCH(H5O_TYPE_NAMED_DATATYPE, "Type", NULL, datatype_list2);

    iter1 = (iter_t *)malloc(sizeof(iter_t));
    iter2 = (iter_t *)malloc(sizeof(iter_t));
    initialize_dset_info(&dsets1, iter1);
    initialize_dset_info(&dsets2, iter2);

    iter1->fname = strdup(fname1);
    iter2->fname = strdup(fname2);
    iter1->fid   = file1_id;
    iter2->fid   = file2_id;
    iter1->gid = iter2->gid = -1;
    iter1->symlink_target = iter2->symlink_target = FALSE;
    iter1->symlink_list                           = &symlink_list1;
    iter2->symlink_list                           = &symlink_list2;
    iter1->symlink_list->dangle_link = iter2->symlink_list->dangle_link = FALSE;

    iter1->base_len = iter2->base_len = 0;
    iter1->name_start = iter2->name_start = 0;

    /* Initialize list of visited symbolic links */
    symlink_list1.nused = symlink_list1.nalloc = 0;
    symlink_list2.nused = symlink_list2.nalloc = 0;
    symlink_list1.objs = symlink_list2.objs = NULL;
    oname                                   = root_name;

    if (gather_obj_info(file1_id, oname, iter1) < 0) {
        printf("gather_obj_info returned an error\n");
    }
    if (gather_obj_info(file2_id, oname, iter2) < 0) {
        printf("gather_obj_info returned an error\n");
    }

done:
    if (file1_id > 0) {
        H5Fclose(file1_id);
        iter1->fid = -1;
    }
    if (file2_id > 0) {
        H5Fclose(file2_id);
        iter2->fid = -1;
    }
    return ret_value;
}
#endif
#endif /* __NOT_USED__ */

void static diff_single(int mpi_rank, int mpi_size)
{
    printf("entered diff_single\n");
}

static diff_err_t
diff_collection(void)
{
    int         k, i = 0;
    int         mpi_rank, mpi_size;
    int         count, min_remaining;
    int         all_diffs;
    int         local_diffs;
    size_t      n_diffs          = 0;
    size_t      collective_limit = 1024 * 1024;
    hid_t       file1_id, file2_id;
    hid_t       obj1, obj2;
    hid_t       fapl = H5P_DEFAULT;
    const char *fname1, *fname2;
    diff_err_t  ret_value = H5DIFF_NO_ERR;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    fname1 = dsets1.iter->fname;
    fname2 = dsets2.iter->fname;

    /* FIXME:  Should move the file open functions inside the while
     * (elements > collective_limit) statement...
     */
    if ((file1_id = h5tools_fopen(fname1, H5F_ACC_RDONLY, H5P_DEFAULT, "mpio", NULL, (size_t)0)) < 0) {
        parallel_print("h5diff: <%s>: unable to open file\n", fname1);
        H5TOOLS_GOTO_ERROR(H5DIFF_ERR, "<%s>: unable to open file\n", fname1);
    }
    else
        dsets1.iter->fid = file1_id;

    if ((file2_id = h5tools_fopen(fname2, H5F_ACC_RDONLY, H5P_DEFAULT, "mpio", NULL, (size_t)0)) < 0) {
        parallel_print("h5diff: <%s>: unable to open file\n", fname2);
        H5TOOLS_GOTO_ERROR(H5DIFF_ERR, "<%s>: unable to open file\n", fname2);
    }
    else
        dsets2.iter->fid = file2_id;

    while (dsets1.known_dsets[i].total_elements > collective_limit) {
        int        ds_ndims;
        hsize_t    ds_dimensions[3] = {0, 0, 0};
        H5O_info_t objInfo1, objInfo2;
        hid_t      dataspace1, dataspace2;
        hsize_t    nattrs1 = 0, nattrs2;

        if ((obj1 = H5Oopen_by_token(file1_id, dsets1.known_dsets[i].token)) < 0) {
            const char *OName =
                (dsets1.known_dsets[i].name != NULL ? dsets1.known_dsets[i].name : "<unknown>");
            H5TOOLS_GOTO_ERROR(H5DIFF_ERR, "Unable to open object %s\n", OName);
        }
        if ((obj2 = H5Oopen_by_token(file2_id, dsets2.known_dsets[i].token)) < 0) {
            const char *OName =
                (dsets1.known_dsets[i].name != NULL ? dsets1.known_dsets[i].name : "<unknown>");
            H5TOOLS_GOTO_ERROR(H5DIFF_ERR, "Unable to open object %s\n", OName);
        }
        if (H5Oget_info(obj1, &objInfo1, H5O_INFO_ALL) < 0) {
            printf("H5Oget_info returned an error!\n");
        }
        if (H5Oget_info(obj2, &objInfo2, H5O_INFO_ALL) < 0) {
            printf("H5Oget_info returned an error!\n");
        }
        if ((dataspace1 = H5Dget_space(obj1)) < 0) {
            printf("H5Dget_space returned an error!\n");
        }
        if ((ds_ndims = H5Sget_simple_extent_ndims(dataspace1)) < 0) {
            printf("H5Dget_simple_extent_ndims returned an error!\n");
        }
        if (H5Sget_simple_extent_dims(dataspace1, ds_dimensions, NULL) < 0) {
            printf("H5Dget_simple_extent_dims returned an error!\n");
        }
        if (mpi_rank == 0) {
            printf("\n----- %2d ------\nName: %s\n", i, dsets1.known_dsets[i].name);
            printf("dataspace1 dims = [%lld, %lld, %lld]\n", ds_dimensions[0], ds_dimensions[1],
                   ds_dimensions[2]);
            fflush(stdout);
        }

        if (objInfo1.type == objInfo2.type) {
            hid_t  dtype1, dtype2;
            size_t type_size1, type_size2;
            if (objInfo1.type == H5O_TYPE_DATASET) {
                if (dsets1.known_dsets[i].total_elements == dsets2.known_dsets[i].total_elements) {
                    dtype1 = H5Dget_type(obj1);
                    dtype2 = H5Dget_type(obj2);
                    if ((dtype1 > 0) && (dtype2 > 0)) {
                        type_size1 = H5Tget_size(dtype1);
                        type_size2 = H5Tget_size(dtype2);
                        if (type_size1 != type_size2) {
                            printf("Dataset type sizes match (%ld != %ld)\n", type_size1, type_size2);
                            fflush(stdout);
                            n_diffs++;
                        }
                    }
                    if (compare_dims(&dsets1.known_dsets[i], &dsets2.known_dsets[i]) != 0) {
                        printf("dataset dimensions don't match\n");
                        fflush(stdout);
                        n_diffs++;
                    }
                    n_diffs +=
                        dset_compare(obj1, obj2, dsets1.known_dsets[i].ndims, dsets1.known_dsets[i].dims,
                                     dsets1.known_dsets[i].total_elements, mpi_size, mpi_rank);
                }
                else {
                    printf("dataset element counts are different\n");
                    fflush(stdout);
                    n_diffs++;
                }
            }
        }
        else {
            printf("Object types differ\n");
            fflush(stdout);
            n_diffs++;
        }
        if (mpi_rank == 0) {
            printf("opens ok, attrs = %lld, %lld\n", objInfo1.num_attrs, objInfo2.num_attrs);
            fflush(stdout);
        }
        H5Oclose(obj1);
        H5Oclose(obj2);
        i++;
    } /* End While(total_elements > collective_limit) */

    k = i + mpi_rank;

    while (k < (int)dsets1.dset_count) {
        int        ds_ndims;
        hsize_t    ds_dimensions[3] = {0, 0, 0};
        char       ds_name[MAX_NAME];
        H5O_info_t objInfo1, objInfo2;
        hid_t      dataspace1, dataspace2;
        hsize_t    nattrs1 = 0, nattrs2;

        if (mpi_rank > 0)
            memset(ds_name, 0, MAX_NAME);

        if ((obj1 = H5Oopen_by_token(file1_id, dsets1.known_dsets[k].token)) < 0) {
            const char *OName =
                (dsets1.known_dsets[k].name != NULL ? dsets1.known_dsets[k].name : "<unknown>");
            H5TOOLS_GOTO_ERROR(H5DIFF_ERR, "Unable to open object %s\n", OName);
        }
        if ((obj2 = H5Oopen_by_token(file2_id, dsets2.known_dsets[k].token)) < 0) {
            const char *OName =
                (dsets1.known_dsets[k].name != NULL ? dsets1.known_dsets[k].name : "<unknown>");
            H5TOOLS_GOTO_ERROR(H5DIFF_ERR, "Unable to open object %s\n", OName);
        }
        if (mpi_rank > 0) {
            ssize_t len = H5Iget_name(obj1, ds_name, MAX_NAME);
            if (len > 0) {
                dsets1.known_dsets[k].name = strdup(ds_name);
            }
        }
        if (dsets2.known_dsets[k].name && dsets2.known_dsets[k].name[0] == 0) {
            ssize_t len = H5Fget_name(obj2, ds_name, MAX_NAME);
            if (len > 0) {
                dsets2.known_dsets[k].name = strdup(ds_name);
            }
        }
        if (H5Oget_info(obj1, &objInfo1, H5O_INFO_ALL) < 0) {
            printf("H5Oget_info returned an error!\n");
            fflush(stdout);
        }
        if (H5Oget_info(obj2, &objInfo2, H5O_INFO_ALL) < 0) {
            printf("H5Oget_info returned an error!\n");
            fflush(stdout);
        }
        if ((dataspace1 = H5Dget_space(obj1)) < 0) {
            printf("H5Dget_space returned an error!\n");
            fflush(stdout);
        }
        if ((ds_ndims = H5Sget_simple_extent_ndims(dataspace1)) < 0) {
            printf("H5Dget_simple_extent_ndims returned an error!\n");
            fflush(stdout);
        }
        if (H5Sget_simple_extent_dims(dataspace1, ds_dimensions, NULL) < 0) {
            printf("H5Dget_simple_extent_dims returned an error!\n");
            fflush(stdout);
        }
        printf("\n----- %2d ------\nName: %s\n", k, dsets1.known_dsets[k].name);
        printf("dataspace1 dims = [%lld, %lld, %lld]\n", ds_dimensions[0], ds_dimensions[1],
               ds_dimensions[2]);
        fflush(stdout);

        if (objInfo1.type == objInfo2.type) {
            hid_t  dtype1, dtype2;
            size_t type_size1, type_size2;
            if (objInfo1.type == H5O_TYPE_DATASET) {
                if (dsets1.known_dsets[i].total_elements == dsets2.known_dsets[i].total_elements) {
                    dtype1 = H5Dget_type(obj1);
                    dtype2 = H5Dget_type(obj2);
                    if ((dtype1 > 0) && (dtype2 > 0)) {
                        type_size1 = H5Tget_size(dtype1);
                        type_size2 = H5Tget_size(dtype2);
                        if (type_size1 != type_size2) {
                            printf("Dataset type sizes match (%ld != %ld)\n", type_size1, type_size2);
                            fflush(stdout);
                            n_diffs++;
                        }
                    }
                    if (compare_dims(&dsets1.known_dsets[i], &dsets2.known_dsets[i]) != 0) {
                        printf("dataset dimensions don't match\n");
                        fflush(stdout);
                        n_diffs++;
                    }
                    n_diffs +=
                        dset_compare_seq(obj1, obj2, dsets1.known_dsets[i].ndims, dsets1.known_dsets[i].dims,
                                         dsets1.known_dsets[i].total_elements);
                    n_diffs += dset_attr_compare(obj1, obj2, &objInfo1, &objInfo2);
                }
                else {
                    printf("dataset element counts are different\n");
                    fflush(stdout);
                    n_diffs++;
                }
            }
        }
        else {
            printf("Object types differ\n");
            fflush(stdout);
            n_diffs++;
        }
        if (mpi_rank == 0) {
            printf("opens ok, attrs = %ld, %ld\n", objInfo1.num_attrs, objInfo2.num_attrs);
            fflush(stdout);
        }
        H5Oclose(obj1);
        H5Oclose(obj2);
        k += mpi_size;
        usleep(0);
    }
    MPI_Barrier(MPI_COMM_WORLD);

done:
    local_diffs = (int)n_diffs;
    MPI_Allreduce(&local_diffs, &all_diffs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (file1_id > 0)
        H5Fclose(file1_id);
    if (file2_id > 0)
        H5Fclose(file2_id);

#if 0
    printf("[%d] Exiting %s\n", mpi_rank, __FUNCTION__);
    fflush(stdout);
#endif

    if (all_diffs != 0)
        ret_value = H5DIFF_ERR_DIFF;
    return ret_value;
}

/*-------------------------------------------------------------------------
 * Function: main
 *
 * Purpose: ph5diff main program
 *
 * Return: An exit status of 0 means no differences were found, 1 means some
 *   differences were found.
 *
 * Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
 *
 * Date: May 9, 2003
 *
 * Comments:
 *
 * This function drives the diff process and will do a serial or parallel diff depending
 * on the value of the global variable g_Parallel (default is 0), set to 1 when the program
 * is run as "ph5diff"
 *-------------------------------------------------------------------------
 */

int
main(int argc, const char *argv[])
{
    int         nID       = 0;
    const char *fname1    = NULL;
    const char *fname2    = NULL;
    const char *objname1  = NULL;
    const char *objname2  = NULL;
    diff_err_t  ret_value = H5DIFF_NO_ERR;
    diff_opt_t  opts;

    h5tools_setprogname(PROGRAMNAME);
    h5tools_setstatus(EXIT_SUCCESS);

    /* Initialize h5tools lib */
    h5tools_init();

    outBuffOffset = 0;
    g_Parallel    = 1;

    MPI_Init(&argc, (char ***)&argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &nID);
    MPI_Comm_size(MPI_COMM_WORLD, &g_nTasks);

    parse_h5diff_command_line(argc, argv, &fname1, &fname2, &objname1, &objname2, &opts);

    if (g_nTasks == 1) {
        HDprintf("Only 1 task available...doing serial diff\n");

        g_Parallel = 0;
        h5diff(fname1, fname2, objname1, objname2, &opts);
        print_info(&opts);
    }
    else {
        g_Parallel = 1;
        size_t ret = ph5diff(fname1, fname2, objname1, objname2, &opts);
        // print_info(&opts);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    exit(EXIT_SUCCESS);
    return 0;
}

/*-------------------------------------------------------------------------
 * Function: ph5diff_worker
 *
 * Purpose: worker process of ph5diff
 *
 * Return: none
 *
 * Programmer: Leon Arber
 * Date: January 2005
 *
 * Comments:
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static void
ph5diff_worker(int nID)
{
    hid_t file1_id = H5I_INVALID_HID;
    hid_t file2_id = H5I_INVALID_HID;

    while (1) {
        MPI_Status Status;

        MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);

        /* Check for filenames */
        if (Status.MPI_TAG == MPI_TAG_PARALLEL) {
            char filenames[2][MAX_FILENAME];

            /* Retrieve filenames */
            MPI_Recv(filenames, MAX_FILENAME * 2, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);

            /* disable error reporting */
            H5E_BEGIN_TRY
            {
                /* Open the files */
                if ((file1_id = H5Fopen(filenames[0], H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
                    HDprintf("h5diff Task [%d]: <%s>: unable to open file\n", nID, filenames[0]);
                    MPI_Abort(MPI_COMM_WORLD, 0);
                }
                if ((file2_id = H5Fopen(filenames[1], H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
                    HDprintf("h5diff Task [%d]: <%s>: unable to open file\n", nID, filenames[1]);
                    MPI_Abort(MPI_COMM_WORLD, 0);
                }
                /* enable error reporting */
            }
            H5E_END_TRY;
        }
        /* Check for work */
        else if (Status.MPI_TAG == MPI_TAG_ARGS) {
            struct diff_mpi_args args;
            struct diffs_found   diffs;
            int                  i;

            /* Make certain we've received the filenames and opened the files already */
            if (file1_id < 0 || file2_id < 0) {
                HDprintf("ph5diff_worker: ERROR: work received before/without filenames\n");
                break;
            }

            /* Recv parameters for diff from manager task */
            MPI_Recv(&args, sizeof(args), MPI_BYTE, 0, MPI_TAG_ARGS, MPI_COMM_WORLD, &Status);

            /* Do the diff */
            diffs.nfound  = diff(file1_id, args.name1, file2_id, args.name2, &(args.opts), &(args.argdata));
            diffs.not_cmp = args.opts.not_cmp;

            /* If print buffer has something in it, request print token.*/
            if (outBuffOffset > 0) {
                MPI_Send(NULL, 0, MPI_BYTE, 0, MPI_TAG_TOK_REQUEST, MPI_COMM_WORLD);

                /* Wait for print token. */
                MPI_Recv(NULL, 0, MPI_BYTE, 0, MPI_TAG_PRINT_TOK, MPI_COMM_WORLD, &Status);

                /* When get token, send all of our output to the manager task and then return the token */
                for (i = 0; i < outBuffOffset; i += PRINT_DATA_MAX_SIZE)
                    MPI_Send(outBuff + i, PRINT_DATA_MAX_SIZE, MPI_BYTE, 0, MPI_TAG_PRINT_DATA,
                             MPI_COMM_WORLD);

                /* An overflow file exists, so we send it's output to the manager too and then delete it */
                if (overflow_file) {
                    char out_data[PRINT_DATA_MAX_SIZE];
                    int  tmp;

                    HDmemset(out_data, 0, PRINT_DATA_MAX_SIZE);
                    i = 0;

                    rewind(overflow_file);
                    while ((tmp = getc(overflow_file)) >= 0) {
                        *(out_data + i++) = (char)tmp;
                        if (i == PRINT_DATA_MAX_SIZE) {
                            MPI_Send(out_data, PRINT_DATA_MAX_SIZE, MPI_BYTE, 0, MPI_TAG_PRINT_DATA,
                                     MPI_COMM_WORLD);
                            i = 0;
                            HDmemset(out_data, 0, PRINT_DATA_MAX_SIZE);
                        }
                    }

                    if (i > 0)
                        MPI_Send(out_data, PRINT_DATA_MAX_SIZE, MPI_BYTE, 0, MPI_TAG_PRINT_DATA,
                                 MPI_COMM_WORLD);

                    fclose(overflow_file);
                    overflow_file = NULL;
                }

                HDfflush(stdout);
                HDmemset(outBuff, 0, OUTBUFF_SIZE);
                outBuffOffset = 0;

                MPI_Send(&diffs, sizeof(diffs), MPI_BYTE, 0, MPI_TAG_TOK_RETURN, MPI_COMM_WORLD);
            }
            else
                MPI_Send(&diffs, sizeof(diffs), MPI_BYTE, 0, MPI_TAG_DONE, MPI_COMM_WORLD);
        }
        /* Check for leaving */
        else if (Status.MPI_TAG == MPI_TAG_END) {
            MPI_Recv(NULL, 0, MPI_BYTE, 0, MPI_TAG_END, MPI_COMM_WORLD, &Status);
            break;
        }
        else {
            HDprintf("ph5diff_worker: ERROR: invalid tag (%d) received\n", Status.MPI_TAG);
            break;
        }
    }

    return;
}

/*-------------------------------------------------------------------------
 * Function: print_manager_output
 *
 * Purpose: special function that prints any output accumulated by the
 *      manager task.
 *
 * Return: none
 *
 * Programmer: Leon Arber
 *
 * Date: Feb 7, 2005
 *
 *-------------------------------------------------------------------------
 */
void
print_manager_output(void)
{
    /* If there was something we buffered, let's print it now */
    if ((outBuffOffset > 0) && g_Parallel) {
        HDprintf("%s", outBuff);

        if (overflow_file) {
            int tmp;
            rewind(overflow_file);
            while ((tmp = HDgetc(overflow_file)) >= 0)
                HDputchar(tmp);
            fclose(overflow_file);
            overflow_file = NULL;
        }

        HDfflush(stdout);
        HDmemset(outBuff, 0, OUTBUFF_SIZE);
        outBuffOffset = 0;
    }
    else if ((outBuffOffset > 0) && !g_Parallel) {
        HDfprintf(stderr, "h5diff error: outBuffOffset>0, but we're not in parallel!\n");
    }
}

/*-------------------------------------------------------------------------
 * Function: h5diff_exit
 *
 * Purpose: dismiss phdiff worker processes and exit
 *
 * Return: none
 *
 * Programmer: Albert Cheng
 * Date: Feb 6, 2005
 *
 * Comments:
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
void
ph5diff_exit(int status)
{
    /* if in parallel mode, dismiss workers, close down MPI, then exit */
    if (g_Parallel) {
        if (g_nTasks > 1) {
            phdiff_dismiss_workers();
            MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Finalize();
        status = EXIT_SUCCESS; /* Reset exit status, since some mpiexec commands generate output on failure
                                  status */
    }

    h5tools_close();

    /* Always exit(0), since MPI implementations do weird stuff when they
     *  receive a non-zero exit value. - QAK
     */
    HDexit(0);
}
