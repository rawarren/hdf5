/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the COPYING file, which can be found at the root of the source code       *
 * distribution tree, or in https://www.hdfgroup.org/licenses.               *
 * If you do not have access to either file, you may request a copy from     *
 * help@hdfgroup.org.                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * Programmer:  Bill Wendling
 *              Tuesday, 6. March 2001
 *
 * Purpose:     Support functions for the various tools.
 */
#ifndef H5TOOLS_UTILS_H
#define H5TOOLS_UTILS_H

#include "hdf5.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ``parallel_print'' information */
#define PRINT_DATA_MAX_SIZE 512
#define OUTBUFF_SIZE        (PRINT_DATA_MAX_SIZE * 4)

H5TOOLS_DLLVAR int           g_nTasks;
H5TOOLS_DLLVAR unsigned char g_Parallel;
H5TOOLS_DLLVAR char          outBuff[];
H5TOOLS_DLLVAR unsigned      outBuffOffset;
H5TOOLS_DLLVAR FILE *overflow_file;

/* Maximum size used in a call to malloc for a dataset */
H5TOOLS_DLLVAR hsize_t H5TOOLS_MALLOCSIZE;
/* size of hyperslab buffer when a dataset is bigger than H5TOOLS_MALLOCSIZE */
H5TOOLS_DLLVAR hsize_t H5TOOLS_BUFSIZE;

/*struct taken from the dumper. needed in table struct*/

typedef struct dspace_table_t {
    hid_t       fid;
    hid_t       space_id;
	hsize_t     size[H5TOOLS_DUMP_MAX_RANK];
    hsize_t     maxsize[H5TOOLS_DUMP_MAX_RANK];
	int         ndims;
} ds_table_t;

typedef struct obj_t {
    H5O_token_t obj_token;
    char *      objname;
    hbool_t     displayed; /* Flag to indicate that the object has been displayed */
    hbool_t     recorded;  /* Flag for named datatypes to indicate they were found in the group hierarchy */
} obj_t;

/*struct for the tables that the find_objs function uses*/
typedef struct table_t {
    hid_t  fid;
    size_t size;
    size_t nobjs;
    obj_t *objs;
} table_t;

/* List of table structures.  There is one table structure for each file */
typedef struct h5tools_table_items_t {
    unsigned long fileno;      /* File number that these tables refer to */
    hid_t         oid;         /* ID of an object in this file, held open so fileno is consistent */
    table_t *     group_table; /* Table of groups */
    table_t *     dset_table;  /* Table of datasets */
    table_t *     type_table;  /* Table of datatypes */
    ds_table_t *   ds_table;    /* Table of dataset dataspace dimensions (order aligned) */
} h5tools_table_items_t;

typedef struct {
    size_t                nalloc;
    size_t                nused;
    h5tools_table_items_t *tables;
} h5tools_table_list_t;

extern h5tools_table_list_t h5tools_table_list;
extern table_t *h5dump_type_table;

/*this struct stores the information that is passed to the find_objs function*/
typedef struct find_objs_t {
    hid_t    fid;
    table_t *group_table;
    table_t *type_table;
    table_t *dset_table;
    ds_table_t *dspace_table;
} find_objs_t;

H5TOOLS_DLLVAR table_t *h5tools_group_table;
H5TOOLS_DLLVAR table_t *h5tools_dset_table;
H5TOOLS_DLLVAR table_t *h5tools_type_table;
H5TOOLS_DLLVAR unsigned h5tools_nCols; /*max number of columns for outputting  */

/* Definitions of useful routines */
H5TOOLS_DLL void   indentation(unsigned);
H5TOOLS_DLL void   print_version(const char *progname);
H5TOOLS_DLL void   parallel_print(const char *format, ...);
H5TOOLS_DLL herr_t parse_tuple(const char *start, int sep, char **cpy_out, unsigned *nelems,
                               char ***ptrs_out);
H5TOOLS_DLL void   error_msg(const char *fmt, ...);
H5TOOLS_DLL void   warn_msg(const char *fmt, ...);
H5TOOLS_DLL void   help_ref_msg(FILE *output);
H5TOOLS_DLL void   free_table(table_t *table);
H5TOOLS_DLL void   h5tools_init_prefix(char **prfx, size_t prfx_len);
H5TOOLS_DLL void   h5tools_table_list_free(void);
H5TOOLS_DLL ssize_t h5tools_table_list_visited(unsigned long file_no);
H5TOOLS_DLL ssize_t h5tools_table_list_add(hid_t oid, unsigned long file_no);



#ifdef H5DUMP_DEBUG
H5TOOLS_DLL void dump_tables(find_objs_t *info);
#endif /* H5DUMP_DEBUG */
H5TOOLS_DLL herr_t init_objs(hid_t fid, find_objs_t *info, table_t **group_table, table_t **dset_table,
                             table_t **type_table, ds_table_t **dspace_table);
H5TOOLS_DLL obj_t *search_obj(table_t *temp, const H5O_token_t *obj_token);
#ifndef H5_HAVE_TMPFILE
H5TOOLS_DLL FILE *tmpfile(void);
#endif

/*************************************************************
 *
 * candidate functions to be public
 *
 *************************************************************/

/* This code is layout for common code among tools */
typedef enum toolname_t {
    TOOL_H5DIFF,
    TOOL_H5LS,
    TOOL__H5DUMP /* add as necessary */
} h5tool_toolname_t;

/* this struct can be used to differntiate among tools */
typedef struct {
    h5tool_toolname_t toolname;
    int               msg_mode;
} h5tool_opt_t;

/* obtain link info from H5tools_get_symlink_info() */
typedef struct {
    H5O_type_t trg_type;     /* OUT: target type */
    char *     trg_path;     /* OUT: target obj path. This must be freed
                              *      when used with H5tools_get_symlink_info() */
    H5O_token_t   obj_token; /* OUT: target object token */
    unsigned long fileno;    /* OUT: File number that target object is located in */
    H5L_info2_t   linfo;     /* OUT: link info */
    h5tool_opt_t  opt;       /* IN: options */
} h5tool_link_info_t;

/* Definitions of routines */
H5TOOLS_DLL int H5tools_get_symlink_info(hid_t file_id, const char *linkpath, h5tool_link_info_t *link_info,
                                         hbool_t get_obj_type);
H5TOOLS_DLL const char *h5tools_getprogname(void);
H5TOOLS_DLL void        h5tools_setprogname(const char *progname);
H5TOOLS_DLL int         h5tools_getstatus(void);
H5TOOLS_DLL void        h5tools_setstatus(int d_status);
H5TOOLS_DLL int         h5tools_getenv_update_hyperslab_bufsize(void);

H5TOOLS_DLL void        h5tools_handle_paths(hid_t fid, const char *path_name, void H5_ATTR_UNUSED *data,
                                             int H5_ATTR_UNUSED pe, const char H5_ATTR_UNUSED *display_name);
H5TOOLS_DLL void        h5tools_handle_attributes(hid_t fid, const char *attr, void H5_ATTR_UNUSED *data,
                                             int H5_ATTR_UNUSED pe, const char H5_ATTR_UNUSED *display_name);
H5TOOLS_DLL void        h5tools_handle_datasets(hid_t fid, const char *dset, void *data,
                                             int H5_ATTR_UNUSED pe, const char H5_ATTR_UNUSED *display_name);
H5TOOLS_DLL void        h5tools_handle_groups(hid_t fid, const char *group, void H5_ATTR_UNUSED *data,
                                             int H5_ATTR_UNUSED pe, const char H5_ATTR_UNUSED *display_name);
H5TOOLS_DLL void        h5tools_handle_links(hid_t fid, const char *links, void H5_ATTR_UNUSED *data,
                                             int H5_ATTR_UNUSED pe, const char H5_ATTR_UNUSED *display_name);
H5TOOLS_DLL void        h5tools_handle_datatypes(hid_t fid, const char *type, void H5_ATTR_UNUSED *data,
                                             int H5_ATTR_UNUSED pe, const char H5_ATTR_UNUSED *display_name);

#ifdef H5_HAVE_ROS3_VFD
H5TOOLS_DLL herr_t h5tools_parse_ros3_fapl_tuple(const char *tuple_str, int delim,
                                                 H5FD_ros3_fapl_t *fapl_config_out);
H5TOOLS_DLL int    h5tools_populate_ros3_fapl(H5FD_ros3_fapl_t *fa, const char **values);
#endif /* H5_HAVE_ROS3_VFD */
#ifdef H5_HAVE_LIBHDFS
H5TOOLS_DLL herr_t h5tools_parse_hdfs_fapl_tuple(const char *tuple_str, int delim,
                                                 H5FD_hdfs_fapl_t *fapl_config_out);
#endif

#ifdef __cplusplus
}
#endif

#endif /* H5TOOLS_UTILS_H */
