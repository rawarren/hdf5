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

#include <string.h>
#include "testpar.h"
#include "mpi.h"

#define KB 1024U
#define VECTOR_LEN 16

#define NFILENAME 4
const char *FILENAME[NFILENAME + 1]=
{
    "subfiling_file",
    "subfile_read_write_bench",
    "subfile_multifile_0_file",
    "reloc_t_multifile_1_file",
    NULL
};

#define FILENAME_BUF_SIZE 1024

#define DATASETNAME1    "Data1"
#define RANK            2
#define ROW_FACTOR      200     /* Nominal row factor for dataset size */
#define COL_FACTOR      5000    /* Nominal column factor for dataset size */
#define NUM_DSETS       5
#define COUNT           ROW_FACTOR * COL_FACTOR


/* File_Access_type bits */
#define FACC_DEFAULT    0x0     /* default */
#define FACC_MPIO       0x1     /* MPIO */
#define FACC_SPLIT      0x2     /* Split File */
#define FACC_SUBFILE    0x4     /* Subfiling */

int nerrors = 0;
int mpi_size;
int mpi_rank;


/*
 * Create the appropriate File access property list
 */
static hid_t
create_faccess_plist(MPI_Comm comm, MPI_Info info, int l_facc_type)
{
    hid_t ret_pl = -1;
    herr_t ret;                 /* generic return value */

    /* need the rank for error checking macros */
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    ret_pl = H5Pcreate (H5P_FILE_ACCESS);

    if (l_facc_type == FACC_DEFAULT)
        return (ret_pl);

    if (l_facc_type == FACC_MPIO){
        /* set Parallel access with communicator */
        ret = H5Pset_fapl_mpio(ret_pl, comm, info);
        ret = H5Pset_all_coll_metadata_ops(ret_pl, TRUE);
        ret = H5Pset_coll_metadata_write(ret_pl, TRUE);
        return(ret_pl);
    }

    if (l_facc_type == (FACC_MPIO | FACC_SPLIT)){
        hid_t mpio_pl;

        mpio_pl = H5Pcreate (H5P_FILE_ACCESS);
        /* set Parallel access with communicator */
        ret = H5Pset_fapl_mpio(mpio_pl, comm, info);

        /* setup file access template */
        ret_pl = H5Pcreate (H5P_FILE_ACCESS);

        /* set Parallel access with communicator */
        ret = H5Pset_fapl_split(ret_pl, ".meta", mpio_pl, ".raw", mpio_pl);
        H5Pclose(mpio_pl);
        return(ret_pl);
    }
    if (l_facc_type == FACC_SUBFILE ){
        H5FD_subfiling_fapl_t fa_in = {H5FD_CURR_SUBFILING_FAPL_T_VERSION};
        ret = H5Pset_fapl_subfiling(ret_pl, &fa_in);
        return (ret_pl);
    }
    if (l_facc_type == (FACC_SUBFILE | FACC_MPIO)){
        H5FD_subfiling_fapl_t fa_in = {H5FD_CURR_SUBFILING_FAPL_T_VERSION};
        ret = H5Pset_fapl_subfiling(ret_pl, &fa_in);
        ret = H5Pset_all_coll_metadata_ops(ret_pl, TRUE);
        ret = H5Pset_coll_metadata_write(ret_pl, TRUE);
    }
    /* unknown file access types */
    return (ret_pl);
}


static herr_t
test_subfiling(void)
{
    hid_t       fid = -1;                   /* file ID                      */
    hid_t       fapl_id = -1;               /* file access property list ID */
    hid_t       fapl_id_out = -1;           /* from H5Fget_access_plist     */
    hid_t       driver_id = -1;             /* ID for this VFD              */
    unsigned long driver_flags = 0;         /* VFD feature flags            */
    char        filename[1024];             /* filename                     */
    void        *os_file_handle = NULL;     /* OS file handle               */
    hsize_t     file_size;                  /* file size                    */
    H5FD_subfiling_fapl_t fa_in = {H5FD_CURR_SUBFILING_FAPL_T_VERSION};
    H5FD_subfiling_fapl_t fa_out;

    if (mpi_rank == 0) {
        TESTING("subfiling file driver");
    }
    /* Set property list and file name for subfiling driver. */
    if((fapl_id = H5Pcreate(H5P_FILE_ACCESS)) < 0) {
        if (mpi_rank == 0) {
            TEST_ERROR;
        }
        goto error;
    }       

    if(H5Pset_fapl_subfiling(fapl_id, &fa_in) < 0) {
        if (mpi_rank == 0) {
            TEST_ERROR;
        }
        goto error;
    }       

    /* get and verify the H5FD_subfiling_fapl_t */
    if(H5Pget_fapl_subfiling(fapl_id, &fa_out) < 0) {
        if (mpi_rank == 0) {
            TEST_ERROR;
        }
        goto error;
    }

    if(fa_out.version != H5FD_CURR_SUBFILING_FAPL_T_VERSION) {
        if (mpi_rank == 0) {
            TEST_ERROR;
        }
        goto error;
    }

    h5_fixname(FILENAME[0], fapl_id, filename, sizeof(filename));

    /* Check that the VFD feature flags are correct */
    if ((driver_id = H5Pget_driver(fapl_id)) < 0) {
        if (mpi_rank == 0) {
            TEST_ERROR;
        }
        goto error;
    }

    if (H5FDdriver_query(driver_id, &driver_flags) < 0) {
        if (mpi_rank == 0) {
            TEST_ERROR;
        }
        goto error;
    }

    if(!(driver_flags & H5FD_FEAT_AGGREGATE_METADATA))      TEST_ERROR
    if(!(driver_flags & H5FD_FEAT_ACCUMULATE_METADATA))     TEST_ERROR
    if(!(driver_flags & H5FD_FEAT_DATA_SIEVE))              TEST_ERROR
    if(!(driver_flags & H5FD_FEAT_AGGREGATE_SMALLDATA))     TEST_ERROR
    if(!(driver_flags & H5FD_FEAT_POSIX_COMPAT_HANDLE))     TEST_ERROR
    if(!(driver_flags & H5FD_FEAT_SUPPORTS_SWMR_IO))        TEST_ERROR

    if((fid = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id)) < 0) {
        if (mpi_rank == 0) {
            TEST_ERROR;
        }
        goto error;
    }

    /* Retrieve the access property list... */
    if((fapl_id_out = H5Fget_access_plist(fid)) < 0) {
        if (mpi_rank == 0) {
            TEST_ERROR;
        }
        goto error;
    }

    /* Check that the driver is correct */
    if(H5FD_SUBFILING != H5Pget_driver(fapl_id_out)) {
        if (mpi_rank == 0) {
            TEST_ERROR;
        }
        goto error;
    }

    /* get and verify the H5FD_subfiling_fapl_t again */
    if(H5Pget_fapl_subfiling(fapl_id_out, &fa_out) < 0) {
        if (mpi_rank == 0) {
            TEST_ERROR;
        }
        goto error;
    }

    if(fa_out.version != H5FD_CURR_SUBFILING_FAPL_T_VERSION) {
        if (mpi_rank == 0) {
            TEST_ERROR;
        }
        goto error;
    }

    /* ...and close the property list */
    if(H5Pclose(fapl_id_out) < 0) {
        if (mpi_rank == 0) {
            TEST_ERROR;
        }
        goto error;
    }

    /* Check that we can get an operating-system-specific handle from
     * the library.
     *
     * Not sure that this will be meaningful in the subfiling case.
     */
    if(H5Fget_vfd_handle(fid, H5P_DEFAULT, &os_file_handle) < 0) {
        if (mpi_rank == 0) {
            TEST_ERROR;
        }
        goto error;
    }

    if(os_file_handle == NULL)
        FAIL_PUTS_ERROR("NULL os-specific vfd/file handle was returned from H5Fget_vfd_handle");


    /* There is no garantee the size of metadata in file is constant.
     * Just try to check if it's reasonable.
     *
     * Currently it should be around 2 KB.
     */
    if(H5Fget_filesize(fid, &file_size) < 0) {
        if (mpi_rank == 0) {
            TEST_ERROR;
        }
        goto error;
    }


    if(file_size < 1 * KB || file_size > 4 * KB)
        FAIL_PUTS_ERROR("suspicious file size obtained from H5Fget_filesize");

    /* Close and delete the file */
    if(H5Fclose(fid) < 0) {
        if (mpi_rank == 0) {
            TEST_ERROR;
        }
        goto error;
    }

    h5_delete_test_file(FILENAME[0], fapl_id);

    /* Close the fapl */
    if(H5Pclose(fapl_id) < 0) {
        if (mpi_rank == 0) {
            TEST_ERROR;
        }
        goto error;
    }

    if (mpi_rank == 0) {
        PASSED();
    }
    return 0;

error:
    H5E_BEGIN_TRY {
        H5Pclose(fapl_id);
        H5Pclose(fapl_id_out);
        H5Fclose(fid);
    } H5E_END_TRY;

    return -1;

} /* end test_subfiling() */



/*-------------------------------------------------------------------------
 * Function:    test_vector_io__setup_v
 *
 * Purpose:     Construct and initialize a vector of I/O requests used 
 *              to test vector I/O.  Note that while the vectors are 
 *              allocated and initialized, they are not assigned 
 *              base addresses.
 *
 *              All arrays parameters are presumed to be of length 
 *              count.
 *
 * Return:      Return TRUE if sucessful, and FALSE if any errors 
 *              are encountered.
 *
 * Programmer:  John Mainzer
 *              6/21/20
 *
 * Modifications:
 *
 *              None.
 *
 *-------------------------------------------------------------------------
 */

static hbool_t
test_vector_io__setup_v(int64_t elements_per_rank, int32_t count, H5FD_mem_t types[],
                        haddr_t addrs[], size_t sizes[], void * write_bufs[],
                        void * read_bufs[], char base_fill_char)
{
    int32_t i;
    int32_t j;
    int32_t last_count;
    int64_t vector_elements = elements_per_rank / VECTOR_LEN;
    int64_t last_vector = elements_per_rank % VECTOR_LEN;
    hbool_t result = TRUE; /* will set to FALSE on failure */
    char fill_char = base_fill_char;
    H5FD_mem_t mem_types[6] = {H5FD_MEM_SUPER, H5FD_MEM_BTREE, H5FD_MEM_DRAW, 
                               H5FD_MEM_GHEAP, H5FD_MEM_LHEAP, H5FD_MEM_OHDR};

    /* set the arrays of pointers to the write and read buffers to NULL,
     * so that we can release memory on failure.
     */
    for ( i = 0; i < count; i++ ) {
        write_bufs[i] = NULL;
        read_bufs[i] = NULL;
    }

    last_count = ((last_vector == 0) ? count : count - 1);
    for ( i = 0; i < count; i++ ) {

        types[i] = mem_types[i % 6];

        addrs[i] = HADDR_UNDEF;

        sizes[i] = (size_t)vector_elements;

        write_bufs[i] = HDmalloc(sizes[i] + 1);
        read_bufs[i] = HDmalloc(sizes[i] + 1);

        if ( ( NULL == write_bufs[i] ) || ( NULL == read_bufs[i] ) ) {
            HDfprintf(stderr, "%s: can't malloc read / write bufs.\n", FUNC);
            result = FALSE;
            goto done;
        }

        for ( j = 0; j < sizes[i]; j++ ) {
            ((char *)(write_bufs[i]))[j] = fill_char;
            ((char *)(read_bufs[i]))[j] = '\0';
        }

        ((char *)(write_bufs[i]))[sizes[i]] = '\0';
        ((char *)(read_bufs[i]))[sizes[i]] = '\0';

        fill_char++;
    }

    /* The last vector (if there is one) has the leftover elements */
    if (last_vector) {
        types[i] = mem_types[i % 6];

        addrs[i] = HADDR_UNDEF;

        sizes[i] = (size_t)last_vector;

        write_bufs[i] = HDmalloc(sizes[i] + 1);
        read_bufs[i] = HDmalloc(sizes[i] + 1);
        if ( ( NULL == write_bufs[i] ) || ( NULL == read_bufs[i] ) ) {
            HDfprintf(stderr, "%s: can't malloc read / write bufs.\n", FUNC);
            result = FALSE;
            goto done;
        }

        for ( j = 0; j < sizes[i]; j++ ) {
            ((char *)(write_bufs[i]))[j] = fill_char;
            ((char *)(read_bufs[i]))[j] = '\0';
        }

        ((char *)(write_bufs[i]))[sizes[i]] = '\0';
        ((char *)(read_bufs[i]))[sizes[i]] = '\0';

        fill_char++;
    }

done:
    if ( ! result ) { /* free buffers */

        for ( i = 0; i < count; i++ ) {

            if ( write_bufs[i] ) {

                HDfree(write_bufs[i]);
                write_bufs[i] = NULL;
            }

            if ( read_bufs[i] ) {

                HDfree(read_bufs[i]);
                read_bufs[i] = NULL;
            }
        }
    }

    return(result);

} /* end test_vector_io__setup_v() */


static int
vector_io_bench(char *filename, int64_t elements_per_rank)
{
    int errors = 0;
    unsigned flags = 0;                  /* file open flags              */
    H5FD_t *lf = NULL;                   /* VFD struct ptr               */
    uint32_t   i;                        /* index                        */
    uint32_t   j;                        /* index                        */
    uint32_t   count = VECTOR_LEN;       /* length of vectors            */
    H5FD_mem_t types_0[VECTOR_LEN];      /* types vector                 */
    H5FD_mem_t types_1[VECTOR_LEN];      /* types vector                 */
    H5FD_mem_t types_2[VECTOR_LEN];      /* types vector                 */
    haddr_t    addrs_0[VECTOR_LEN];      /* addresses vector             */
    haddr_t    addrs_1[VECTOR_LEN];      /* addresses vector             */
    haddr_t    addrs_2[VECTOR_LEN];      /* addresses vector             */
    size_t     sizes_0[VECTOR_LEN];      /* sizes vector                 */
    size_t     sizes_1[VECTOR_LEN];      /* sizes vector                 */
    size_t     sizes_2[VECTOR_LEN];      /* sizes vector                 */
    void *     write_bufs_0[VECTOR_LEN]; /* write bufs vector            */
    void *     write_bufs_1[VECTOR_LEN]; /* write bufs vector            */
    void *     write_bufs_2[VECTOR_LEN]; /* write bufs vector            */
    void *     read_bufs_0[VECTOR_LEN];  /* read bufs vector             */
    void *     read_bufs_1[VECTOR_LEN];  /* read bufs vector             */
    void *     read_bufs_2[VECTOR_LEN];  /* read bufs vector             */

    double     t_start, t_end;
    double     t_total_write = 0.0;
    double     t_total_read = 0.0;
    double     bw = 0.0;
    size_t     total_size = 0;
    haddr_t    addr_offset = 0;

    H5FD_subfiling_fapl_t fa_in = {H5FD_CURR_SUBFILING_FAPL_T_VERSION};
    hid_t fapl_id;

    if((fapl_id = H5Pcreate(H5P_FILE_ACCESS)) < 0)
        TEST_ERROR;

    if(H5Pset_fapl_subfiling(fapl_id, &fa_in) < 0)
        TEST_ERROR;

    /* setup the test vectors -- note that addresses are not set until 
     * we allocate space via the file driver.
     */
    if ( ! ( test_vector_io__setup_v(elements_per_rank, count, types_0,
                 addrs_0, sizes_0, write_bufs_0, read_bufs_0, 'a'))) {
        if (mpi_rank == 0) {
             TEST_ERROR;
        }
        else goto error;
    }

    flags = H5F_ACC_RDWR | H5F_ACC_CREAT | H5F_ACC_TRUNC;

    if ( NULL == (lf = H5FDopen(filename, flags, fapl_id, HADDR_UNDEF)))
        TEST_ERROR;
    
    addr_offset = elements_per_rank * mpi_rank;
    for ( i = 0; i < count; i++ ) {
        addrs_0[i] = H5FDalloc(lf, types_0[i], H5P_DEFAULT, (hsize_t)(sizes_0[i]));
        addrs_0[i] += addr_offset; /* Update the logical file offsets */
    }

    /* For this particular test, we use a simple approach which 
     * takes a contiguous block of storage and breaking the IO into
     * 'count' vectors. The size of the local contiguous block is 'elements_per_rank'
     * and each MPI rank will start at a multiple of this.  The global contiguous
     * block is 'total_size'. Each vector element is described by an entry in the 
     * 'types_0[]','addrs_0[]', sizes_0[], and 'write_bufs_0[]' or 'read_bufs_0[]'
     * vector variables.
     */
    total_size = elements_per_rank * mpi_size;

    /* Make the file EOA the same for all parallel ranks */
    addr_offset = lf->base_addr + total_size;
    if ( H5FDset_eoa(lf, types_0[0], addr_offset) != SUCCEED)
        TEST_ERROR;

    for(i = 0; i < 10; i++) {
        /* Time vector writes for the benchmark */
        t_start = MPI_Wtime();
        if ( H5FDwrite_vector(lf, H5P_DEFAULT, count, &(types_0[0]), &(addrs_0[0]), 
                          &(sizes_0[0]), &(write_bufs_0[0])) < 0 )
            TEST_ERROR;

        t_end   = MPI_Wtime();
        t_total_write += t_end - t_start;

        /* Time vector reads for the benchmark */
        t_start = t_end;
        if ( H5FDread_vector(lf, H5P_DEFAULT, count, &(types_0[0]), &(addrs_0[0]), 
                       &(sizes_0[0]), &(read_bufs_0[0])) < 0 )
            TEST_ERROR;
        t_end   = MPI_Wtime();
        t_total_read += t_end - t_start;

    }
    if ( H5FDclose(lf) < 0 )
        TEST_ERROR;

    /* Close the fapl */
    if(H5Pclose(fapl_id) < 0)
        TEST_ERROR;
        
    if (mpi_rank == 0) {        /* The calculation must account for the number of loops (10) */
        bw = (total_size * 10) / t_total_write / (1024.0 * 1024.0);
        printf("Write Perf: %lf BW/[MBs] %ld Bytes AvgTime[sec] %lf\n", bw, total_size, t_total_write);
        bw = (total_size * 10) / t_total_read / (1024.0 * 1024.0);
        printf("Read  Perf: %lf BW/[MBs] %ld Bytes AvgTime[sec] %lf\n", bw, total_size, t_total_read);
        fflush(stdout);
    }
    return 0;

error:
    return -1;
}


/*
 *------------------------------------------------
 * This test can be run with the following
 * input arguments:
 * testName <directory_path> <data_size_per_rank> 
 *------------------------------------------------
 */


int
main( int argc, char **argv)
{
    int ret, nerrs = 0;
    int mpi_provides, require = MPI_THREAD_MULTIPLE;
    int64_t elements_per_rank = COUNT;
    hid_t fcpl, fapl;
    char data_pathname[FILENAME_BUF_SIZE];
    char data_filename[FILENAME_BUF_SIZE];

    if ( (MPI_Init_thread(&argc, &argv, require, &mpi_provides)) != MPI_SUCCESS) {
       HDfprintf(stderr, "FATAL: Unable to initialize MPI\n");
       HDexit(EXIT_FAILURE);
    }

    if ( (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank)) != MPI_SUCCESS) {
        HDfprintf(stderr, "FATAL: MPI_Comm_rank returned an error\n");
        HDexit(EXIT_FAILURE);
    }

    if ( (MPI_Comm_size(MPI_COMM_WORLD, &mpi_size)) != MPI_SUCCESS) {
        HDfprintf(stderr, "FATAL: MPI_Comm_size returned an error\n");
        HDexit(EXIT_FAILURE);
    }

    if (argc > 2) {
        DIR *directory = opendir(argv[2]);
        if (directory != NULL) {
            int len = strlen(argv[2]);
            closedir(directory);
            strcpy(data_pathname, argv[2]);
            if (data_pathname[len-1] == '/') {
                data_pathname[len-1] = '\0';
            }
            sprintf(data_filename,"%s/%s.h5", data_pathname, FILENAME[0]);
        }
        else {
            nerrs++;
            if (mpi_rank == 0) {
                puts("Error: input arg is not a directory name");
            }
            goto done;
        }
    } else {
        sprintf(data_filename,"%s/%s.h5", getcwd(data_pathname, FILENAME_BUF_SIZE), FILENAME[0]);
    }

    if (argc > 1) {
        int64_t value_check = atoll(argv[1]);
        if (value_check <= 0) {
            if (mpi_rank == 0) {
                puts("WARNING: Input arg is not a valid data size. Using defaults!");
            }
        }
        else {
            elements_per_rank = value_check;
        }
    }
    if (mpi_rank == 0) {
        printf("Test filename = %s\n", data_filename);
    }

    H5open();
    
    if (test_subfiling() < 0) {
        nerrs++;
        goto done;
    }

    ret = vector_io_bench(data_filename, elements_per_rank);

    MPI_Barrier(MPI_COMM_WORLD);

    if ( mpi_rank == 0 ) {           /* only process 0 reports */
        const char *header = "File open and read/write tests";

        HDfprintf(stdout, "===================================\n");
        if ( nerrs > 0 ) {
            HDfprintf(stdout, "***%s detected %d failures***\n", header, nerrs);
        }
        else {
            HDfprintf(stdout, "%s finished with no failures\n", header);
        }
        HDfprintf(stdout, "===================================\n");
    }

    /* close HDF5 library */
    if (H5close() != SUCCEED) {
        HDfprintf(stdout, "H5close() failed. (Ignoring)\n");
    }

done:
    /* MPI_Finalize must be called AFTER H5close which may use MPI calls */
    MPI_Finalize();

    /* cannot just return (nerrs) because exit code is limited to 1byte */
    return((nerrs > 0) ? EXIT_FAILURE : EXIT_SUCCESS );
}
