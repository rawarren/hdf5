#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "hdf5.h"

#include "mpi.h"

int
main(int argc, char **argv)
{

	int i, mpi_size, mpi_rank;
	int loop_count = 20;
	int mpi_provides, require = MPI_THREAD_MULTIPLE;
	hid_t subfile_id = 1;
	double t_start, t_end, t_open, t_close;
	const char *h5_filename = "unused.h5";
	FILE *h5file;

	MPI_Init_thread(&argc, &argv, require, &mpi_provides);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	if (argc > 1) {
		int check_value = atoi(argv[1]);
		if (check_value > 0) {
			loop_count = check_value;
		}
	}

    H5open();
	
	h5file = fopen(h5_filename, "w+");
	for(i=0; i < loop_count; i++) {
		t_start = MPI_Wtime();
		sf_open_subfiles(subfile_id, h5_filename, NULL, O_CREAT|O_TRUNC|O_RDWR);
		t_end = MPI_Wtime();
		t_open = t_end - t_start;
		t_start = t_end;
		sf_close_subfiles(subfile_id);
		t_end = MPI_Wtime();
		t_close = t_end - t_start;
		if (mpi_rank == 0) {
			printf("loop_count(%d) ... time = %lf = %lf + %lf seconds\n", i, t_open + t_close, t_open, t_close);
			fflush(stdout);
		}
	}

	if (h5file) {
		fclose(h5file);
		if (mpi_rank == 0)
			unlink(h5_filename);
	}

	MPI_Barrier(MPI_COMM_WORLD);

    H5close();

	MPI_Finalize();
	return 0;
}
