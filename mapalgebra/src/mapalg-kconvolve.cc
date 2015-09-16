#include <iostream>
#include <cstring>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include "ogr_spatialref.h"
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "util.h"
#include "data.h"

#include "mpi.h"


using std::string;

void parse_options(int *argc, char** argv, int rank, int np);
void show_help();

bool perf = false, csv=false, debug=false;
/****************************************************************
 *
 * Parse the command line options specific to Map Algebra
 *
 * parse_options should be called after MPI_Init which will
 * remove MPI options after processing. Similarly parse_options
 * should remove any options found and resulting argv should have
 * the method as argv[1].
 *
 * Note: modifying the argv is a questionable practice, but MPI_Init
 * does it so I guess we can too.
 */
void parse_options(int *argc, char** argv, int rank, int np) {
    std::string option;
    int opt_count = 0;

    for(int i = 1; i<*argc; i++) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            if (rank == 0)
                show_help();
            MPI_Finalize();
            exit(0);
        }

        else if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--debug")) {
            debug = true;
            opt_count++;

        }

        else if (!strcmp(argv[i], "--perf")) {
            perf = true;
            opt_count++;

        }

        else if (!strcmp(argv[i], "--csv")) {
            csv = true;
            opt_count++;

        }

        else if (*argv[i] == '-') {
            if (rank == 0) {
                fprintf(stderr, "Error: unknown option '%s'\n\n", argv[i]);
                show_help();
            }
            MPI_Finalize();
            exit(1);
        }
    }

    // update argc and argv for return
    *argc -= opt_count;
    memmove(argv+1, argv+1+opt_count, *argc * sizeof (char *));
}


void show_help() {
    const char *help_info =
        "Usage: mapalg-kernel [options] <raster> <kernel> <x size> <y size> <out>\n"
        "\nOptions:\n"
        "    -h | --help   : this help section\n"
        "    -d | --debug  : debug information to stderr\n"
        "    --perf        : performance information written to stdout\n"
        "    --csv         : performance written to stdout in CSV format\n"
        "\nArguments:\n"
        "    raster        : input raster filename\n"
        "    kernel        : kernel filename\n"
        "    xsize         : kernel size in x direction\n"
        "    ysize         : kernel size in y direction\n"
        "    out           : output raster filename\n"
        ;

    printf(help_info);
}

// assumption: float value; kernel row size and column size are odd numbers
int main(int argc, char *argv[]) {
    // get mpi rank and size
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    parse_options(&argc, argv, rank, size);

    double t0, t1, t2, t3, t4, t5; // timing info

    t0 = get_timemark();

    // step 1 : read input raster
    // user defined arguments
    char * input_file_name;
    char * output_file_name;
    char * kernel_file_name;
    int kernel_x_size, kernel_y_size; // dimension of input kernel matrix

    // data only root process knows and handles

    // off set list and data chunk size list
    int *y_off_list = NULL, *y_size_list = NULL;
    float * result_buffer = NULL;      // final output buffer
    int result_buffer_size;            // final output buffer size
    int *buffer_size_list = NULL, *displacement_list = NULL, *mode_list = NULL;
    int *scatter_displacement_list = NULL, *scatter_buffer_size_list = NULL;
    double georef[6]; // georef data structure for a raster
    char prj[2048];  // store projection wkt
    double nodata; // value for missing data
    GDALDatasetH inDataset; // input data set
    float *raster=NULL;
    int scatter_raster_size;

    // data every proc knows
    float *kernel = NULL;
    int x_size = 0; // number of columns for each process
    int y_size = 0; // number of rows for each process
    int mode = 0;
    int raster_x_size = 0, raster_y_size = 0; // size of input raster dat
    float * proc_chunk = NULL; //data chunk for each proc
    float * chunk = NULL; // tmp pointer to assess each proc's raster chunk
    float* proc_result=NULL; // each procs's output buffer
    int proc_result_x_size;
    int proc_result_y_size;
    int bandInd = 1;

    // command line arguments
    if (argc != 6) {
        if (rank == 0) {
            fprintf(stderr, "Error: incorrect number of arguments.\n\n");
            show_help();
        }
        MPI_Finalize();
        exit(1);
    }

    input_file_name = argv[1];
    kernel_file_name = argv[2];
    kernel_x_size = atoi(argv[3]);
    kernel_y_size = atoi(argv[4]);
    output_file_name = argv[5];


    // read the kernel file
    kernel = (float *)malloc(sizeof(float)*kernel_x_size*kernel_y_size);
    // kernel = NULL;   WHY IS THIS HERE?
    kernelReader(kernel_file_name, kernel_x_size, kernel_y_size, kernel);

    // step 1 : root proc read process data
    if(rank == 0) {
        inDataset = raster_open(input_file_name, georef, prj, &nodata,
                                &raster_x_size, &raster_y_size, bandInd);

        // root process print raster info
        raster_info(inDataset, georef, prj, nodata, raster_x_size, raster_y_size);
        result_buffer_size = raster_x_size * raster_y_size;
        result_buffer = (float *)malloc(sizeof(float) * (result_buffer_size));

        // now root process determine data chunk size for each process
        get_chunk(rank, size, y_size_list, y_off_list, buffer_size_list,
                  displacement_list, scatter_buffer_size_list,
                  scatter_displacement_list, mode_list, raster_x_size,
                  raster_y_size, kernel_y_size, scatter_raster_size);
        x_size = raster_x_size; // get x size of each procs

        // now root procs read the input raster
        raster = (float *)malloc(sizeof(float) * scatter_raster_size);
        for(int i = 0; i< size; i++) {
            chunk = raster + scatter_displacement_list[i];
            raster_read(inDataset, chunk, 0, y_off_list[i], x_size,
                        y_size_list[i]);
        }
    }

    // root process broadcast row size list and mode list
    MPI_Scatter(y_size_list, 1, MPI_INT, &y_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(mode_list, 1, MPI_INT, &mode, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&x_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // each proc malloc enough memory for input data chunk
    proc_chunk = (float *)malloc(sizeof(float) * x_size * y_size);

    t1 = get_timemark();

    // step 2 : transfer data chunks to procs
    MPI_Scatterv(raster, scatter_buffer_size_list, scatter_displacement_list,
                 MPI_FLOAT, proc_chunk, x_size * y_size, MPI_FLOAT, 0,
                 MPI_COMM_WORLD);

    // each proc process assigned data chunk from input raster
    t2 = get_timemark();

    // step 3: kernel convolution computation
    // each process process assigned chunk data and store the result to buffer
    process(proc_chunk, y_size, x_size, kernel, kernel_y_size, kernel_x_size,
            proc_result, proc_result_y_size, proc_result_x_size, mode);
    int proc_result_buf_length = proc_result_x_size * proc_result_y_size;

    t3 = get_timemark();

    // step 4: transter results for writing
    MPI_Gatherv(proc_result, proc_result_buf_length, MPI_FLOAT, result_buffer,
                buffer_size_list, displacement_list,  MPI_FLOAT, 0,
                MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    t4 = get_timemark();

    if(debug) {
        // test
        if(rank == 0) {
            for(int i = 0; i < 10; i++) {
                for(int j=0; j<10; j++) {
                    std::cout<<" "<<result_buffer[i*raster_x_size+j];
                }
                std::cout<<"\n";
            }
        }
    }
    // step 5 : root process write accumulated result to output file
    if(rank == 0) {
        GDALDatasetH outDataset;
        outDataset = raster_create(output_file_name, raster_x_size,
                                   raster_y_size, georef, prj, nodata);
        raster_write(outDataset, result_buffer, 0, 0,
                     raster_x_size, raster_y_size);
        raster_close(outDataset);
    }
    t5 = get_timemark();

    // compute timing
    if(rank == 0) {
        jobstat.Tread = t1 - t0;
        jobstat.Tcommdata = t2 - t1;
        jobstat.Tcompute = t3 - t2;
        jobstat.Tcommresult = t4 - t3;
        jobstat.Twrite = t5 - t4;
        jobstat.Ttotal = t5 - t0;
        if(perf) {
            if (csv)
                print_jobstat_csv();
            else
                print_jobstat();
        }
    }

    // step 6: free all resource, clean up
    if(rank == 0) {
        free(y_off_list);
        free(y_size_list);
        free(buffer_size_list);
        free(mode_list);
        free(raster);
        free(result_buffer);
        free(displacement_list);
        free(scatter_displacement_list);
        free(scatter_buffer_size_list);
    }

    free(kernel);
    free(proc_chunk);
    free(proc_result);
    MPI_Finalize();
}
