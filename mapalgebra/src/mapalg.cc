/** mapalg.cc: map algebra
 * Author: Mingze Gao <mgao16@illinois.edu>
 * Date: 09/09/2015
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <errno.h>

#include <mpi.h>
#include <gdal.h>
#include <cpl_conv.h>
#include <cpl_string.h>

#include "util.h"
#include "data.h"

// options
bool perf=false, csv=false, debug=false;

// each process is assigned with one block
typedef struct {
    int maxx, maxy;
    int bsizex, bsizey;
    int boffsetx, boffsety;
    double nodata;
    float* data;
} Block;

// block topo information
typedef struct {
    int *offsetx, *offsety;
    int *sizex, *sizey;
    int maxx, maxy;
} Position;

// only root process handle the Raster
typedef struct {
    double georef[6]; // georef data struct for a raster
    char prj[2048];   // store projection wkt
    double nodata;
    int x, y; // size of raster on x and y dim
    float *data;
    Position p;
} Raster;


void show_help();
bool is_float(char* arg);
void get_position(Position &p, int x, int y, int np);
void populate_raster(Raster &raster, char* raster_file, int rank, int np);
void transfer_data(Block &block,Raster &raster);
void gather_data(Block &block, Raster &raster);
void write_raster(Raster &raster, char* output_file_name, int np);
void clean_raster(Raster &raster, int rank);
void clean_block(Block &b);
void show_perf(double t0, double t1, double t2,double t3, 
        double t4, double t5);
void distribute_function_call(int argc, char** argv, 
        float (*f)(float, float), int rank, int np);
void distance(char* raster1, float xp, float yp, float zp,  
        char* output_file,int rank, int np);
void setnull_proc(char* raster1, float value, char *con_op, float con_val,
        char* output_file, int rank, int np);
void setnull2_proc(char* raster1, char * raster2, char *con_op, float con_val,
        char* output_file, int rank, int np);
void con_proc(char* raster1, char *con_op, float con_val, float value,
        char* output_file, int rank, int np);

// write performation information to stdou
void show_perf(double t0, double t1, double t2,double t3, double t4, double t5) {
    jobstat.Tread = t1 - t0;
    jobstat.Tcommdata = t2 - t1;
    jobstat.Tcompute = t3 - t2;
    jobstat.Tcommresult = t4 - t3;
    jobstat.Twrite = t5 - t4;
    jobstat.Ttotal = t5 - t0;
    if(csv)
        print_jobstat_csv();
    else
        print_jobstat();
}

// free allocated memory of a Raster
void clean_raster(Raster &raster, int rank) {
    // only root raster has real data
    if(rank == 0) {
        free(raster.data);
        free(raster.p.sizex);
        free(raster.p.sizey);
        free(raster.p.offsetx);
        free(raster.p.offsety);
    }
}

void clean_block(Block &b) {
    free(b.data);
}

// create output raster file and write result
void write_raster(Raster &raster, char* output_file_name, int np) {
    GDALDatasetH rout;
    float *block;
    rout = raster_create(output_file_name, raster.x, raster.y,
                         raster.georef, raster.prj, raster.nodata);
    for(int i=0; i<np; i++) {
        block = raster.data + i * (raster.p.maxx *raster.p.maxy);
        raster_write(rout, block, raster.p.offsetx[i], raster.p.offsety[i],
                     raster.p.sizex[i], raster.p.sizey[i]);
    }
    raster_close(rout);
}


// root proc gathers results
void gather_data(Block &block, Raster &raster) {
    MPI_Gather(block.data, block.maxx*block.maxy, MPI_FLOAT, raster.data,
               block.maxx*block.maxy, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}

// root proc transfer blocks to other proc
void transfer_data(Block &block,Raster &raster) {
    block.maxx = raster.p.maxx;
    block.maxy = raster.p.maxy;
    block.nodata = raster.nodata;
    MPI_Bcast(&(block.maxx), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(block.maxy), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(block.nodata), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(raster.p.offsetx, 1, MPI_INT, &(block.boffsetx),
                1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(raster.p.offsety, 1, MPI_INT, &(block.boffsety),
                1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(raster.p.sizex, 1, MPI_INT, &(block.bsizex),
                1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(raster.p.sizey, 1, MPI_INT, &(block.bsizey),
                1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    block.data = (float *) malloc(sizeof(float)*block.maxx*block.maxy);
    MPI_Scatter(raster.data, block.maxx * block.maxy, MPI_FLOAT, block.data,
                block.maxx*block.maxy, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}
// root process determine block sizes
void get_position(Position &p, int x, int y, int np) {
    p.maxx=0, p.maxy=0;
    p.offsetx = (int *) malloc(sizeof(int) * np);
    memset(p.offsetx, 0, sizeof(int) * np);
    p.offsety = (int *) malloc(sizeof(int) * np);
    memset(p.offsety, 0, sizeof(int) * np);
    p.sizex = (int *) malloc(sizeof(int) * np);
    memset(p.sizex, 0, sizeof(int) * np);
    p.sizey = (int *) malloc(sizeof(int) * np);
    memset(p.sizey, 0, sizeof(int) * np);
    for (int i=0; i<np; i++) {
        get_block(i, np, x, y, &(p.offsetx[i]), &(p.offsety[i]),
                  &(p.sizex[i]), &(p.sizey[i]));
    }
    // find max sizex and sizey
    for (int i=0; i<np; i++) {
        if (p.maxx < p.sizex[i]) p.maxx = p.sizex[i];
        if (p.maxy < p.sizey[i]) p.maxy = p.sizey[i];
    }
}

// root process reads dataset
void populate_raster(Raster &raster, char* raster_file, int rank, int np) {
    GDALDatasetH rin;
    float* block;
    rin = raster_open(raster_file, raster.georef, raster.prj,
                      &(raster.nodata), &(raster.x), &(raster.y));
    if(debug && rank == 0) {
        raster_info(rin, raster.georef, raster.prj, raster.nodata,
        raster.x, raster.y);
    }
    if(rank == 0) {
        get_position(raster.p, raster.x, raster.y, np);
        raster.data = (float *)malloc(sizeof(float)
                                      * raster.p.maxx * raster.p.maxy * np);
        memset(raster.data, 0, sizeof(float)
                                      * raster.p.maxx * raster.p.maxy * np);
        for(int i=0; i<np; i++) {
            block = raster.data + i * (raster.p.maxx * raster.p.maxy);
            raster_read(rin, block, raster.p.offsetx[i], raster.p.offsety[i],
                        raster.p.sizex[i], raster.p.sizey[i]);
        }
    }
    raster_close(rin);
}

// check if argument can be convert to float number
bool is_float(char* arg) {
    float d = strtof(arg, NULL);  // try to convert str to float
    if (d == (0.0F)) {
        return false;  // failed in conversion
    } else {
        if (errno == ERANGE) {
            perror("Invalid float number"); // float is out of range
            MPI_Finalize();
            exit(1);
        }
        return true;
    }
}

// raster raster operations
void raster_raster(char* raster1, char* raster2, char* output_file, float (*f)(float, float), int rank, int np) {

    Raster r1, r2;
    Block b1, b2;
    double t0, t1, t2, t3, t4, t5;
    t0 = get_timemark();
    // step 1 : read input raster
    populate_raster(r1, raster1, rank, np);
    populate_raster(r2, raster2, rank, np);
    t1 = get_timemark();
    // step 2 : transfer data blocks to procs
    transfer_data(b1, r1);
    transfer_data(b2, r2);
    if(debug)
        fprintf(stderr, "rank %d: max[x,y]=%d,%d nodata=%.5lf size=%d,%d\n", rank,
                b1.maxx, b1.maxy, b1.nodata, b1.bsizex, b1.bsizey);
    t2 = get_timemark();
    // map algebra operation
    int index;
    for(int i=0; i<b1.bsizey; i++) {
        for(int j=0; j<b1.bsizex; j++) {
            index = i * b1.bsizex + j;
            if(fabs((double)(b1.data[index]) - b1.nodata) < 0.001)
                b1.data[index] = b1.nodata;
            else
                b1.data[index] = (*f)(b1.data[index], b2.data[index]);
        }
    }
    t3 = get_timemark();

    // step4 : transfer results for writing
    gather_data(b1, r1);
    t4 = get_timemark();

    // step5 : write output
    if (rank == 0)
        write_raster(r1, output_file, np);
    t5 = get_timemark();

    if (perf && rank == 0)
        show_perf(t0, t1, t2, t3, t4, t5);

    // step6 : clean up
    clean_raster(r1, rank);
    clean_raster(r2, rank);
    clean_block(b1);
    clean_block(b2);
}

void setnull_proc(char* raster1, float value, char *con_op, float con_val,
                  char* output_file, int rank, int np) {
    Raster r1;
    Block b1;
    double t0, t1, t2, t3, t4, t5;
    t0 = get_timemark();
    // step1 : read input raster
    if (rank == 0) {
        populate_raster(r1, raster1, rank, np);
    }
    t1 = get_timemark();

    // step2 : transfer data blocks to procs
    transfer_data(b1, r1);
    if (debug)
        fprintf(stderr, "rank %d: max[x,y]=%d,%d nodata=%.5lf size=%d,%d\n", rank,
                b1.maxx, b1.maxy, b1.nodata, b1.bsizex, b1.bsizey);
    t2 = get_timemark();

    // step3 : map algebra operation
    int index;
    for (int i=0; i<b1.bsizey; i++) {
        for (int j=0; j<b1.bsizex; j++) {
            index = i * b1.bsizex + j;
            if (fabs((float)(b1.data[index]) - b1.nodata) < 0.001)
                continue;
            else {
                if (strcmp(con_op,"gt")==0) {
                    if (b1.data[index] > con_val)
                        b1.data[index] = b1.nodata;
                    else
                        b1.data[index] = value;
                } else if (strcmp(con_op,"ge")==0) {
                    if (b1.data[index] >= con_val)
                        b1.data[index] = b1.nodata;
                    else
                        b1.data[index] = value;
                } else if (strcmp(con_op,"lt")==0) {
                    if (b1.data[index] < con_val)
                        b1.data[index] = b1.nodata;
                    else
                        b1.data[index] = value;
                } else if (strcmp(con_op,"le")==0) {
                    if (b1.data[index] <= con_val)
                        b1.data[index] = b1.nodata;
                    else
                        b1.data[index] = value;
                }
            }
        }
    }

    t3 = get_timemark();
    // step4 : transfer results for writing
    gather_data(b1, r1);
    t4 = get_timemark();
    // step5 : write output
    if (rank == 0) {
        write_raster(r1, output_file, np);
    }
    t5 = get_timemark();
    if (perf && rank==0)
        show_perf(t0, t1, t2, t3, t4, t5);
    clean_raster(r1, rank);
    clean_block(b1);
}

void setnull2_proc(char* raster1, char * raster2, char *con_op, float con_val,
                   char* output_file, int rank, int np) {
    Raster r1, r2;
    Block b1, b2;
    double t0, t1, t2, t3, t4, t5;

    t0 = get_timemark();

    // step 1 : read input raster
    populate_raster(r1, raster1, rank, np);
    populate_raster(r2, raster2, rank, np);
    t1 = get_timemark();

    // step 2 : transfer data blocks to procs
    transfer_data(b1, r1);
    transfer_data(b2, r2);
    if(debug)
        fprintf(stderr, "rank %d: max[x,y]=%d,%d nodata=%.5lf size=%d,%d\n",
                rank, b1.maxx, b1.maxy, b1.nodata, b1.bsizex, b1.bsizey);
    t2 = get_timemark();

    // map algebra operation
    int index;
    for(int i=0; i<b1.bsizey; i++) {
        for(int j=0; j<b1.bsizex; j++) {
            index = i * b1.bsizex + j;
            if (fabs((float)(b1.data[index]) - b1.nodata) < 0.001)
                b1.data[index] = b1.nodata;
            else {
                if (strcmp(con_op,"gt")==0) {
                    if (b1.data[index] > con_val)
                        b1.data[index] = b1.nodata;
                    else
                        b1.data[index] = b2.data[index];
                } else if (strcmp(con_op,"ge")==0) {
                    if (b1.data[index] >= con_val)
                        b1.data[index] = b1.nodata;
                    else
                        b1.data[index] = b2.data[index];
                } else if (strcmp(con_op,"lt")==0) {
                    if (b1.data[index] < con_val)
                        b1.data[index] = b1.nodata;
                    else
                        b1.data[index] = b2.data[index];
                } else if (strcmp(con_op,"le")==0) {
                    if (b1.data[index] <= con_val)
                        b1.data[index] = b1.nodata;
                    else
                        b1.data[index] = b2.data[index];
                }
            }
        }
    }
    t3 = get_timemark();

    // step4 : transfer results for writing
    gather_data(b1, r1);
    t4 = get_timemark();

    // step5 : write output
    if (rank == 0)
        write_raster(r1, output_file, np);
    t5 = get_timemark();
    if (perf && rank == 0)
        show_perf(t0, t1, t2, t3, t4, t5);

    // step6 : clean up
    clean_raster(r1, rank);
    clean_raster(r2, rank);
    clean_block(b1);
    clean_block(b2);
}

void con_proc(char* raster1, char *con_op, float con_val, float value,
              char* output_file, int rank, int np) {
    Raster r1;
    Block b1;
    double t0, t1, t2, t3, t4, t5;
    t0 = get_timemark();
    // step1 : read input raster
    if (rank == 0) {
        populate_raster(r1, raster1, rank, np);
    }
    t1 = get_timemark();
    // step2 : transfer data blocks to procs
    transfer_data(b1, r1);
    if(debug)
        fprintf(stderr, "rank %d: max[x,y]=%d,%d nodata=%.5lf size=%d,%d\n", rank,
                b1.maxx, b1.maxy, b1.nodata, b1.bsizex, b1.bsizey);
    t2 = get_timemark();
    // step3 : map algebra operation
    int index;
    for (int i=0; i<b1.bsizey; i++) {
        for (int j=0; j<b1.bsizex; j++) {
            index = i * b1.bsizex + j;
            if (fabs((float)(b1.data[index]) - b1.nodata) < 0.001)
                continue;
            else {
                if (strcmp(con_op,"gt")==0) {
                    if (b1.data[index] > con_val)
                        b1.data[index] = value;
                } else if (strcmp(con_op,"ge")==0) {
                    if (b1.data[index] >= con_val)
                        b1.data[index] = value;
                } else if (strcmp(con_op,"lt")==0) {
                    if (b1.data[index] < con_val)
                        b1.data[index] = value;
                } else if (strcmp(con_op,"le")==0) {
                    if (b1.data[index] <= con_val)
                        b1.data[index] = value;
                }
            }
        }
    }

    t3 = get_timemark();
    // step4 : transfer results for writing
    gather_data(b1, r1);
    t4 = get_timemark();
    // step5 : write output
    if (rank == 0) {
        write_raster(r1, output_file, np);
    }
    t5 = get_timemark();
    if (perf && rank==0)
        show_perf(t0, t1, t2, t3, t4, t5);
    clean_raster(r1, rank);
    clean_block(b1);
}

void distance(char* raster1, float xp, float yp, float zp,  char* output_file,
              int rank, int np) {
    Raster r1;
    Block b1;
    double t0, t1, t2, t3, t4, t5;
    t0 = get_timemark();
    // step1 : read input raster
    if (rank == 0) {
        populate_raster(r1, raster1, rank, np);
    }
    t1 = get_timemark();
    // step2 : transfer data blocks to procs
    transfer_data(b1, r1);
    if(debug)
        fprintf(stderr, "rank %d: max[x,y]=%d,%d nodata=%.5lf size=%d,%d\n", rank,
                b1.maxx, b1.maxy, b1.nodata, b1.bsizex, b1.bsizey);
    t2 = get_timemark();
    // step3 : map algebra operation
    int index;
    float xc, yc, zc;
    for(int i=0; i<b1.bsizey; i++) {
        for(int j=0; j<b1.bsizex; j++) {
            index = i * b1.bsizex + j;
            if(fabs((double)(b1.data[index]) - b1.nodata) < 0.001) continue;
            xc = r1.georef[0] + (b1.boffsetx + j) * r1.georef[1];
            yc = r1.georef[3] + (b1.boffsety + i) * r1.georef[5];
            zc = b1.data[index];
            b1.data[index] = get_dist(xc, yc, zc, xp, yp, zp);
        }
    }
    t3 = get_timemark();
    // step4 : transfer results for writing
    gather_data(b1, r1);
    t4 = get_timemark();
    // step5 : write output
    if (rank == 0) {
        write_raster(r1, output_file, np);
    }
    t5 = get_timemark();
    if (perf && rank==0)
        show_perf(t0, t1, t2, t3, t4, t5);
    clean_raster(r1, rank);
    clean_block(b1);
}

// raster number operations
void raster_number(char* raster1, float num, char* output_file, float (*f)(float, float),
                   bool regular_order, int rank, int np) {

    Raster r1;
    Block b1;
    double t0, t1, t2, t3, t4, t5;

    t0 = get_timemark();

    // step1 : read input raster
    populate_raster(r1, raster1, rank, np);
    t1 = get_timemark();

    // step2 : transfer data blocks to procs
    transfer_data(b1, r1);
    if(rank == 0 && debug)
        fprintf(stderr, "rank %d: max[x,y]=%d,%d nodata=%.5lf size=%d,%d\n", rank,
                b1.maxx, b1.maxy, b1.nodata, b1.bsizex, b1.bsizey);
    t2 = get_timemark();

    // step3 : map algebra operation
    int index;
    for(int i=0; i<b1.bsizey; i++) {
        for(int j=0; j<b1.bsizex; j++) {
            index = i * b1.bsizex + j;
            if(fabs((double)(b1.data[index]) - b1.nodata) < 0.001)
                b1.data[index] = b1.nodata;
            else {
                if (regular_order)
                    b1.data[index] = (*f)(b1.data[index], num);
                else
                    b1.data[index] = (*f)(num, b1.data[index]);
            }
        }
    }
    t3 = get_timemark();

    // step4 : transfer results for writing
    gather_data(b1, r1);
    t4 = get_timemark();

    // step5 : write output
    if (rank == 0) {
        write_raster(r1, output_file, np);
    }
    t5 = get_timemark();

    // wrap-up
    if (perf && rank==0)
        show_perf(t0, t1, t2, t3, t4, t5);
    clean_raster(r1, rank);
    clean_block(b1);
}

inline float plus(float x, float y) {
    return  x+y;
}

inline float minus(float x, float y) {
    return x-y;
}

inline float multiplies(float x, float y) {
    return x*y;
}

inline float divides(float x, float y) {
    return x/y;
}

inline float power(float base, float exponent) {
    return pow(base, exponent);
}

inline float rms_func(float x, float y) {
    return pow(0.5 * (x * x + y * y), 0.5);
}

void distribute_function_call(int argc, char** argv,
                              float (*f)(float, float), int rank, int np) {

    // operation with two number
    if(is_float(argv[2]) && is_float(argv[3])) {

        if (rank == 0) {
            fprintf(stdout, "%f\n", 
                    (*f)(strtof(argv[2],NULL),strtof(argv[3], NULL)));
        }

    }

    // operation of two raster
    else if((! is_float(argv[2])) && (! is_float(argv[3]))) {
        raster_raster(argv[2], argv[3], argv[4], f, rank, np);
    }

    else if (is_float(argv[2])) {
        raster_number(argv[3], strtof(argv[2], NULL), argv[4], f, false,
                      rank, np);
    }

    // operation of raster with number
    else {
        raster_number(argv[2], strtof(argv[3], NULL), argv[4], f, true,
                      rank, np);
    }
}

void sum(int argc, char** argv, int rank, int np) {
    if(argc != 5) {
        if (rank == 0) {
            fprintf(stderr, "Invalid number of argument\n");
            fprintf(stderr, "Usage: add <raster|number> "
                             "<raster|number> <out>\n");
        }
        MPI_Finalize();
        exit(1);
    }
    distribute_function_call(argc, argv, plus, rank, np);
}

void diff(int argc, char** argv, int rank, int np) {
    if(argc != 5) {
        if (rank==0) {
            fprintf(stderr, "Invalid number of argument\n");
            fprintf(stderr, "Usage: diff <raster|number> "
                            "<raster|number> <out>\n");
        }
        exit(1);
    }
    distribute_function_call(argc, argv, minus, rank, np);
}

void mult(int argc, char** argv, int rank, int np) {
    if(argc != 5) {
        if (rank == 0) {
            fprintf(stderr, "Invalid number of argument\n");
            fprintf(stderr, "Usage: mult <raster|number> "
                            "<raster|number> <out>\n");
        }
        MPI_Finalize();
        exit(1);
    }
    distribute_function_call(argc, argv, multiplies, rank, np);
}

void div(int argc, char** argv, int rank, int np) {
    if(argc != 5) {
        if (rank==0) {
            fprintf(stderr, "Invalid number of argument\n");
            fprintf(stderr, "Usage: div <raster|number> "
                            "<raster|number> <out>\n");
        }
        MPI_Finalize();
        exit(1);
    }
    distribute_function_call(argc, argv, divides, rank, np);
}

void exp(int argc, char** argv, int rank, int np) {
    if(argc != 5) {
        if (rank==0) {
            fprintf(stderr, "Invalid number of argument\n");
            fprintf(stderr, "Usage: exp <raster|number> "
                            "<raster|number> <out>\n");
        }
        MPI_Finalize();
        exit(1);
    }
    distribute_function_call(argc, argv, power, rank, np);
}

void sqrt(int argc, char**argv, int rank, int np) {
    if(argc != 4) {
        if (rank==0) {
            fprintf(stderr, "Invalid number of argument\n");
            fprintf(stderr, "Usage: [options] sqrt <raster|number> <out>\n");
        }
        MPI_Finalize();
        exit(1);
    }
    char point5[] = "0.5";
    char** updated_argv = argv;
    updated_argv[4] = argv[3];
    updated_argv[3] = point5;

    distribute_function_call(5, updated_argv, power, rank, np);
}

void rms(int argc, char** argv, int rank, int np) {
    if (argc < 5) {
        if (rank==0)  {
            fprintf(stderr, "Invalid number of argument\n");
            fprintf(stderr, "Usage: [options] rms <raster1> "
                            "<raster2> <out>\n");
        }
        MPI_Finalize();
        exit(1);
    }

    distribute_function_call(argc, argv, rms_func, rank, np);
}

void dist(int argc, char**argv, int rank, int np) {
    if(argc != 7) {
        if (rank == 0) {
            fprintf(stderr, "Invalid number of argument\n");
            fprintf(stderr, "Usage: dist <raster> <x> <y> <z> <out>\n");
        }
        MPI_Finalize();
        exit(1);
    }

    else if (! (is_float(argv[3]) && is_float(argv[4]) && is_float(argv[5]))) {
        if (rank == 0) {
            fprintf(stderr, "Invalid argument\n");
            fprintf(stderr, "Usage: dist <raster> <x> <y> <z> <out>\n");
        }
        MPI_Finalize();
        exit(1);
    } else {
        distance(argv[2], strtof(argv[3],NULL), strtof(argv[4],NULL),
                 strtof(argv[5],NULL), argv[6], rank, np);
    }
}

void setnull(int argc, char**argv, int rank, int np) {
    char set_null_usage [] = "Usage: setnull <raster> <false val> "
                              "<condition operator> <condition val> <out>\n";
    if(argc != 7) {
        if (rank == 0) {
            fprintf(stderr, "Invalid number of argument\n");
            fprintf(stderr, set_null_usage);
        }
        MPI_Finalize();
        exit(1);
    }

    else if (! (is_float(argv[3]) && is_float(argv[5]))) {
        if (rank == 0) {
            fprintf(stderr, "Invalid argument\n");
            fprintf(stderr, set_null_usage);
        }
        MPI_Finalize();
        exit(1);

    } else if (strcmp(argv[4],"gt")!=0 && strcmp(argv[4],"ge")!=0
               && strcmp(argv[4],"lt")!=0 && strcmp(argv[4],"le")!=0) {
        printf("The condition operator is not valid!\n");
        MPI_Finalize();
        exit(0);

    } else {
        float con_val = atof(argv[5]);
        char *con_op = argv[4];
        char * infn = argv[2];
        float value = atof(argv[3]);
        char * outfn = argv[6];
        setnull_proc(infn, value, con_op, con_val, outfn, rank, np);
    }


}
void setnull2(int argc, char**argv, int rank, int np) {
    char  set_null_usage [] = "Usage: setnull2 <raster1> <raster2> "
                              "<condition operator> <condition val> <out>\n";
    if(argc != 7) {
        if (rank == 0) {
            fprintf(stderr, "Invalid number of argument\n");
            fprintf(stderr, set_null_usage);
        }
        MPI_Finalize();
        exit(1);

    } else if (! is_float(argv[5])) {
        if (rank == 0) {
            fprintf(stderr, "Invalid argument %s, must be a float number\n",
                    argv[5]);
            fprintf(stderr, set_null_usage);
        }
        MPI_Finalize();
        exit(1);

    } else if (strcmp(argv[4],"gt")!=0 && strcmp(argv[4],"ge")!=0
               && strcmp(argv[4],"lt")!=0 && strcmp(argv[4],"le")!=0) {
        printf("The condition operator is not valid!\n");
        MPI_Finalize();
        exit(0);

    } else {
        float con_val = atof(argv[5]);
        char *con_op = argv[4];
        char * infn = argv[2];
        char * infn2 = argv[3];
        char * outfn = argv[6];
        setnull2_proc(infn, infn2, con_op, con_val, outfn, rank, np);
    }
}

void con(int argc, char**argv, int rank, int np) {
    char  usage [] = "Usage: con <raster> <condition operator> "
                     "<condition val> <value> <out>\n";

    if(argc != 7) {
        if (rank == 0) {
            fprintf(stderr, "Invalid number of argument\n");
            fprintf(stderr, usage);
        }
        MPI_Finalize();
        exit(1);

    } else if (! (is_float(argv[4]) && is_float(argv[5]))) {
        if (rank == 0) {
            fprintf(stderr, "Invalid argument\n");
            fprintf(stderr, usage);
        }
        MPI_Finalize();
        exit(1);

    } else if (strcmp(argv[3],"gt")!=0 && strcmp(argv[3],"ge")!=0
               && strcmp(argv[3],"lt")!=0 && strcmp(argv[3],"le")!=0) {
        printf("The condition operator is not valid!\n");
        MPI_Finalize();
        exit(0);

    } else {
        float con_val = atof(argv[4]);
        char *con_op = argv[3];
        char * infn = argv[2];
        float value = atof(argv[5]);
        char * outfn = argv[6];
        con_proc(infn, con_op, con_val, value, outfn, rank, np);
    }

}

/****************************************************************
 *
 * Parse the command line options specific to Map Algebra
 *
 * parse_options should be called after MPI_Init which will
 * MPI options after processing. Similarly parse_options should
 * remove any options found and resulting argv should have
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
            opt_count += 1;
        }

        else if (!strcmp(argv[i], "--perf")) {
            perf = true;
            opt_count += 1;
        }

        else if (!strcmp(argv[i], "--csv")) {
            csv = true;
            opt_count += 1;
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


/****************************************************************
 *
 * Method Switch checks the remaining command line argumment and
 * passes control to appropriate method to handle the reading of
 * raster files, specific processing, and writing of results.
 */
void method_switch(int argc, char** argv, int rank, int np) {
    std::string method;
    if(argc < 2) {
        if (rank == 0) {
            fprintf(stderr, "Invalid number of argument\n");
        }
        MPI_Finalize();
        exit(1);
    }
    method = argv[1];

    if (method == "add") {
        sum(argc, argv, rank, np);
    }

    else if (method == "diff") {
        diff(argc, argv, rank, np);
    }

    else if (method == "mult") {
        mult(argc, argv, rank, np);
    }

    else if (method == "div") {
        div(argc, argv, rank, np);
    }

    else if (method == "exp") {
        exp(argc, argv, rank, np);
    } else if (method == "rms") {
        rms(argc, argv, rank, np);
    }

    else if (method == "dist") {
        dist(argc, argv, rank, np);
    } else if (method == "sqrt") {
        sqrt(argc, argv, rank, np);
    }

    else if (method == "setnull") {
        setnull(argc, argv, rank, np);
    }

    else if (method == "setnull2") {
        setnull2(argc, argv, rank, np);
    }

    else if (method == "con") {
        con(argc, argv, rank, np);
    }

    else {
        if (rank == 0) {
            fprintf(stderr, "Unknown method %s\n", argv[1]);
        }
        MPI_Finalize();
        exit(1);
    }
}

// print help information
void show_help() {
    const char *help_info =
        "Usage: mapalg [options] <method> <args> [<args>...] <output>\n"
        "\nOptions:\n"
        "  --help, -h   : this help section\n"
        "  --debug, -d  : debug information written to stderr\n"
        "  --perf       : performance information written to stdout\n"
        "  --csv        : performance written to stdout in CSV format\n"
        "\nMethods:\n"
        "  add | diff | mult | div | exp "
        "<raster|number> <raster|number> <output>\n"
        "               : outputs the basic cell-wise algebraic operations\n"
        "                 between input1 and input2 where inputs may be \n"
        "                 either a raster file or a number\n\n"
        "  sqrt <raster>\n"
        "               : outputs the square root of input raster\n\n"
        "  rms <raster1> <raster2>\n"
        "               : outputs the root mean square of input rasters\n\n"
        "  dist <raster> <x> <y> <z>\n"
        "               : outputs the distance between each cell and point\n"
        "                 where x, y, and z are given in the units of the\n"
        "                 input raster\n\n"
        "  setnull <raster> <false val> <con op> <con val>\n"
        "               : output valid cell values to null and the false\n"
        "                 value is from a constant\n\n"
        "  setnull2 <raster> <raster> <con op> <con val>\n"
        "               : set valid cell values to null and the false value "
        "                 is from a raster\n\n"
        "  con <raster> <con op> <con val> <value>\n"
        "               : set valid cell values to a constant\n\n";

    fprintf(stdout, help_info);
}


// assumption: float value; same dimensions
int main(int argc, char** argv) {
    static int rank, np;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    parse_options(&argc, argv, rank, np);
    method_switch(argc, argv, rank, np);

    MPI_Finalize();
}
