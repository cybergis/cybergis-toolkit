#ifndef UTIL_H
#define UTIL_H
/** util.h: common utility functions needed by other functions
 * Author: Yan Y. Liu <yanliu@illinois.edu>
 * Date: 06/18/2014
 */

typedef struct {
	double Tread;
	double Tcommdata;
	double Tcommdata_nd;
	double Tcompute;
	double Tcommresult;
	double Twrite;
	double Ttotal;
	double Ttotal_nd;
} Jobstat;
extern Jobstat jobstat;

int parseArg(int argc, char** argv, char *inputFileName, char *outputFileName, char *kernelFileName,
             int &kernelXSize, int &kernelYSize, int &band);
int process(float* pImage, int img_row, int img_col, float* pKernel, int ker_row, int ker_col,
						float*& pRet, int& ret_row, int& ret_col, int mode);
int kernelReader(char* filename, int rows, int cols, float*& ret);
float get_dist(float x1, float y1, float z1, float x2, float y2, float z2);
double get_timemark();
void print_jobstat_cuda() ;
void print_jobstat() ;
void print_jobstat_csv();
int get_best_dim(int np, int *rowSize, int *colSize);
int get_block(int rank, int np, int x, int y, int *offsetx, int *offsety, int * sizex, int *sizey);
int get_chunk(int rank, int size, int *&y_size_list, int *&y_off_list, int *&buffer_size_list, 
	      int *&displacement_list, int *&scatter_buffer_size_list, int *&scatter_displacement_list, 
	      int *&mode_list, int raster_x_size, int raster_y_size, int kernel_y_size, int &scatter_raster_size);

#endif
