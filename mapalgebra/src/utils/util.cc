#ifndef UTIL_CC
#define UTIL_CC
/** util.cc: common utility functions needed by other functions
 * Author: Yan Y. Liu <yanliu@illinois.edu>
 * Date: 06/18/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <cstring>
#include <iostream>
#include <string>
#include <fstream>
#include <cblas.h>

#include "util.h"

using std::string;
using std::ifstream;

#define ERR(e) {fprintf(stderr, "Error: %s\n", e); return-1;}

Jobstat jobstat; // stat 


/*	* read a kernel file, the file format should be in csv format,
 *	* each line indicates a row of the kernel matrix.
 *	* 
 *	* Parameters
 *	*	filename: char pointer of the filename
 *	*	rows: number of rows in the matrix
 *	*	cols: number of columns in the matirx
 *	*	ret: the result matrix, should be set to NULL when calling
 *	*
 *	* Return value
 *	* 	-1: error happens
 *	*	0: success
 */
int kernelReader(char* filename, int rows, int cols, float*& ret){

  if(ret!=NULL)
    ERR("Error: in reader, array pointer is not NULL in input.");
  if(1!=rows%2||1!=cols%2)
    ERR("Error: in reader, row/column number is not an odd number.");
  if(rows!=cols)
    ERR("Error: in reader, kernel is not a square.");

  ifstream file (filename);
  if(!file.is_open())
    ERR("ERROR: in reader, invalid input file name.");
  ret=(float*)malloc(sizeof(float)*rows*cols);
  string line;
  int ind=0;
  for(int i=0;i<rows;i++) {
    if(file.eof()){
      free(ret);
      ret=NULL;
      ERR("ERROR: in reader, invalid kernel row size.");
    }

    getline(file, line);
    //line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
    size_t curPos=0, nextPos=0;
    for(int j=0;j<cols-1;j++){
      nextPos=line.find(",", curPos);
      if(nextPos==string::npos){
        free(ret);
        ret=NULL;
        ERR("ERROR: in reader, invalid kernel column size.");
      }
      float tmp=atof(line.substr(curPos, nextPos-curPos).c_str());
      ret[ind++]=tmp;
      curPos=nextPos+1;
    }
    float tmp=atof(line.substr(curPos).c_str());
    ret[ind++]=tmp;
  }

  return 0;
}

int kernelWriter(float* kernel, int rows, int cols) {

  // let user check the parsed input kernel matrix
  fprintf(stdout, "please check your kernel matrix\n");
  for(int i = 0; i < rows; i++) {
      for(int j =0; j< cols; j++) {
	  fprintf(stdout, "%0.2f ", kernel[i*rows + j]);
      }
      fprintf(stdout, "\n");
  }
  return 0;
}

// distance calcuation: considers elevation. must have the same unit for all
// coordinates
float get_dist(float x1, float y1, float z1, float x2, float y2, float z2)
{
	return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}

// get current system time
double get_timemark()
{
	struct timeval tsec;
	struct timezone tzone;
	gettimeofday(&tsec, &tzone);
	return (double)(tsec.tv_sec + tsec.tv_usec/1000000.0);
}

void print_jobstat_cuda() {
    fprintf(stdout, "====================\n");
    fprintf(stdout, "||    Job Stat    ||\n");
    fprintf(stdout, "====================\n");
    fprintf(stdout, "\nNO DEVICE INIT re-times first CUDA call\n\n");
    fprintf(stdout, "Reading time: %.5lf seconds\n", jobstat.Tread);
    fprintf(stdout, "Data distribution time: %.5lf seconds\n", jobstat.Tcommdata);
    fprintf(stdout, "Data distribution time (NO DEVICE INIT): %.5lf seconds\n",
            jobstat.Tcommdata_nd);
    fprintf(stdout, "Running time: %.5lf seconds\n", jobstat.Tcompute);
    fprintf(stdout, "Result collection time: %.5lf seconds\n",
            jobstat.Tcommresult);
    fprintf(stdout, "Writing time: %.5lf seconds\n", jobstat.Twrite);
    fprintf(stdout, "TOTAL TIME: %.5lf seconds\n", jobstat.Ttotal);
    fprintf(stdout, "TOTAL TIME (NO DEVICE INIT): %.5lf seconds\n",
            jobstat.Ttotal_nd);
}

void print_jobstat() 
{
	fprintf(stdout, "====================\n");
	fprintf(stdout, "||    Job Stat    ||\n");
	fprintf(stdout, "====================\n");
	fprintf(stdout, "Reading time: %.5lf seconds\n", jobstat.Tread);
	fprintf(stdout, "Data distribution time: %.5lf seconds\n", jobstat.Tcommdata);
	fprintf(stdout, "Running time: %.5lf seconds\n", jobstat.Tcompute);
	fprintf(stdout, "Result collection time: %.5lf seconds\n", jobstat.Tcommresult);
	fprintf(stdout, "Writing time: %.5lf seconds\n", jobstat.Twrite);
	fprintf(stdout, "TOTAL TIME: %.5lf seconds\n", jobstat.Ttotal);
}

void print_jobstat_csv()
{
	fprintf(stdout, "====================\n");
	fprintf(stdout, "||    Job Stat    ||\n");
	fprintf(stdout, "====================\n");
	fprintf(stdout, "Reading time, Data distribution time, Running time, Result collection time, Writing time, TOTAL TIMR\n");
	fprintf(stdout, "%.5lf seconds, %.5lf seconds, %.5lf seconds, %.5lf seconds, %.5lf seconds, %.5lf seconds\n", jobstat.Tread, jobstat.Tcommdata, jobstat.Tcompute, jobstat.Tcommresult, jobstat.Twrite, jobstat.Ttotal);  
}
// find the rectangle size closest to the square root of np
int get_best_dim(int np, int *rowSize, int *colSize)
{
	int i = 0; // loop counter
	int p = (int) sqrt(np);
	int q = np / p;
	while (p > 1 &&(p * q != np || p > q)) {
		p --;
		q = np / p;
		i++;
	}
	*rowSize = p;
	*colSize = q;
	return i;
}

int get_chunk(int rank, int size, int *&y_size_list, int *&y_off_list, int *&buffer_size_list, int *&displacement_list, int *&scatter_buffer_size_list, int *&scatter_displacement_list, int *&mode_list, int raster_x_size, int raster_y_size, int kernel_y_size, int &scatter_raster_size)
{
  // calculate overlapping rows
  int overlap_row_num = (kernel_y_size - 1)/2;
  // number of rows that each process is assigned with
  int proc_y_size = raster_y_size/size;
  // array  of row offset for each process
  y_off_list = (int *)malloc(sizeof(int) * size);
  // array of row number for each process
  y_size_list = (int *)malloc(sizeof(int) * size);
  // array of result buffer size of each process, only matters for root process
  buffer_size_list = (int *)malloc(sizeof(int) * size);    
  // displacement between each result buffer, only matters for root process
  displacement_list = (int *)malloc(sizeof(int) * size);
  // array of chunk buffer size of each process, only matters for root process
  scatter_buffer_size_list = (int *)malloc(sizeof(int) * size);    
  // displacement between each chunk  buffer, only matters for root process
  scatter_displacement_list = (int *)malloc(sizeof(int) * size);
  mode_list = (int *)malloc(sizeof(int) * size);
  for(int i =0; i< size; i++)
    {
      if(i!=size-1){
	displacement_list[i] = proc_y_size * i * raster_x_size;
	buffer_size_list[i] = proc_y_size * raster_x_size;
      }
      else{
	displacement_list[i] = proc_y_size * i * raster_x_size;
	buffer_size_list[i] = (raster_y_size-proc_y_size*(size-1))*raster_x_size;
      }
      
      if(size==1){//only one job running
	mode_list[0]=4;
	y_off_list[0]=0;
	y_size_list[0]=raster_y_size;
      }
      else if(i == 0)
	{
	  mode_list[i] = 1;
	  y_off_list[i] = 0;
	  y_size_list[i] = proc_y_size + overlap_row_num;
	}
      else if(i == size - 1)
	{
	  mode_list[i] = 3;
	  y_off_list[i] = i*proc_y_size - overlap_row_num; 
	  y_size_list[i] = raster_y_size - proc_y_size*(size-1) + overlap_row_num;
	}
      else
	{
	  mode_list[i] = 2;
	  y_off_list[i] = i*proc_y_size - overlap_row_num;
	  y_size_list[i] = proc_y_size + 2*overlap_row_num;
	}
      if(i == 0)
	{
	  scatter_displacement_list[i] = 0;
	}
      else
	{ 
	  scatter_displacement_list[i] = scatter_displacement_list[i-1] + scatter_buffer_size_list[i-1];	  
	}
      scatter_buffer_size_list[i] = y_size_list[i] * raster_x_size;
#ifdef DEBUG
      std::cout<<i<<" y_off "<<y_off_list[i]<<" y_size "<<y_size_list[i]<<" x size "<<raster_x_size<<std::endl;
#endif
    }
  scatter_raster_size = scatter_displacement_list[size-1] + scatter_buffer_size_list[size-1];
  return 1;
}

// calc the data block to be read given the rank of the process
// rank, np: rank and number of procs
// x, y: size of raster 
// offsetx, offsety: output. 
// sizex, sizey: output. useful for procs positioned at the end of each dim
int get_block(int rank, int np, int x, int y, int *offsetx, int *offsety, int * sizex, int *sizey)
{
	int dimx, dimy, blockx, blocky;
	int cellsx, cellsy, rcellsx, rcellsy; // num cells on 2 dims and remainders
	get_best_dim(np, &dimy, &dimx); // number of blocks on each dimension
	blockx = rank % dimx; // block id on x dim
	blocky = rank / dimx; // block id on y dim
#ifdef DEBUG
	fprintf(stderr, "np %d gets %d x %d blocks; rank %d blocky=%d blockx=%d\n", np, dimy, dimx, rank, blocky, blockx);
#endif
	cellsx = x / dimx;
	cellsy = y / dimy;
	rcellsx = x % dimx;
	rcellsy = y % dimy;
	*offsetx = blockx * cellsx;
	*offsety = blocky * cellsy;
	*sizex = cellsx;
	*sizey = cellsy;
	if (blockx == dimx - 1) { // end block on x dim
		*sizex = *sizex + rcellsx;
	}
	if (blocky == dimy - 1) { // end block on y dim
		*sizey = *sizey + rcellsy;
	}
	return 1;
}

/* Function to do submatrix addition, Let A,B be the full matrices, a,b be their
 * sub-matrices, the function calculates a=a*beta*b.  
 * 	lda: the column number of A
 * 	ldb: the column number of B
 * 	A: the float array in row vector form of A
 * 	B: the float array in row vector form of B
 */ 
void fsm_add(int m, int n, float* A, int lda, float beta, float* B, int ldb){
	for(int y=0;y<m;y++){
		for(int x=0;x<n;x++){
			A[y*lda+x]+=beta*B[y*ldb+x];
		}
	}
};

/*	* Process an image(matrix) with a kernel(matrix). The parameters pRet, ret_row, ret_col
	* are modified to show the output matrix.
	* 
	* Parameters
	*	pImage: 	array of a image in row vector form.
	*	img_row:	number of rows for the input image.
	*	img_col:	number of columns for the input image.
	*	pKernel:	array of a kernel in row vector form.
	*	ker_row:	number of rows for the kernel matrix.
	*	ker_col:	number of columns for the kernel matrix.
	*	pRet: 		array of output matrix in row vector form. Must be NULL in input. Set to NULL
					if the function is in error.
	*	ret_row:	number of rows of the output matrix.
	*	ret_col: 	number of columns of the output matrix.
	*	mode: 		different ways to produce the result matrix.
	*					1, the image matrix is the uppest one. 
	*					2, the image matrix is the middle one. 
	*					3, the image matrix is the bottom one. 
	*					4, the image matrix is not seperated(used when only one job is running)
	*
	* Return value
	*	-1 for fail, 0 for success.
	* 	If success, pRet would be an allocated result array in row vector form, ret_row/ret_col
	*	would be the size of the result matrix. User should free pRet manually. 
*/
int process(float* pImage, int img_row, int img_col, float* pKernel, int ker_row, int ker_col,
						float*& pRet, int& ret_row, int& ret_col, int mode){
	//check for inputs
	if(ker_col==0||ker_row==0)
		ERR("invalid kernel size 0!");
	if(img_col==0||img_row==0)
		ERR("invalid image size 0!");
	if(img_col<ker_col||img_row<ker_row)
		ERR("invalid image size, smaller than kernel!");
	if(1!=ker_col%2||1!=ker_row%2)
		ERR("kernel height or width is not odd!");
	if(pRet!=NULL)
		ERR("input matrix pointer is not NULL!");
	
	//check for different modes
	if(mode==1){		//for uppest submatrix of an image
		ret_col=img_col;
		ret_row=img_row-ker_row/2;
		pRet=(float*)calloc(sizeof(float), ret_col*ret_row);
		for(int x=-ker_col/2;x<=ker_row/2;x++){
			for(int y=-ker_row/2;y<=ker_row/2;y++){		//loop for every value in the kernel
				int img_x, img_y, ret_x, ret_y, h, w;
				if(x<=0 && y<=0){
					img_x=img_y=0;
					ret_x=-x;
					ret_y=-y;
					h=ret_row+y;
					w=ret_col+x;
				}else if(x<=0 && y >0){
					img_x=0;
					img_y=y;
					ret_x=-x;
					ret_y=0;
					h=ret_row;
					w=ret_col+x;
				}else if(x>0 && y<=0){
					img_x=x;
					img_y=0;
					ret_x=0;
					ret_y=-y;
					h=ret_row+y;
					w=ret_col-x;
				}else{
					img_x=x;
					img_y=y;
					ret_x=ret_y=0;
					h=ret_row;
					w=ret_col-x;
				}
				//for each value in kernel matrix, do a submatrix element-wise addition
				int i=x+ker_col/2;
				int j=y+ker_row/2;
				fsm_add(h, w, &pRet[ret_y*ret_col+ret_x], ret_col, pKernel[j*ker_col+i], &pImage[img_y*img_col+img_x], img_col);
			}	
		}
	}else if(mode==2){		//for middle submatrix of an image
		ret_col=img_col;
		ret_row=img_row-ker_row/2*2;
		pRet=(float*)calloc(sizeof(float), ret_col*ret_row);
		for(int x=-ker_col/2;x<=ker_row/2;x++){
			for(int y=-ker_row/2;y<=ker_row/2;y++){	
				int img_x, img_y, ret_x, ret_y, h, w;
				if(x<=0 && y<=0){
					img_x=0;
					img_y=ker_row/2+y;
					ret_x=-x;
					ret_y=0;
					h=ret_row;
					w=ret_col+x;
				}else if(x<=0 && y >0){
					img_x=0;
					img_y=ker_row/2+y;
					ret_x=-x;
					ret_y=0;
					h=ret_row;
					w=ret_col+x;
				}else if(x>0 && y<=0){
					img_x=x;
					img_y=ker_row/2+y;
					ret_x=0;
					ret_y=0;
					h=ret_row;
					w=ret_col-x;
				}else{
					img_x=x;
					img_y=ker_row/2+y;
					ret_x=ret_y=0;
					h=ret_row;
					w=ret_col-x;
				}
				int i=x+ker_col/2;
				int j=y+ker_row/2;
				fsm_add(h, w, &pRet[ret_y*ret_col+ret_x], ret_col, pKernel[j*ker_col+i], &pImage[img_y*img_col+img_x], img_col);
			}
		}
	}else if(mode==3){		//for bottom submatrix of an image
		ret_col=img_col;
		ret_row=img_row-ker_row/2;
		pRet=(float*)calloc(sizeof(float), ret_col*ret_row);
		for(int x=-ker_col/2;x<=ker_row/2;x++){
			for(int y=-ker_row/2;y<=ker_row/2;y++){		
				int img_x, img_y, ret_x, ret_y, h, w;
				if(x<=0 && y<=0){
					img_x=0;
					img_y=ker_row/2+y;
					ret_x=-x;
					ret_y=0;
					h=ret_row;
					w=ret_col+x;
				}else if(x<=0 && y >0){
					img_x=0;
					img_y=ker_row/2+y;
					ret_x=-x;
					ret_y=0;
					h=ret_row-y;
					w=ret_col+x;
				}else if(x>0 && y<=0){
					img_x=x;
					img_y=ker_row/2+y;
					ret_x=0;
					ret_y=0;
					h=ret_row;
					w=ret_col-x;
				}else{
					img_x=x;
					img_y=ker_row/2+y;
					ret_x=ret_y=0;
					h=ret_row-y;
					w=ret_col-x;
				}
				int i=x+ker_col/2;
				int j=y+ker_row/2;
				fsm_add(h, w, &pRet[ret_y*ret_col+ret_x], ret_col, pKernel[j*ker_col+i], &pImage[img_y*img_col+img_x], img_col);
			}
		}
	}else if(mode==4){		//for a full image, used when only one job is running
		ret_col=img_col;
		ret_row=img_row;
		pRet=(float*)calloc(sizeof(float), ret_col*ret_row);
		for(int x=-ker_col/2;x<=ker_row/2;x++){
			for(int y=-ker_row/2;y<=ker_row/2;y++){		
				int img_x, img_y, ret_x, ret_y, h, w;
				if(x<=0){
					img_x=0;
					ret_x=-x;
					w=ret_col+x;
				}else{
					img_x=x;
					ret_x=0;
					w=ret_col-x;
				}
				if(y<=0){
					img_y=0;
					ret_y=-y;
					h=ret_row+y;
				}else{
					img_y=y;
					ret_y=0;
					h=ret_row-y;
				}
				int i=x+ker_col/2;
				int j=y+ker_row/2;
				fsm_add(h, w, &pRet[ret_y*ret_col+ret_x], ret_col, pKernel[j*ker_col+i], &pImage[img_y*img_col+img_x], img_col);
			}
		}
	}
	else{
		ERR("invalid mode, should be 1(uppest), 2(middle), or 3(downnest).");
	}
	return 0;
}
#endif
