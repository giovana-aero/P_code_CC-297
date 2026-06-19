#ifndef _lib_2d_arrays_
#define _lib_2d_arrays_

void copy_2d_array(int m,int n,double A[m][n],double B[m][n]);
void eye(int m,double A[m][m]);
void fill_2d_array(int m,int n,double A[m][n]);
void flip_2d_array(int m,int n,double A[m][n]);
double max_2d_array(int m,int n,double A[m][n]);
void ones_2d_array(int m,int n,double A[m][n]);
void print_2d_array(int m,int n,double A[m][n]);
void print_2d_array_to_file(int m,int n,double A[m][n],char *filename,
                            int invert_y);
void read_2d_array_from_file(int m,int n,double A[m][n],char* filename);
void zeros_2d_array(int m,int n,double A[m][n]);

#endif