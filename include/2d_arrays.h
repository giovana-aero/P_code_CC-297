#ifndef _lib_2d_arrays_
#define _lib_2d_arrays_

void copy_2d_array(int m,int n,double A[m][n],double B[m][n]);
void fill_2d_array(int m,int n,double A[m][n]);
void ones_2d_array(int m,int n,double A[m][n]);
void print_2d_array(int m,int n,double A[m][n]);
void print_2d_array_to_file(int m,int n,double A[m][n],char *filename,
                            int invert_y);
void zeros_2d_array(int m,int n,double A[m][n]);

#endif