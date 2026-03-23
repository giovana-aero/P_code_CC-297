#ifndef _lib_1d_arrays_
#define _lib_1d_arrays_

void copy_1d_array_range(int i_1,int i_f,double *v1,double *v2);
void fill_1d_array(int n,double *v);
void ones_1d_array(int n,double *v);
void print_1d_array(int n,double *v);
void print_1d_array_int(int n,int *v);
void print_1d_array_to_file(int n,double *v,char *filename);
void zeros_1d_array(int n,double *v);

#endif