#include<stdlib.h>
#include"../../include/2d_arrays.h"

int main(){
  char filename[] = "../../P02/results/eom_x_initial.dat";
  int m = 15;
  int n = 93;

  double (*A)[n] = calloc(m,sizeof *A);


  read_2d_array_from_file(m,n,A,filename);
  print_2d_array(m,n,A);

  free(A);

  return 0;
}