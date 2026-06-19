#include<stdio.h>
#include<stdlib.h>
#include"../include/2d_arrays.h"

/*
Copy array A to B
- gigiaero, 19/03/2026, 2200 hours
*/
void copy_2d_array(int m,int n,double A[m][n],double B[m][n]){
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      B[j][i] = A[j][i];
    
  }
}

/*
-gigiaero, 31/03/2026, 0925 hours
*/
void eye(int m,double A[m][m]){
  for(int i=0;i<m;i++)
    A[i][i] = 1.;
}

/* 
Fills a 2d array with increasing values
- gigiaero, 10/03/2026, 1344 hours
*/
void fill_2d_array(int m,int n,double A[m][n]){
  int k = 0;
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++){
      A[j][i] = (double) k;
      k++;
    }
  }
}

/*
- gigiaero, 10/06/2026, 1620 hours
*/
void flip_2d_array(int m,int n,double A[m][n]){
  double (*tmp)[n] = calloc(m,sizeof *tmp);

  copy_2d_array(m,n,A,tmp);

  for(int j1=0,j2=m-1;j1<m;j1++,j2--){
    for(int i1=0,i2=n-1;i1<n;i1++,i2--)
      A[j1][i1] = tmp[j2][i2];
  }

  free(tmp);
}

/*
- gigiaero, 20/05/2026, 0842 hours
*/
double max_2d_array(int m,int n,double A[m][n]){
  double k = 0.;
  
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++){
      if(A[j][i] > k)
        k = A[j][i];
    }
  }

  return k;
}

/*
- gigiaero, 10/03/2026, 1345 hours
*/
void ones_2d_array(int m,int n,double A[m][n]){
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      A[j][i] = 1.;
  }
}

/*
- gigiaero, 10/03/2026, 1346 hours
*/
void print_2d_array(int m,int n,double A[m][n]){
  // if(invert_y){
  //   for(int j=m-1;j>=0;j--){
  //     for(int i=0;i<n;i++)
  //       printf("%f ",A[j][i]);
  //     putchar('\n');
  //   }
  // }
  // else{
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      printf("%f ",A[j][i]);
    putchar('\n');
  }
  // }
}

/*
- gigiaero, 20/03/2026, 0058 hours
*/
void print_2d_array_to_file(int m,int n,double A[m][n],char *filename,
                            int invert_y){
  FILE *output;

  output = fopen(filename,"w");

  if(invert_y){
    for(int j=m-1;j>=0;j--){
      for(int i=0;i<n;i++)
        fprintf(output,"%.6f ",A[j][i]);
      
      fprintf(output,"\n");
    }  
  }
  else{
    for(int j=0;j<m;j++){
      for(int i=0;i<n;i++)
        fprintf(output,"%.6f ",A[j][i]);
      
      fprintf(output,"\n");
    }
  }

  fclose(output);
}

/*
- gigiaero, 28/05/2026, 2030 hours
*/
void read_2d_array_from_file(int m,int n,double A[m][n],char* filename){
  FILE *input;

  int i = 0,j = 0;

  input = fopen(filename,"r");
  puts(filename);

  while(1){
    if(fscanf(input,"%lf",&A[j][i]) == EOF)
      break;

    i++;
    
    if(i == n){
      i = 0;
      j++;
    }
  }

  fclose(input);
}

/*
- gigiaero, 10/03/2026, 1346 hours
*/
void zeros_2d_array(int m,int n,double A[m][n]){
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      A[j][i] = 0.;
  }
}
