#include<stdio.h>

void print_2d_array(int m,int n,double A[m][n]){
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++)
      printf("%f ",A[i][j]);
    putchar('\n');
  }
}

void ones_2d(int m,int n,double A[m][n]){
  int k = 0;
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++)
      A[i][j] = 1.;
  }
}

void zeros_2d(int m,int n,double A[m][n]){
  int k = 0;
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++)
      A[i][j] = 0.;
  }
}

void fill_2d_array(int m,int n,double A[m][n]){
  int k = 0;
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
      A[i][j] = (double) k;
      k++;
    }
  }
}

int main(){

  int m = 5; // lines
  int n = 7; // columns

  double A[m][n];

  // int k = 0;
  // for(int i=0;i<m;i++){
  //   for(int j=0;j<n;j++){
  //     A[i][j] = (double) k;
  //     k++;
  //   }
  // }

  

  // for(int i=0;i<m;i++){
  //   for(int j=0;j<n;j++)
  //     printf("%f ",A[i][j]);
  //   putchar('\n');
  // }

  ones_2d(m,n,A);
  print_2d_array(m,n,A);
  putchar('\n');
  
  zeros_2d(m,n,A);
  print_2d_array(m,n,A);
  putchar('\n');

  fill_2d_array(m,n,A);
  print_2d_array(m,n,A);
  putchar('\n');

  return 0;
  
}