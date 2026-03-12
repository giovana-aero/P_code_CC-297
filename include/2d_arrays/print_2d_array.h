/*
- gigiaero, 10/03/2026, 1346 hours
*/

void print_2d_array(int m,int n,double A[m][n]){
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++)
      printf("%f ",A[i][j]);
    putchar('\n');
  }
}
