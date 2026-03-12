/*
- gigiaero, 10/03/2026, 1346 hours
*/

void zeros_2d(int m,int n,double A[m][n]){
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++)
      A[i][j] = 0.;
  }
}