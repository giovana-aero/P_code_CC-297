/* 
Fills a 2d array with increasing values
- gigiaero, 10/03/2026, 1344 hours
*/

void fill_2d_array(int m,int n,double A[m][n]){
  int k = 0;
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
      A[i][j] = (double) k;
      k++;
    }
  }
}
