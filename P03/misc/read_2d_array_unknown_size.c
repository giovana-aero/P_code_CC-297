#include<stdio.h>
// #include<stdlib.h>

int main(){
  char filename[] = "../../P02/results/eom_x_initial.dat";
  FILE *input;
  double x;

  input = fopen(filename,"r");

  int i;
  int j = 0;
  int n = 0;


  while(1){
    if(fscanf(input,"%lf",&x) == EOF){
      puts("break!");
      j++;
      break;
    }
    n++;
    printf("%f\n",x);
    // getchar();
  }

  fclose(input);

  input = fopen(filename,"r");

  while(!feof(input)){
    if(fgetc(input) == '\n')
      j++; 
  }

  j--; // disregard last line break

  i = n/j;

  printf("%d | %d | %d\n",i,j,n);

  return 0;
}