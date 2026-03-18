#include<stdio.h>
#include"../../include/1d_arrays.h"

void operator1(double *x){
  *x = (*x)*2;
}

void operator2(double *x){
  *x = (*x)*5;
}

void operator3(double *x){
  *x = (*x)*10;
}

double operator0(double x){
  return x + 1.;
}

void loop(void (*f)(double*),int n,double *v){
  for(int i=0;i<n;i++)
    (*f)(&v[i]);
}

void execute(int op,int n,double *v){

  void (*f)(double *);

  switch(op){
    case 1:
      f = operator1;
      break;

    case 2:
      f = operator2;
      break;

    case 3:
      f = operator3;
      break;
  }

  for(int i=0;i<n;i++)
    (*f)(&v[i]);
}

void execute2(int op,int n,double *v){

  void (*f)(double *);

  switch(op){
    case 1:
      f = operator1;
      break;

    case 2:
      f = operator2;
      break;

    case 3:
      f = operator3;
      break;
  }

  loop(f,n,v);
}

int main(){

  int n = 5;
  int op = 1;
  double v[n];
  
  fill_1d_array(n,v);
  print_1d_array(n,v);
  execute(op,n,v);
  print_1d_array(n,v);
  putchar('\n');
  fill_1d_array(n,v);
  print_1d_array(n,v);
  execute2(op,n,v);
  print_1d_array(n,v);
  

  // double (*f)(double);
  // double x = 3.;
  // f = operator0;
  // x = f(3);
  // printf("%f\n",x);


  return 0;
}