#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"../include/2d_arrays.h"
#include"../include/ecm_lib.h"
#include"../include/eom_lib.h"

#define div_ref 1e100
#define pi 3.1415926535897932384626433

void evaluate_delta_form_ecm(sim_prmtrs *config,msh_prmtrs *msh,
                             control_prmtrs *c_prmtrs,int init_only){
  int n = msh->IMAX+2*msh->JMAX;
  int m = msh->JMAX;
  double (*x)[n] = calloc(m,sizeof *x);
  double (*y)[n] = calloc(m,sizeof *y);

  initialize_mesh_ecm(m,n,x,y,msh);

  if(config->save_i_c){
    char *filename = malloc(sizeof(char)*200);

    sprintf(filename,"%s%s",config->casename,"_x_initial.dat");
    print_2d_array_to_file(m,n,x,filename,0);
    sprintf(filename,"%s%s",config->casename,"_y_initial.dat");
    print_2d_array_to_file(m,n,y,filename,0);

    free(filename);
  }

  if(!init_only){
    // puts("ADI, elliptical C mesh");
    solve_adi_2d_rectangular_ecm(msh->JMAX,msh->IMAX,x,y,config,c_prmtrs);
  }

  free(x);
  free(y);
}

void half_ellipse(double *x,double *y,double *prmtrs,int n,int invert_th){
  double *th = malloc(sizeof(double)*n);

  if(invert_th)
    linspace(th,3.*pi/2.,pi/2.,n);
  else
    linspace(th,pi/2.,3.*pi/2.,n);

  for(int i=0;i<n;i++){
    x[i] = prmtrs[2] + prmtrs[0]*cos(th[i]);
    y[i] = prmtrs[3] + prmtrs[1]*sin(th[i]);
  }
    
  free(th);
}

void initialize_mesh_ecm(int m,int n,double x[m][n],double y[m][n],
                       msh_prmtrs *msh){
  int chord_n = (msh->IMAX+1)/2;
  double *x_axis = malloc(sizeof(double)*chord_n);
  double *tmp_x = malloc(sizeof(double)*msh->IMAX);
  double *tmp_y = malloc(sizeof(double)*msh->IMAX);
  double *tmp_x_wake = malloc(sizeof(double)*(m+1));
  cosspace(x_axis,0.,msh->c,chord_n,0);

  msh->end_prmtrs[2] = msh->c;
  
  switch(msh->af_type){
    case 1:
      init_af_bi_air(tmp_x,tmp_y,x_axis,chord_n,msh);
      break;

    case 2:
      init_af_naca4(tmp_x,tmp_y,x_axis,chord_n,msh);
      break;

    case 3:
      init_af_cst(tmp_x,tmp_y,x_axis,chord_n,msh);
      break;

    default:
      puts("initialize_mesh: invalid af_type");
      exit(13);
  }

  for(int i=0,k=m;i<msh->IMAX;i++,k++){
    x[0][k] = tmp_x[i];
    y[0][k] = tmp_y[i];
    x[m-1][k] = tmp_x[i];
    y[m-1][k] = tmp_y[i];
  }

  half_ellipse(tmp_x,tmp_y,msh->end_prmtrs,msh->IMAX,1);

  for(int i=0,k=m;i<msh->IMAX;i++,k++){
    x[m-1][k] = tmp_x[i];
    y[m-1][k] = tmp_y[i];
  }

  cosspace(tmp_x_wake,msh->end_prmtrs[2],msh->c+msh->end_prmtrs[0],m+1,1);
  
  for(int i=0,k1=m,k2=msh->IMAX+m;i<m;i++,k1--,k2++){
    x[0][i] = tmp_x_wake[k1];
    x[0][k2] = tmp_x_wake[i+1];
    y[0][i] = 0.;
    y[0][k2] = 0.;
    x[m-1][i] = tmp_x_wake[k1];
    x[m-1][k2] = tmp_x_wake[i+1];
    y[m-1][i] = -msh->end_prmtrs[0];
    y[m-1][k2] = msh->end_prmtrs[0];
  }

  init_type3(m,n,x,y,msh);

  free(tmp_x);
  free(tmp_y);
  free(tmp_x_wake);
  free(x_axis);
}