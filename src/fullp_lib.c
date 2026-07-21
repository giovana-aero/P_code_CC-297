#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"../include/2d_arrays.h"
#include"../include/eom_lib.h"
#include"../include/fullp_lib.h"
#include"../include/num_methods.h"
#include"../include/strings.h"

#define div_ref 1e100
#define gamma 1.4
#define pi 3.14159265358979323846261003

void check_j(int m,int j){
  if(j >= m || j < 0){
    puts("check_j");
    // exit(100);
  }
}

/*
- gigiaero, 04/06/2026, 1621 hours
*/
void beta_initial_config(sim_prmtrs *config,fp_beta_prmtrs *fpb_prmtrs,
                         double beta){
  if(!config->alpha_seq)
    config->M = 1;
  
  fpb_prmtrs->L2 = malloc(sizeof(double)*(config->M+1));
  fpb_prmtrs->Linf = malloc(sizeof(double)*(config->M+1));

  fpb_prmtrs->beta_high = beta + 1.;
  fpb_prmtrs->beta_low = beta - 1.;
}

/*
- gigiaero, 04/06/2026, 1700 hours
*/
void beta_switch(double *beta,double beta_sub,double beta_super,int supersonic){
  if(supersonic)
    *beta = beta_super;
  else
    *beta = beta_sub;

  // if(supersonic)
  //   puts("supersonic");
  // else
  //   puts("subsonic");
}

/*
supersonic beta, that is
- gigiaero, 04/06/2026, 1714 hours
*/
void beta_update(double *beta,double L2_now,double Linf_now,int M,
                 fp_beta_prmtrs *fpb_prmtrs,int iter){
  double R;

  if(iter == 1){
    for(int i=0;i<M+1;i++){
      fpb_prmtrs->L2[i] = L2_now;
      fpb_prmtrs->Linf[i] = Linf_now;
    }
  }
  else{
    for(int i=M-1;i>=0;i--){
      fpb_prmtrs->L2[i+1] = fpb_prmtrs->L2[i];
      fpb_prmtrs->Linf[i+1] = fpb_prmtrs->Linf[i];
    }

    fpb_prmtrs->L2[0] = L2_now;
    fpb_prmtrs->Linf[0] = Linf_now;
  }

  R = fpb_prmtrs->L2[0]/fpb_prmtrs->L2[M] + 
      fpb_prmtrs->Linf[0]/fpb_prmtrs->Linf[M];
  
  if(R < 2.)
    (*beta) *= .98;

  else if(R > 2.1)
    (*beta) *= 1.1;

  if((*beta) > fpb_prmtrs->beta_high)
    *beta = fpb_prmtrs->beta_high;
  
  else if((*beta) < fpb_prmtrs->beta_low)
    *beta = fpb_prmtrs->beta_low;
}

/*
- gigiaero, 03/06/2026, 2105 hours
NOTE: adapted to af2 only for now

refactored
- gigiaero, 22/06/2026, 2256 hours
*/
double calc_Ai(int m,int n,double Ai[m][n],double rho_til[m][n],
               double A1_ih[m][n],double J_ih[m][n]){
  // for(int j=0+af2_flip;j<m-1+af2_flip;j++){
  for(int j=1;j<m;j++){
    for(int i=0;i<n-1;i++)
      Ai[j][i] = rho_til[j][i]*A1_ih[j][i]/J_ih[j][i];
    
    Ai[j][n-1] = Ai[j][0];
  }
}

/*
- gigiaero, 18/06/2026, 2109 hours
NOTE: adapted to af2 only for now

refactored
- gigiaero, 22/06/2026, 2301 hours
*/
double calc_Aj(int m,int n,double Aj[m][n],double rho_bar[m][n],
               double A3_jh[m][n],double J_jh[m][n]){
  for(int j=0;j<m;j++){
    for(int i=0;i<n-1;i++)
      Aj[j][i] = rho_bar[j][i]*A3_jh[j][i]/J_jh[j][i];
    
    Aj[j][n-1] = Aj[j][0];
  }
}

/*
- gigiaero, 12/07/2026, 2311 hours
*/
void calc_circulation(double *Gamma,int n,double *phi,double Rg){
  double Gamma_til,Gamma_bar,Gamma_estimate;

  Gamma_til = phi[1] - phi[n-2];

  Gamma_bar = Rg*Gamma_til + (1. - Rg)*Gamma[1];
  
  Gamma_estimate = 2*Gamma[1] - Gamma[2];

  Gamma[2] = Gamma[1];

  Gamma[1] = Gamma[0];

  Gamma[0] = (Gamma_bar + Gamma_estimate)*.5;
}

/*
- gigiaero, 31/05/2026, 1143h hours

refactored
- gigiaero, 21/06/2026, 1337 hours

refactored (again)
- gigiaero, 21/06/2026, 2202 hours

yeah ┐(￣ー￣)┌
- gigiaero, 23/06/2026, 1028
*/
void calc_contraUV(int m,int n,double contraUV[m][n],double phi[m][n],
                   double A1[m][n],double A2[m][n],double A3[m][n],int op,
                   int axis){
  double dphi_dksi,dphi_deta;
  int j_end;

  if(op == 2)
    j_end = m-1;
  else
    j_end = m;

  for(int j=0;j<j_end;j++){
    for(int i=0;i<n-1;i++){
      dphi_dksi = get_dphi_dksi(m,n,phi,i,j,op);
      dphi_deta = get_dphi_deta(m,n,phi,A2[j][i],A3[j][i],dphi_dksi,i,j,op);

      switch(axis){
        case 1: // contraU
          contraUV[j][i] = dphi_dksi*A1[j][i] + dphi_deta*A2[j][i];
          break;

        case 2: // contraV
          contraUV[j][i] = dphi_dksi*A2[j][i] + dphi_deta*A3[j][i];
          break;

        default:
          puts("calc_contraUV: invalid axis");
          exit(15);
      }
    }

    contraUV[j][n-1] = contraUV[j][0];
  }

  // Boundary conditions for the i,j+1/2 case
  if(op == 2){
    switch(axis){
      case 1:
        for(int i=0;i<n;i++)          
          contraUV[m-1][i] = contraUV[m-2][i];
      
      break;

      case 2:
        for(int i=0;i<n;i++)          
          contraUV[m-1][i] = -contraUV[m-2][i];

      break;

      default:
        puts("calc_contraUV: invalid axis");
        exit(15);
    }
  }
}

/*
- gigiaero, 23/06/2026, 0927 hours
*/
void calc_J_A_metrics(int m,int n,double J[m][n],double A1[m][n],
                      double A2[m][n],double A3[m][n],double x[m][n],
                      double y[m][n],int op){
  double dx_dksi,dx_deta;
  double dy_dksi,dy_deta;
  double ksix,ksiy;
  double etax,etay;

  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++){
      dx_dksi = get_dxy_dksi(m,n,x,i,j,op);
      dx_deta = get_dxy_deta(m,n,x,i,j,op);
      dy_dksi = get_dxy_dksi(m,n,y,i,j,op);
      dy_deta = get_dxy_deta(m,n,y,i,j,op);

      J[j][i] = 1./(dx_dksi*dy_deta - dx_deta*dy_dksi);

      ksix = J[j][i]*dy_deta;
      ksiy = -J[j][i]*dx_deta;
      etax = -J[j][i]*dy_dksi;
      etay = J[j][i]*dx_dksi;

      A1[j][i] = ksix*ksix + ksiy*ksiy;
      A2[j][i] = ksix*etax + ksiy*etay;
      A3[j][i] = etax*etax + etay*etay;
    }

    J[j][n-1] = J[j][0];
    A1[j][n-1] = A1[j][0];
    A2[j][n-1] = A2[j][0];
    A3[j][n-1] = A3[j][0];
  }
}

// /*
// - gigiaero, 22/06/2026, 2127 hours
// */
// double calc_A_metrics(int m,int n,double x[m][n],double y[m][n],int i,int j,
//                       int A_num,int op){
//   double ksix,ksiy;
//   double etax,etay;
//   double J_val = calc_J(m,n,x,y,i,j,op);

//   switch(A_num){
//     case 1: // A1
//       ksix = J_val*get_dxy_deta(m,n,y,i,j,op);
//       ksiy = -J_val*get_dxy_deta(m,n,x,i,j,op);

//       return ksix*ksix + ksiy*ksiy;

//       break;

//     case 2: // A2
//       ksix = J_val*get_dxy_deta(m,n,y,i,j,op);
//       ksiy = -J_val*get_dxy_deta(m,n,x,i,j,op);
//       etax = -J_val*get_dxy_dksi(m,n,y,i,j,op);
//       etay = J_val*get_dxy_dksi(m,n,x,i,j,op);

//       return ksix*etax + ksiy*etay;

//       break;

//     case 3: // A3
//       etax = -J_val*get_dxy_dksi(m,n,y,i,j,op);
//       etay = J_val*get_dxy_dksi(m,n,x,i,j,op);

//       return etax*etax + etay*etay;

//       break;
      
//     default:
//       puts("calc_A1: invalid A_num");
//       exit(15);
//   }
// }

/*
- gigiaero, 01/06/2026, 0956 hours

refactored
- gigiaero, 22/06/2026, 2217 hours

refactored, again
- gigiaero, 23/06/2026, 0949 hours
*/
double calc_rho(int m,int n,double rho[m][n],double phi[m][n],double A2[m][n],
                double A3[m][n],double contraU[m][n],double contraV[m][n],
                int op){
  double dphi_dksi,dphi_deta;
  double rho_primitive;
  int j_end;

  if(op == 2)
    j_end = m-1;
  else
    j_end = m;

  for(int j=0;j<j_end;j++){
    for(int i=0;i<n-1;i++){
      dphi_dksi = get_dphi_dksi(m,n,phi,i,j,op);
      dphi_deta = get_dphi_deta(m,n,phi,A2[j][i],A3[j][i],dphi_dksi,i,j,op);

      rho_primitive = 1 - (gamma - 1.)/(gamma + 1.)*(contraU[j][i]*dphi_dksi + 
                      contraV[j][i]*dphi_deta);
      
      if(rho_primitive < 0.)
        rho_primitive = .01;
      
      rho[j][i] = pow(rho_primitive,1./(gamma - 1.));
    }

    rho[j][n-1] = rho[j][0];
  }

  // Boundary conditions for the i,j+1/2 case
  if(op == 2){
    for(int i=0;i<n;i++)
      rho[m-1][i] = rho[m-2][i];
  }
}

/*
- gigiaero, 29/05/2026, 1448 hours
*/
void evaluate_delta_form_fullp(int m,int n,sim_prmtrs *config,
                               fullp_prmtrs *fp_prmtrs,char *fname_msh_x,
                               char *fname_msh_y){
  double (*x)[n] = calloc(m,sizeof *x);
  double (*y)[n] = calloc(m,sizeof *y);
  double (*J)[n] = calloc(m,sizeof *J);
  double (*A1)[n] = calloc(m,sizeof *A1);
  double (*A2)[n] = calloc(m,sizeof *A2);
  double (*A3)[n] = calloc(m,sizeof *A3);
  double (*J_ih)[n] = calloc(m,sizeof *J_ih);
  double (*A1_ih)[n] = calloc(m,sizeof *A1_ih);
  double (*A2_ih)[n] = calloc(m,sizeof *A2_ih);
  double (*A3_ih)[n] = calloc(m,sizeof *A3_ih);
  double (*J_jh)[n] = calloc(m,sizeof *J_jh);
  double (*A1_jh)[n] = calloc(m,sizeof *A1_jh);
  double (*A2_jh)[n] = calloc(m,sizeof *A2_jh);
  double (*A3_jh)[n] = calloc(m,sizeof *A3_jh);
  double (*phi)[n] = calloc(m,sizeof *phi);

  save_prmtrs_sim(config);
  save_prmtrs_fullp(config->casename,fp_prmtrs,fname_msh_x,fname_msh_y);

  read_2d_array_from_file(m,n,x,fname_msh_x);
  read_2d_array_from_file(m,n,y,fname_msh_y);

  // if(config->Ntype == 3){
  //   flip_2d_array(m,n,x);
  //   flip_2d_array(m,n,y);
  //   // af2_flip = 1;
  // }

  // Center of the external boundary must be positioned at the c/4 point
  if(fp_prmtrs->lift){
    for(int j=0;j<m;j++){
      for(int i=0;i<n;i++)
        x[j][i] -= .25;
    }
  }

  char *filename = malloc(sizeof(char)*200);
  sprintf(filename,"%s_x_mesh.dat",config->casename);
  print_2d_array_to_file(m,n,x,filename,0);
  sprintf(filename,"%s_y_mesh.dat",config->casename);
  print_2d_array_to_file(m,n,y,filename,0);

  if(config->Ntype == 3){
    flip_2d_array(m,n,x);
    flip_2d_array(m,n,y);
  }

  // print_2d_array_to_file(m,n,x,"mat_mesh_x.dat",0);
  // print_2d_array_to_file(m,n,y,"mat_mesh_y.dat",0);

  calc_J_A_metrics(m,n,J_ih,A1_ih,A2_ih,A3_ih,x,y,1);
  calc_J_A_metrics(m,n,J_jh,A1_jh,A2_jh,A3_jh,x,y,2);
  calc_J_A_metrics(m,n,J,A1,A2,A3,x,y,3);

  initialize_fullp(m,n,phi,x,y,fp_prmtrs);

  if(config->save_i_c){
    sprintf(filename,"%s_iter_%010d.dat",config->casename,0);
    print_2d_array_to_file(m,n,phi,filename,0);
  }

  alpha_sequence_aH(m,n,x,y,config,1);

  switch(config->Ntype){
    case 1:
      puts("fullp slor pending");
      // puts("SLOR, full potential flow over airfoil");
      // solve_slor_2d_rectangular_fullp();
      break;
    
    case 2:
      puts("fullp adi pending");
      // puts("ADI, full potential flow over airfoil");
      // solve_adi_2d_rectangular_fullp(m,n,phi,J,A1,A2,A3,config,fp_prmtrs);
      break;
    
    case 3:
      puts("AF2, full potential flow over airfoil");
      solve_af2_2d_rectangular_fullp(m,n,phi,J,A1,A2,A3,J_ih,A1_ih,A2_ih,A3_ih,
                                     J_jh,A1_jh,A2_jh,A3_jh,x[0],y[0],config,
                                     fp_prmtrs);
      break;
      
    default:
      puts("evaluate_delta_form_fullp: Invalid Ntype");
      exit(32);
  }

  double (*u)[n] = calloc(m,sizeof *u);
  double (*v)[n] = calloc(m,sizeof *v);
  double (*q)[n] = calloc(m,sizeof *q);
  double (*contraU)[n] = calloc(m,sizeof *contraU);
  double (*contraV)[n] = calloc(m,sizeof *contraV);
  double (*rho)[n] = calloc(m,sizeof *rho);
  double (*p)[n] = calloc(m,sizeof *p);
  double (*cp)[n] = calloc(m,sizeof *cp);

  calc_contraUV(m,n,contraU,phi,A1,A2,A3,3,1);
  calc_contraUV(m,n,contraV,phi,A1,A2,A3,3,2);
  calc_rho(m,n,rho,phi,A2,A3,contraU,contraV,3);
  get_u_v_fullp(m,n,u,v,phi,x,y,J,A1,A3,3);
  get_q_fullp(m,n,q,phi,contraU,contraV,A2,A3,3);
  get_p_cp_fullp(m,n,p,cp,rho,fp_prmtrs->Ma);

  if(config->Ntype == 3){
    flip_2d_array(m,n,contraU);
    flip_2d_array(m,n,contraV);
    flip_2d_array(m,n,rho);
    flip_2d_array(m,n,u);
    flip_2d_array(m,n,v);
    flip_2d_array(m,n,q);
    flip_2d_array(m,n,p);
    flip_2d_array(m,n,cp);
  }

  sprintf(filename,"%s_u.dat",config->casename);
  print_2d_array_to_file(m,n,u,filename,0);
  sprintf(filename,"%s_v.dat",config->casename);
  print_2d_array_to_file(m,n,v,filename,0);
  sprintf(filename,"%s_q.dat",config->casename);
  print_2d_array_to_file(m,n,q,filename,0);
  sprintf(filename,"%s_contraU.dat",config->casename); // é de fato útil salvar as velocidades contravariantes?
  print_2d_array_to_file(m,n,contraU,filename,0);
  sprintf(filename,"%s_contraV.dat",config->casename);
  print_2d_array_to_file(m,n,contraV,filename,0);
  sprintf(filename,"%s_rho.dat",config->casename);
  print_2d_array_to_file(m,n,rho,filename,0);
  sprintf(filename,"%s_p.dat",config->casename);
  print_2d_array_to_file(m,n,p,filename,0);
  sprintf(filename,"%s_cp.dat",config->casename);
  print_2d_array_to_file(m,n,cp,filename,0);

  // // // apagar depois
  // double (*dphi_dksi)[n] = calloc(m,sizeof *dphi_dksi);
  // for(int j=0;j<m;j++){
  //   for(int i=0;i<n;i++)
  //   dphi_dksi[j][i] = get_dphi_dksi(m,n,phi,i,j,3);
  // }
  // sprintf(filename,"mat_dphi_dksi.dat");
  // print_2d_array_to_file(m,n,dphi_dksi,filename,0);

  // sprintf(filename,"%s_A1.dat",config->casename);
  // print_2d_array_to_file(m,n,A1,filename,0);
  // sprintf(filename,"%s_A2.dat",config->casename);
  // print_2d_array_to_file(m,n,A2,filename,0);
  // sprintf(filename,"%s_A3.dat",config->casename);
  // print_2d_array_to_file(m,n,A3,filename,0);
  // sprintf(filename,"%s_J.dat",config->casename);
  // print_2d_array_to_file(m,n,J,filename,0);

  // double (*x_ih)[n-1] = calloc(m,sizeof *x_ih);
  // double (*y_ih)[n-1] = calloc(m,sizeof *y_ih);
  // double (*x_jh)[n] = calloc(m-1,sizeof *x_jh);
  // double (*y_jh)[n] = calloc(m-1,sizeof *y_jh);
  // get_half_meshes(m,n,x,y,x_ih,y_ih,x_jh,y_jh);

  // double (*phi_ih)[n-1] = calloc(m,sizeof *phi_ih);
  // double (*u_ih)[n-1] = calloc(m,sizeof *u_ih);
  // double (*v_ih)[n-1] = calloc(m,sizeof *v_ih);
  // double (*q_ih)[n-1] = calloc(m,sizeof *q_ih);
  // initialize_fullp(m,n-1,phi_ih,x_ih,y_ih,fp_prmtrs);
  // get_u_v_potential_fullp(m,n-1,phi_ih,x_ih,y_ih,J_ih,A2_ih,A3_ih,u_ih,v_ih,q_ih,1);

  // double (*phi_jh)[n] = calloc(m-1,sizeof *phi_jh);
  // double (*u_jh)[n] = calloc(m-1,sizeof *u_jh);
  // double (*v_jh)[n] = calloc(m-1,sizeof *v_jh);
  // double (*q_jh)[n] = calloc(m-1,sizeof *q_jh);
  // initialize_fullp(m-1,n,phi_jh,x_jh,y_jh,fp_prmtrs);
  // get_u_v_potential_fullp(m-1,n,phi_jh,x_jh,y_jh,J_jh,A2_jh,A3_jh,u_jh,v_jh,q_jh,2);

  // sprintf(filename,"%s_potential_ih.dat",config->casename);
  // print_2d_array_to_file(m,n-1,phi_ih,filename,0);
  // sprintf(filename,"%s_potential_jh.dat",config->casename);
  // print_2d_array_to_file(m-1,n,phi_jh,filename,0);

  // sprintf(filename,"%s_x_mesh_ih.dat",config->casename);
  // print_2d_array_to_file(m,n-1,x_ih,filename,0);
  // sprintf(filename,"%s_y_mesh_ih.dat",config->casename);
  // print_2d_array_to_file(m,n-1,y_ih,filename,0);
  // sprintf(filename,"%s_x_mesh_jh.dat",config->casename);
  // print_2d_array_to_file(m-1,n,x_jh,filename,0);
  // sprintf(filename,"%s_y_mesh_jh.dat",config->casename);
  // print_2d_array_to_file(m-1,n,y_jh,filename,0);

  // sprintf(filename,"%s_q_ih.dat",config->casename);
  // print_2d_array_to_file(m,n-1,q_ih,filename,0);
  // sprintf(filename,"%s_q_jh.dat",config->casename);
  // print_2d_array_to_file(m-1,n,q_jh,filename,0);
  // sprintf(filename,"%s_u_ih.dat",config->casename);
  // print_2d_array_to_file(m,n-1,u_ih,filename,0);
  // sprintf(filename,"%s_u_jh.dat",config->casename);
  // print_2d_array_to_file(m-1,n,u_jh,filename,0);
  // sprintf(filename,"%s_v_ih.dat",config->casename);
  // print_2d_array_to_file(m,n-1,v_ih,filename,0);
  // sprintf(filename,"%s_v_jh.dat",config->casename);
  // print_2d_array_to_file(m-1,n,v_jh,filename,0);

  // free(dphi_dksi);
  // free(x_ih);
  // free(y_ih);
  // free(x_jh);
  // free(y_jh);
  // free(phi_ih);
  // free(u_ih);
  // free(v_ih);
  // free(q_ih);
  // free(phi_jh);
  // free(u_jh);
  // free(v_jh);
  // free(q_jh);


  free(filename);
  free(x);
  free(y);
  free(J);
  free(A1);
  free(A2);
  free(A3);
  free(J_ih);
  free(A1_ih);
  free(A2_ih);
  free(A3_ih);
  free(J_jh);
  free(A1_jh);
  free(A2_jh);
  free(A3_jh);
  free(phi);
  free(u);
  free(v);
  free(q);
  free(contraU);
  free(contraV);
  free(rho);
  free(p);
  free(cp);
}

// stopped for now. this code is getting to ugly for my taste and i should place
// my focus elsewhere
// - gigiaero, 15/07/2026, 1408 hours
// void fullp_multigrid(int m[4],int n[4],sim_prmtrs *config,
//                      fullp_prmtrs *fp_prmtrs,char *fname_msh_x[],
//                      char *fname_msh_y[]){
//   // Prepare meshes
//   double (*x1)[n[0]] = calloc(m[0],sizeof *x1);
//   double (*y1)[n[0]] = calloc(m[0],sizeof *y1);
//   double (*x2)[n[1]] = calloc(m[1],sizeof *x2);
//   double (*y2)[n[1]] = calloc(m[1],sizeof *y2);
//   double (*x3)[n[2]] = calloc(m[2],sizeof *x3);
//   double (*y3)[n[2]] = calloc(m[2],sizeof *y3);
//   double (*x4)[n[3]] = calloc(m[3],sizeof *x4);
//   double (*y4)[n[3]] = calloc(m[3],sizeof *y4);

//   save_prmtrs_sim(config);
//   save_prmtrs_fullp_multigrid(config->casename,fp_prmtrs,fname_msh_x,
//                               fname_msh_y);

//   read_2d_array_from_file(m[0],n[0],x1,fname_msh_x[0]);
//   read_2d_array_from_file(m[0],n[0],y1,fname_msh_y[0]);
//   read_2d_array_from_file(m[1],n[1],x2,fname_msh_x[1]);
//   read_2d_array_from_file(m[1],n[1],y2,fname_msh_y[1]);
//   read_2d_array_from_file(m[2],n[2],x3,fname_msh_x[2]);
//   read_2d_array_from_file(m[2],n[2],y3,fname_msh_y[2]);
//   read_2d_array_from_file(m[3],n[3],x4,fname_msh_x[3]);
//   read_2d_array_from_file(m[3],n[3],y4,fname_msh_y[3]);

//   char *filename = malloc(sizeof(char)*200);
//   sprintf(filename,"%s_x_mesh1.dat",config->casename);
//   print_2d_array_to_file(m[0],n[0],x1,filename,0);
//   sprintf(filename,"%s_y_mesh1.dat",config->casename);
//   print_2d_array_to_file(m[0],n[0],y1,filename,0);
//   sprintf(filename,"%s_x_mesh2.dat",config->casename);
//   print_2d_array_to_file(m[1],n[1],x2,filename,0);
//   sprintf(filename,"%s_y_mesh2.dat",config->casename);
//   print_2d_array_to_file(m[1],n[1],y2,filename,0);
//   sprintf(filename,"%s_x_mesh3.dat",config->casename);
//   print_2d_array_to_file(m[2],n[2],x3,filename,0);
//   sprintf(filename,"%s_y_mesh3.dat",config->casename);
//   print_2d_array_to_file(m[2],n[2],y3,filename,0);
//   sprintf(filename,"%s_x_mesh4.dat",config->casename);
//   print_2d_array_to_file(m[3],n[3],x4,filename,0);
//   sprintf(filename,"%s_y_mesh4.dat",config->casename);
//   print_2d_array_to_file(m[3],n[3],y4,filename,0);

//   flip_2d_array(m[0],n[0],x1);
//   flip_2d_array(m[0],n[0],y1);
//   flip_2d_array(m[1],n[1],x2);
//   flip_2d_array(m[1],n[1],y2);
//   flip_2d_array(m[2],n[2],x3);
//   flip_2d_array(m[2],n[2],y3);
//   flip_2d_array(m[3],n[3],x4);
//   flip_2d_array(m[3],n[3],y4);

//   // Prepare general curvilinear coordinate parameters
//   double (*J_1)[n[0]] = calloc(m[0],sizeof *J_1); // Mesh 1
//   double (*A1_1)[n[0]] = calloc(m[0],sizeof *A1_1);
//   double (*A2_1)[n[0]] = calloc(m[0],sizeof *A2_1);
//   double (*A3_1)[n[0]] = calloc(m[0],sizeof *A3_1);
//   double (*J_ih_1)[n[0]] = calloc(m[0],sizeof *J_ih_1);
//   double (*A1_ih_1)[n[0]] = calloc(m[0],sizeof *A1_ih_1);
//   double (*A2_ih_1)[n[0]] = calloc(m[0],sizeof *A2_ih_1);
//   double (*A3_ih_1)[n[0]] = calloc(m[0],sizeof *A3_ih_1);
//   double (*J_jh_1)[n[0]] = calloc(m[0],sizeof *J_jh_1);
//   double (*A1_jh_1)[n[0]] = calloc(m[0],sizeof *A1_jh_1);
//   double (*A2_jh_1)[n[0]] = calloc(m[0],sizeof *A2_jh_1);
//   double (*A3_jh_1)[n[0]] = calloc(m[0],sizeof *A3_jh_1);
//   double (*phi_1)[n[0]] = calloc(m[0],sizeof *phi_1);
//   double (*J_2)[n[1]] = calloc(m[1],sizeof *J_2); // Mesh 2
//   double (*A1_2)[n[1]] = calloc(m[1],sizeof *A1_2);
//   double (*A2_2)[n[1]] = calloc(m[1],sizeof *A2_2);
//   double (*A3_2)[n[1]] = calloc(m[1],sizeof *A3_2);
//   double (*J_ih_2)[n[1]] = calloc(m[1],sizeof *J_ih_2);
//   double (*A1_ih_2)[n[1]] = calloc(m[1],sizeof *A1_ih_2);
//   double (*A2_ih_2)[n[1]] = calloc(m[1],sizeof *A2_ih_2);
//   double (*A3_ih_2)[n[1]] = calloc(m[1],sizeof *A3_ih_2);
//   double (*J_jh_2)[n[1]] = calloc(m[1],sizeof *J_jh_2);
//   double (*A1_jh_2)[n[1]] = calloc(m[1],sizeof *A1_jh_2);
//   double (*A2_jh_2)[n[1]] = calloc(m[1],sizeof *A2_jh_2);
//   double (*A3_jh_2)[n[1]] = calloc(m[1],sizeof *A3_jh_2);
//   double (*phi_2)[n[1]] = calloc(m[1],sizeof *phi_2);
//   double (*J_3)[n[2]] = calloc(m[2],sizeof *J_3); // Mesh 3
//   double (*A1_3)[n[2]] = calloc(m[2],sizeof *A1_3);
//   double (*A2_3)[n[2]] = calloc(m[2],sizeof *A2_3);
//   double (*A3_3)[n[2]] = calloc(m[2],sizeof *A3_3);
//   double (*J_ih_3)[n[2]] = calloc(m[2],sizeof *J_ih_3);
//   double (*A1_ih_3)[n[2]] = calloc(m[2],sizeof *A1_ih_3);
//   double (*A2_ih_3)[n[2]] = calloc(m[2],sizeof *A2_ih_3);
//   double (*A3_ih_3)[n[2]] = calloc(m[2],sizeof *A3_ih_3);
//   double (*J_jh_3)[n[2]] = calloc(m[2],sizeof *J_jh_3);
//   double (*A1_jh_3)[n[2]] = calloc(m[2],sizeof *A1_jh_3);
//   double (*A2_jh_3)[n[2]] = calloc(m[2],sizeof *A2_jh_3);
//   double (*A3_jh_3)[n[2]] = calloc(m[2],sizeof *A3_jh_3);
//   double (*phi_3)[n[2]] = calloc(m[2],sizeof *phi_3);
//   double (*J_4)[n[3]] = calloc(m[3],sizeof *J_4); // Mesh 4
//   double (*A1_4)[n[3]] = calloc(m[3],sizeof *A1_4);
//   double (*A2_4)[n[3]] = calloc(m[3],sizeof *A2_4);
//   double (*A3_4)[n[3]] = calloc(m[3],sizeof *A3_4);
//   double (*J_ih_4)[n[3]] = calloc(m[3],sizeof *J_ih_4);
//   double (*A1_ih_4)[n[3]] = calloc(m[3],sizeof *A1_ih_4);
//   double (*A2_ih_4)[n[3]] = calloc(m[3],sizeof *A2_ih_4);
//   double (*A3_ih_4)[n[3]] = calloc(m[3],sizeof *A3_ih_4);
//   double (*J_jh_4)[n[3]] = calloc(m[3],sizeof *J_jh_4);
//   double (*A1_jh_4)[n[3]] = calloc(m[3],sizeof *A1_jh_4);
//   double (*A2_jh_4)[n[3]] = calloc(m[3],sizeof *A2_jh_4);
//   double (*A3_jh_4)[n[3]] = calloc(m[3],sizeof *A3_jh_4);
//   double (*phi_4)[n[3]] = calloc(m[3],sizeof *phi_4);

//   calc_J_A_metrics(m[0],n[0],J_ih_1,A1_ih_1,A2_ih_1,A3_ih_1,x1,y1,1);
//   calc_J_A_metrics(m[0],n[0],J_jh_1,A1_jh_1,A2_jh_1,A3_jh_1,x1,y1,2);
//   calc_J_A_metrics(m[0],n[0],J_1,A1_1,A2_1,A3_1,x1,y1,3);
//   calc_J_A_metrics(m[1],n[1],J_ih_2,A1_ih_2,A2_ih_2,A3_ih_2,x2,y2,1);
//   calc_J_A_metrics(m[1],n[1],J_jh_2,A1_jh_2,A2_jh_2,A3_jh_2,x2,y2,2);
//   calc_J_A_metrics(m[1],n[1],J_2,A1_2,A2_2,A3_2,x2,y2,3);
//   calc_J_A_metrics(m[2],n[2],J_ih_3,A1_ih_3,A2_ih_3,A3_ih_3,x3,y3,1);
//   calc_J_A_metrics(m[2],n[2],J_jh_3,A1_jh_3,A2_jh_3,A3_jh_3,x3,y3,2);
//   calc_J_A_metrics(m[2],n[2],J_3,A1_3,A2_3,A3_3,x3,y3,3);
//   calc_J_A_metrics(m[3],n[3],J_ih_4,A1_ih_4,A2_ih_4,A3_ih_4,x4,y4,1);
//   calc_J_A_metrics(m[3],n[3],J_jh_4,A1_jh_4,A2_jh_4,A3_jh_4,x4,y4,2);
//   calc_J_A_metrics(m[3],n[3],J_4,A1_4,A2_4,A3_4,x4,y4,3);

//   initialize_fullp(m[0],n[0],phi_1,x1,y1,fp_prmtrs);
//   initialize_fullp(m[1],n[1],phi_2,x2,y2,fp_prmtrs);
//   initialize_fullp(m[2],n[2],phi_3,x3,y3,fp_prmtrs);
//   initialize_fullp(m[3],n[3],phi_4,x4,y4,fp_prmtrs);

//   // Meshes
//   free(x1);
//   free(y1);
//   free(x2);
//   free(y2);
//   free(x3);
//   free(y3);
//   free(x4);
//   free(y4);
//   free(filename);

//   // GCC
//   free(A1_1); // Mesh 1
//   free(A2_1);
//   free(A3_1);
//   free(J_ih_1);
//   free(A1_ih_1);
//   free(A2_ih_1);
//   free(A3_ih_1);
//   free(J_jh_1);
//   free(A1_jh_1);
//   free(A2_jh_1);
//   free(A3_jh_1);
//   free(phi_1);
//   free(A1_2); // Mesh 2
//   free(A2_2);
//   free(A3_2);
//   free(J_ih_2);
//   free(A1_ih_2);
//   free(A2_ih_2);
//   free(A3_ih_2);
//   free(J_jh_2);
//   free(A1_jh_2);
//   free(A2_jh_2);
//   free(A3_jh_2);
//   free(phi_2);
//   free(A1_3); // Mesh 3
//   free(A2_3);
//   free(A3_3);
//   free(J_ih_3);
//   free(A1_ih_3);
//   free(A2_ih_3);
//   free(A3_ih_3);
//   free(J_jh_3);
//   free(A1_jh_3);
//   free(A2_jh_3);
//   free(A3_jh_3);
//   free(phi_3);
//   free(A1_4); // Mesh 4
//   free(A2_4);
//   free(A3_4);
//   free(J_ih_4);
//   free(A1_ih_4);
//   free(A2_ih_4);
//   free(A3_ih_4);
//   free(J_jh_4);
//   free(A1_jh_4);
//   free(A2_jh_4);
//   free(A3_jh_4);
//   free(phi_4);
// }

/*
- gigiaero, 01/06/2026, 1555 hours
*/
double freestream_u(double Ma){
  return sqrt((gamma + 1.)/(gamma - 1. + 2./Ma/Ma));
}

// /*
// - gigiaero, 20/06/2026, 2146 hours
// */
// //NOTA: revisar a velocidade do som aqui
// double get_mach(double contraU,double contraV,double dphi_dksi,
//                 double dphi_deta){
//   double q_sq,a_sq;

//   q_sq = contraU*dphi_dksi + contraV*dphi_deta;

//   // a_sq = .5*(gamma + 1.) - .5*q_sq*(gamma - 1.);

//   a_sq = 1. - (gamma - 1.)/(gamma + 1.)*q_sq;

//   a_sq *= .5*(gamma + 1.);

//   return sqrt(q_sq/a_sq);
// }

/*
- gigiaero, 21/06/2026, 1958 hours
*/
double get_dphi_deta(int m,int n,double phi[m][n],double A2_val,double A3_val,
                     double dphi_dksi,int i,int j,int op){
  check_j(m,j);
  check_j(n,i);
  
  switch(op){
    case 1: // i+1/2,j
      // if(j == 0 || j == m-1) // adaptation for adi?
      if(j == m-1)
        return -A2_val*dphi_dksi/A3_val;

      else if(j == 0){
        if(i == n-1)
          i = 0;

        return (uniform_scheme_der1_o2_forward(m,n,phi,i,j,2) + 
                uniform_scheme_der1_o2_forward(m,n,phi,i+1,j,2))*.5; // tem que ser de primeira ordem?
      }
      else{
        if(i == n-1)
          i = 0;

        return (uniform_scheme_der1_o2_central(m,n,phi,i,j,2) + 
                uniform_scheme_der1_o2_central(m,n,phi,i+1,j,2))*.5;
      }
      break;

    case 2: // i,j+1/2
      if(j == m-1){
        puts("get_dphi_deta: case 2, invalid j");
        exit(207);
      }

      else
        return uniform_scheme_der1_o1_forward(m,n,phi,i,j,2);

      break;

    case 3: // i,j
      // if(j == 0 || j == m-1) // adaptation for adi?
      if(j == m-1)
        return -A2_val*dphi_dksi/A3_val;
        // return uniform_scheme_der1_o1_backward(m,n,phi,i,j,2);

      else if (j == 0){
        check_j(m,j+1);
        check_j(m,j+2);
        return uniform_scheme_der1_o2_forward(m,n,phi,i,j,2);
      }
      else
        return uniform_scheme_der1_o2_central(m,n,phi,i,j,2);
      
      break;

    default:
      puts("get_dphi_deta: invalid op");
      exit(15);
  }
}

/*
- gigiaero, 21/06/2026, 1356 hours
*/
double get_dphi_dksi(int m,int n,double phi[m][n],int i,int j,int op){
  check_j(m,j);
  check_j(n,i);
        
  switch(op){
    case 1: // i+1/2,j
      if(i == n-1)
        i = 0;

      return uniform_scheme_der1_o1_forward(m,n,phi,i,j,1);

      break;

    case 2: // i,j+1/2
      if(j == m-1){ // Approximate to case 3
        if(i == 0 || i == n-1)
          return uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j);

        else
          return uniform_scheme_der1_o2_central(m,n,phi,i,j,1);
      
      }
      else{
        if(i == 0 || i == n-1)
          return (uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j) + 
                  uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j+1))*.5;

        else
          return (uniform_scheme_der1_o2_central(m,n,phi,i,j,1) + 
                  uniform_scheme_der1_o2_central(m,n,phi,i,j+1,1))*.5;
      }
      break;

    case 3: // i,j
      if(i == 0 || i == n-1)
        return uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j);

      else
        return uniform_scheme_der1_o2_central(m,n,phi,i,j,1);
      
      break;

    default:
      puts("get_dphi_dksi: invalid op");
      exit(15);
  }
}

/*
- gigiaero, 22/06/2026, 2103 hours
*/
double get_dxy_deta(int m,int n,double xy[m][n],int i,int j,int op){
  switch(op){
    case 1: // i+1/2,j
      if(i == n-1)
        i = 0;

      if(j == 0)
        return (uniform_scheme_der1_o2_forward(m,n,xy,i,j,2) + 
                uniform_scheme_der1_o2_forward(m,n,xy,i+1,j,2))*.5;
      
      else if(j == m-1)
        return (uniform_scheme_der1_o2_backward(m,n,xy,i,j,2) + 
                uniform_scheme_der1_o2_backward(m,n,xy,i+1,j,2))*.5;

      else
        return (uniform_scheme_der1_o2_central(m,n,xy,i,j,2) + 
                uniform_scheme_der1_o2_central(m,n,xy,i+1,j,2))*.5;
      
      break;

    case 2: // i,j+1/2
      if(j == m-1) // Approximate to case 3
        return uniform_scheme_der1_o2_backward(m,n,xy,i,j,2);
      
      else
        return uniform_scheme_der1_o1_forward(m,n,xy,i,j,2);

      break;

    case 3: // i,j
      if(j == 0)
        return uniform_scheme_der1_o2_forward(m,n,xy,i,j,2);

      else if(j == m-1)
        return uniform_scheme_der1_o2_backward(m,n,xy,i,j,2);

      else
        return uniform_scheme_der1_o2_central(m,n,xy,i,j,2);
      
      break;
      
    default:
      puts("get_dx_deta: invalid op");
      exit(15);
  }
}

/*
- gigiaero, 22/06/2026, 2056 hours
*/
double get_dxy_dksi(int m,int n,double xy[m][n],int i,int j,int op){
  switch(op){
    case 1: // i+1/2,j
      if(i == n-1)
        return uniform_scheme_der1_o1_forward(m,n,xy,0,j,1);

      else
        return uniform_scheme_der1_o1_forward(m,n,xy,i,j,1);

      break;

    case 2: // i,j+1/2
      if(j == m-1){ // Approximate to case 3
        if(i == 0)
          return uniform_scheme_der1_o2_central_prdc_ksi(m,n,xy,j);

        else
          return uniform_scheme_der1_o2_central(m,n,xy,i,j,1);
      }
      else{
        if(i == 0)
          return (uniform_scheme_der1_o2_central_prdc_ksi(m,n,xy,j) + 
                  uniform_scheme_der1_o2_central_prdc_ksi(m,n,xy,j+1))*.5;

        else
          return (uniform_scheme_der1_o2_central(m,n,xy,i,j,1) + 
                  uniform_scheme_der1_o2_central(m,n,xy,i,j+1,1))*.5;
      }
      break;

    case 3: // i,j
      if(i == 0)
        return uniform_scheme_der1_o2_central_prdc_ksi(m,n,xy,j);

      else
        return uniform_scheme_der1_o2_central(m,n,xy,i,j,1);

      break;
      
    default:
      puts("get_dx_dksi: invalid op");
      exit(15);
  }
}

/*
and here comes another effort in debugging
-gigiaero, 10/07/2026, 1439 hours
*/
void get_half_meshes(int m,int n,double x[m][n],double y[m][n],
                     double x_ih[m][n-1],double y_ih[m][n-1],
                     double x_jh[m-1][n],double y_jh[m-1][n]){
  for(int j=0;j<m;j++){
    for(int i=0;i<n-1;i++){
      x_ih[j][i] = (x[j][i] + x[j][i+1])*.5;
      y_ih[j][i] = (y[j][i] + y[j][i+1])*.5;
    }
  }

  for(int j=0;j<m-1;j++){
    for(int i=0;i<n;i++){
      x_jh[j][i] = (x[j][i] + x[j+1][i])*.5;
      y_jh[j][i] = (y[j][i] + y[j+1][i])*.5;
    }
  }
}

/*
- gigiaero, 12/07/2026, 2105 hours
*/
void get_p_cp_fullp(int m,int n,double p[m][n],double cp[m][n],double rho[m][n],
                    double Ma){
  double Uinf = freestream_u(Ma);
  double rho_inf =pow(1. - (gamma - 1.)/(gamma + 1.)*Uinf*Uinf,1./(gamma - 1.));
  double p_inf = (gamma + 1.)*pow(rho_inf,gamma)/2./gamma;

  for(int j=0;j<m;j++){
    for(int i=0;i<n-1;i++){
      p[j][i] = (gamma + 1.)*pow(rho[j][i],gamma)/2./gamma;
      cp[j][i] = (p[j][i] - p_inf)/(.5*rho_inf*Uinf*Uinf);
    }

    p[j][n-1] = p[j][0];
    cp[j][n-1] = cp[j][0];
  }
}

/*
- gigiaero, 12/07/2026, 2013 hours
*/
void get_q_fullp(int m,int n,double q[m][n],double phi[m][n],
                 double contraU[m][n],double contraV[m][n],double A2[m][n],
                 double A3[m][n],int op){
  double dphi_dksi,dphi_deta;
                  
  for(int j=0;j<m;j++){
    for(int i=0;i<n-1;i++){
      dphi_dksi = get_dphi_dksi(m,n,phi,i,j,op);
      dphi_deta = get_dphi_deta(m,n,phi,A2[j][i],A3[j][i],dphi_dksi,i,j,op);

      q[j][i] = sqrt(contraU[j][i]*dphi_dksi + contraV[j][i]*dphi_deta);
    }

    q[j][n-1] = q[j][0];
  }
}

/*
- gigiaero, 01/06/2026, 2204 hours

q calculations moved to a dedicated function
- gigiaero, 12/07/2026, 2010 hours
*/
void get_u_v_fullp(int m,int n,double u[m][n],double v[m][n],double phi[m][n],
                   double x[m][n],double y[m][n],double J[m][n],double A2[m][n],
                   double A3[m][n],int op){
  double dx_dksi,dx_deta;
  double dy_dksi,dy_deta;
  double ksix,ksiy;
  double etax,etay;
  double dphi_dksi,dphi_deta;

  for(int j=0;j<m;j++){
    for(int i=0;i<n-1;i++){
      dx_dksi = get_dxy_dksi(m,n,x,i,j,op);
      dx_deta = get_dxy_deta(m,n,x,i,j,op);
      dy_dksi = get_dxy_dksi(m,n,y,i,j,op);
      dy_deta = get_dxy_deta(m,n,y,i,j,op);

      ksix = J[j][i]*dy_deta;
      ksiy = -J[j][i]*dx_deta;
      etax = -J[j][i]*dy_dksi;
      etay = J[j][i]*dx_dksi;

      dphi_dksi = get_dphi_dksi(m,n,phi,i,j,op);
      dphi_deta = get_dphi_deta(m,n,phi,A2[j][i],A3[j][i],dphi_dksi,i,j,op);

      // if(i == 0 || i == n-1)
      //   dphi_dksi = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j);

      // else
      //   dphi_dksi = uniform_scheme_der1_o2_central(m,n,phi,i,j,1);

      // if(j == 0)
      //   dphi_deta = uniform_scheme_der1_o2_forward(m,n,phi,i,j,2);

      // else if(j == m-1)
      //   dphi_deta = uniform_scheme_der1_o2_backward(m,n,phi,i,j,2);

      // else
      //   dphi_deta = uniform_scheme_der1_o2_central(m,n,phi,i,j,2);

      u[j][i] = dphi_dksi*ksix + dphi_deta*etax;
      v[j][i] = dphi_dksi*ksiy + dphi_deta*etay;

      // q[j][i] = sqrt(pow(u[j][i],2) + pow(v[j][i],2));
    }

    u[j][n-1] = u[j][0];
    v[j][n-1] = v[j][0];
    // q[j][n-1] = q[j][0];
  }

  // // Domain interior
  // for(int j=1;j<m-1;j++){
  //   metric_terms(&ksix,&ksiy,&etax,&etay,m,n,J,x,y,0,j);
  //   u[j][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j)*ksix + 
  //             uniform_scheme_der1_o2_central(m,n,phi,0,j,2)*etax;
  //   v[j][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j)*ksiy + 
  //             uniform_scheme_der1_o2_central(m,n,phi,0,j,2)*etay;
    
  //   for(int i=1;i<n-1;i++){
  //     metric_terms(&ksix,&ksiy,&etax,&etay,m,n,J,x,y,i,j);
  //     u[j][i] = uniform_scheme_der1_o2_central(m,n,phi,i,j,1)*ksix + 
  //               uniform_scheme_der1_o2_central(m,n,phi,i,j,2)*etax;
  //     v[j][i] = uniform_scheme_der1_o2_central(m,n,phi,i,j,1)*ksiy + 
  //               uniform_scheme_der1_o2_central(m,n,phi,i,j,2)*etay;
  //   }
  // }

  // // External boundary and airfoil surface
  // for(int i=1;i<n-1;i++){
  //   metric_terms(&ksix,&ksiy,&etax,&etay,m,n,J,x,y,i,0);
  //   u[0][i] = uniform_scheme_der1_o2_central(m,n,phi,i,0,1)*ksix + 
  //             uniform_scheme_der1_o2_forward(m,n,phi,i,0,2)*etax;
  //   v[0][i] = uniform_scheme_der1_o2_central(m,n,phi,i,0,1)*ksiy + 
  //             uniform_scheme_der1_o2_forward(m,n,phi,i,0,2)*etay;

  //   metric_terms(&ksix,&ksiy,&etax,&etay,m,n,J,x,y,i,m-1);
  //   u[m-1][i] = uniform_scheme_der1_o2_central(m,n,phi,i,m-1,1)*ksix + 
  //               uniform_scheme_der1_o2_backward(m,n,phi,i,m-1,2)*etax;
  //   v[m-1][i] = uniform_scheme_der1_o2_central(m,n,phi,i,m-1,1)*ksiy + 
  //               uniform_scheme_der1_o2_backward(m,n,phi,i,m-1,2)*etay;
  // }

  // // Corners
  // metric_terms(&ksix,&ksiy,&etax,&etay,m,n,J,x,y,0,0);
  // u[0][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,0)*ksix + 
  //           uniform_scheme_der1_o2_forward(m,n,phi,0,0,2)*etax;
  // v[0][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,0)*ksiy + 
  //           uniform_scheme_der1_o2_forward(m,n,phi,0,0,2)*etay;

  // metric_terms(&ksix,&ksiy,&etax,&etay,m,n,J,x,y,0,m-1);
  // u[m-1][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,m-1)*ksix + 
  //             uniform_scheme_der1_o2_backward(m,n,phi,0,m-1,2)*etax;
  // v[m-1][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,m-1)*ksiy + 
  //             uniform_scheme_der1_o2_backward(m,n,phi,0,m-1,2)*etay;

  // // Velocity resultant
  // for(int j=0;j<m;j++){
  //   for(int i=0;i<n-1;i++)
  //     Ve[j][i] = sqrt(pow(u[j][i],2) + pow(v[j][i],2));

  //   u[j][n-1] = u[j][0];
  //   v[j][n-1] = v[j][0];
  //   Ve[j][n-1] = Ve[j][0];
  // }
}

/*
- gigiaero, 0/06/2026, 2057 hours
*/
void initialize_fullp(int m,int n,double phi[m][n],double x[m][n],double y[m][n],
                      fullp_prmtrs *fp_prmtrs){
  double Uinf = freestream_u(fp_prmtrs->Ma);

  fp_prmtrs->alpha *= pi/180.;

  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      phi[j][i] = Uinf*cos(fp_prmtrs->alpha)*x[j][i] + 
                  Uinf*sin(fp_prmtrs->alpha)*y[j][i];
  }
}

/*
- gigiaero, 01/06/2026,1315 hours

refactored
- gigiaero, 22/06/2026, 2252 hours
*/
void L_phi_fullp(int m,int n,double L_phi[m][n],double L_phi_terms_ih[m][n-1],
                 double L_phi_terms_jh[m][n-1]){
  // for(int j=0+af2_flip;j<m-1+af2_flip;j++){
  for(int j=1;j<m;j++){
    L_phi[j][0] = (L_phi_terms_ih[j][0] - L_phi_terms_ih[j][n-2]) + 
                  (L_phi_terms_jh[j][0] - L_phi_terms_jh[j-1][0]);

    for(int i=1;i<n-1;i++)
      L_phi[j][i] = (L_phi_terms_ih[j][i] - L_phi_terms_ih[j][i-1]) + 
                    (L_phi_terms_jh[j][i] - L_phi_terms_jh[j-1][i]);

    L_phi[j][n-1] = L_phi[j][0];
  }
}

/*
- gigiaero, 23/06/2026, 1422 hours
*/
void L_phi_fullp_der_terms_ih(int m,int n,double L_phi_terms_ih[m][n-1],
                              double rho_til[m][n],double contraU_ih[m][n],
                              double J_ih[m][n]){
  // for(int j=0+af2_flip;j<m-1+af2_flip;j++){
  for(int j=1;j<m;j++){
    for(int i=0;i<n-1;i++)
      L_phi_terms_ih[j][i] = rho_til[j][i]*contraU_ih[j][i]/J_ih[j][i];
  }
}

/*
adapted to af2 only for now
- gigiaero, 23/06/2026, 1434 hours
*/
void L_phi_fullp_der_terms_jh(int m,int n,double L_phi_terms_jh[m][n-1],
                              double rho_bar[m][n],double contraV_jh[m][n],
                              double J_jh[m][n]){
  for(int i=0;i<n-1;i++){
    for(int j=0;j<m-1;j++)
      L_phi_terms_jh[j][i] = rho_bar[j][i]*contraV_jh[j][i]/J_jh[j][i];
    
    L_phi_terms_jh[m-1][i] = -L_phi_terms_jh[m-2][i];
  }
}

/*
- gigiaero, 31/05/2026, 1326 hours
*/
double max2(double n1,double n2){
  if(n1 >= n2){
    // puts("zero");
    return n1;}
  else{
    // puts("nonzero");
    return n2;}
}

/*
- gigiaero, 31/05/2026, 1434 hours
*/
double mean4_j(int m,int n,double rho[m][n-1],int i,int j){
  // check_j(m,j-1);
  check_j(m,j);
  check_j(m,j+1);
  check_j(n,i);
  
  if(i == 0)
    return (rho[j][i] + rho[j][n-2] + rho[j+1][n-2] + rho[j+1][i])*.25;
  else{
    check_j(n,i-1);
    return (rho[j][i] + rho[j][i-1] + rho[j+1][i-1] + rho[j+1][i])*.25;}
}

// /*
// - gigiaero, 02/06/2026, 1024 hours
// */
// void metric_terms(double *ksix,double *ksiy,double *etax,double *etay,int m,
//                   int n,double J[m][n],double x[m][n],double y[m][n],
//                   int i,int j){
//   double dy_deta,dx_deta;
//   double dy_dksi,dx_dksi;

//   if(i == 0){
//     dy_dksi = uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j);
//     dx_dksi = uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j);
//   }
//   else{
//     dy_dksi = uniform_scheme_der1_o2_central(m,n,y,i,j,1);
//     dx_dksi = uniform_scheme_der1_o2_central(m,n,x,i,j,1);
//   }

//   if(j == 0){
//     // dy_deta = uniform_scheme_der1_o1_forward(m,n,y,i,j,2);
//     // dx_deta = uniform_scheme_der1_o1_forward(m,n,x,i,j,2);
//     dy_deta = uniform_scheme_der1_o2_forward(m,n,y,i,j,2);
//     dx_deta = uniform_scheme_der1_o2_forward(m,n,x,i,j,2);
//   }
//   else if(j == m-1){
//     // dy_deta = uniform_scheme_der1_o1_backward(m,n,y,i,j,2);
//     // dx_deta = uniform_scheme_der1_o1_backward(m,n,x,i,j,2);
//     dy_deta = uniform_scheme_der1_o2_backward(m,n,y,i,j,2);
//     dx_deta = uniform_scheme_der1_o2_backward(m,n,x,i,j,2);
//   }
//   else{
//     dy_deta = uniform_scheme_der1_o2_central(m,n,y,i,j,2);
//     dx_deta = uniform_scheme_der1_o2_central(m,n,x,i,j,2);
//   }

//   *ksix = J[j][i]*dy_deta;
//   *ksiy = -J[j][i]*dx_deta;
//   *etax = -J[j][i]*dy_dksi;
//   *etay = J[j][i]*dx_dksi;
// }

/*
- gigiaero, 04/06/2026, 1551 hours
*/
double norm_L2(int m,int n,double L_phi[m][n]){
  double sum = 0.;

  for(int j=0;j<m;j++){
    for(int i=0;i<n-1;i++)
      sum += L_phi[j][i]*L_phi[j][i];
  }

  return sqrt(sum);
}

/*
- gigiaero, 01/06/2026, 1048 hours

refactored 
- gigiaero, 22/06/2026, 2227 hours

refactored once more
- gigiaero, 23/06/2026, 1054 hours
*/
void nu_switch(int m,int n,double nu[m][n],double contraUV[m][n],
               double rho[m][n],double C,int axis){
  double C1 = .633939;
  double C2 = 4.9325;
  double rho_val;

  for(int j=0;j<m;j++){
    for(int i=0;i<n-1;i++){
      check_j(m,j);
      check_j(n,i);

      if(contraUV[j][i] >= 0)
        rho_val = rho[j][i];

      else{
        switch(axis){
          case 1:
            // check_j(n,i+1);
            rho_val = rho[j][i+1];
            break;

          case 2:
            // check_j(m,j+1);
            if(j+1 > m-1)
              rho_val = rho[m-1][i];
            
            else
              rho_val = rho[j+1][i];

            break;

          default:
            puts("nu_switch: invalid axis");
            exit(15);
        }
      }

      nu[j][i] = max2(0.,(C1 - rho_val)*C2*C);

      if(nu[j][i] > 1.)
        nu[j][i] = 1.;
    }

    nu[j][n-1] = nu[j][0];
  }
}

/*
- gigiaero, 31/05/2026, 1444 hours

refactored
- gigiaero, 22/06/2026, 2244 hours

refactored
- gigiaeor, 23/06/2026, 1343 hours
*/
double rho_coeffs(int m,int n,double rho_tb[m][n],double rho[m][n],
                  double contraUV[m][n],double nu[m][n],int op,int axis){
  int i_p_r,j_p_s;
  int j_end;
  int rs;

  if(op == 2)
    j_end = m-1;
  else
    j_end = m;

  for(int j=0;j<j_end;j++){
    for(int i=0;i<n-1;i++){
      check_j(m,j);
      check_j(n,i);

      rs = rs_idx(contraUV[j][i]);

      switch(axis){
        case 1: // rho_til
          if(i+rs < 0)
            i_p_r = n - 2;

          else if(i+rs > n-2)
            i_p_r = 0;

          else
            i_p_r = i + rs;

          // check_j(n,i_p_r);
          // check_j(n,i);
          
          rho_tb[j][i] = (1. - nu[j][i])*rho[j][i] + nu[j][i]*rho[j][i_p_r];

          break;

        case 2: // rho_bar
          if(j+rs > m-1) // af2_flip? condição de contorno, página 52
            j_p_s = m-1;
            
          else if(j+rs < 0)
            j_p_s = 0;
          
          else
            j_p_s = j + rs;
            
          check_j(m,j_p_s);
          check_j(m,j);

          rho_tb[j][i] = (1. - nu[j][i])*rho[j][i] + nu[j][i]*rho[j_p_s][i];

          break;

        default:
          puts("rho_coeffs: invalid axis");
          exit(15);
      }
    }

    rho_tb[j][n-1] = rho_tb[j][0];
  }

  // Boundary conditions for the i,j+1/2 case
  if(op == 2){
    for(int i=0;i<n;i++)
      rho_tb[m-1][i] = rho_tb[m-2][i];
  }
}

/*
- gigiaero, 31/05/2026, 1406 hours
*/
int rs_idx(double contraUV){
  if(contraUV < 0)
    return +1;
  else
    return -1;
}

/*
- gigiaero, 15/07/2026, ???? hours
*/
void save_prmtrs_fullp_multigrid(char *casename,fullp_prmtrs *fp_prmtrs,
                                 char *fname_msh_x[],char *fname_msh_y[]){
  char *filename = malloc(sizeof(char)*200);
  FILE *output;

  sprintf(filename,"%s_prmtrs_fullp.dat",casename);

  output = fopen(filename,"w");

  fprintf(output,"fp_prmtrs.alpha = %f;\n",fp_prmtrs->alpha);
  fprintf(output,"fp_prmtrs.Ma = %f;\n",fp_prmtrs->Ma);
  fprintf(output,"fp_prmtrs.C = %f;\n",fp_prmtrs->C);
  fprintf(output,"fp_prmtrs.beta_sub = %f;\n",fp_prmtrs->beta_sub);
  fprintf(output,"fp_prmtrs.beta_super = %f;\n",fp_prmtrs->beta_super);
  fprintf(output,"fp_prmtrs.lift = %d;\n",fp_prmtrs->lift);
  fprintf(output,"fp_prmtrs.Rg = %f;\n\n",fp_prmtrs->Rg);

  fprintf(output,"char *fname_msh_x[] = {\"%s\",\n\"%s\",\n\"%s\",\n\"%s\"};\n",
          fname_msh_x[0],fname_msh_x[1],fname_msh_x[2],fname_msh_x[3]);
  fprintf(output,"char *fname_msh_y[] = {\"%s\",\n\"%s\",\n\"%s\",\n\"%s\"};\n",
          fname_msh_y[0],fname_msh_y[1],fname_msh_y[2],fname_msh_y[3]);

  fclose(output);
  free(filename);
}

// (adi here)
/*
// void solve_adi_2d_rectangular_fullp(int m,int n,double phi[m][n],double J[m][n],
//                                     double A1[m][n],double A2[m][n],
//                                     double A3[m][n],sim_prmtrs *config,
//                                     fullp_prmtrs *fp_prmtrs){
//   // Solver variables
//   double (*L_phi)[n] = calloc(m,sizeof *L_phi);
//   double (*Cij)[n] = calloc(m,sizeof *Cij);
//   double (*rho)[n-1] = calloc(m,sizeof *rho);
//   double (*f)[n] = calloc(m,sizeof *f);
//   double *ak = malloc(sizeof(double)*(n-1));
//   double *bk = malloc(sizeof(double)*(n-1));
//   double *ck = malloc(sizeof(double)*(n-1));
//   double *fk = malloc(sizeof(double)*(n-1));
//   double *uk = malloc(sizeof(double)*(n-1));
//   double *an = malloc(sizeof(double)*(m-3));
//   double *bn = malloc(sizeof(double)*(m-2));
//   double *cn = malloc(sizeof(double)*(m-3));
//   double *fn = malloc(sizeof(double)*(m-2));
//   double *un = malloc(sizeof(double)*(m-2));
//   double Ai,Aip1,Aj,Ajp1;
//   double dphi_dksi,dphi_deta;
//   double contraU,contraV;
//   double beta_sub = fp_prmtrs->beta_sub;
//   double beta_super;
//   double beta;
//   fp_beta_prmtrs fpb_prmtrs;
//   // double beta,beta_n;
//   double omega = config->w;
//   double alpha;
//   int k = 1;
//   int iter = 1;
//   // Save files
//   char *filename_save = malloc(sizeof(char)*200);
//   char *buffer = malloc(sizeof(char)*200);
//   int str_end_idx;
//   // Residuals
//   char *filename_log = malloc(sizeof(char)*200);
//   double res,L2_res;
//   FILE *file_log;

//   // Configure log files
//   sprintf(filename_log,"%s.log",config->casename);
//   file_log = fopen(filename_log,"w");

//   // Prepare string to save simulation data  
//   sprintf(filename_save,"%s_iter_",config->casename);
//   find_str_end(filename_save,&str_end_idx);

//   // a se pensar: calc_contraU e calc_contraV poderiam ser a mesma função
  
//   beta_initial_config(config,&fpb_prmtrs,fp_prmtrs->beta_super);

//   for(iter;iter<=config->max_iter;iter++){         
//     alpha_sequence(&alpha,&k,iter,config);

//     calc_rho(m,n,rho,phi,A1,A2,A3);
//     // print_2d_array_to_file(m,n-1,rho,"mat_rho.dat",0);               
//     L_phi_fullp(m,n,L_phi,phi,J,A1,A2,A3,rho,fp_prmtrs->C);

//     // print_2d_array_to_file(m,n,L_phi,"mat_L_phi.dat",0);

//     res = 0.;
//     for(int j=0;j<m-1;j++){
//       for(int i=0;i<n-1;i++){
//         if(fabs(L_phi[j][i]) > res)
//           res = fabs(L_phi[j][i]);
//       }
//     }

//     L2_res = norm_L2(m,n,L_phi);

//     printf("ADI Iteration %010d | Res %.6e\n",iter,res);

//     fprintf(file_log,"%.6e\n",res);

//     // Test for convergence
//     if(res <= config->eps && iter != 0){
//       puts("<< Convergence! >>");
//       iter++;
//       break;
//     }

//     if(res >= div_ref){
//       puts("- Divergence");
//       iter++;
//       break;
//     }

//     // Solve for corrections - step 1 (ksi)
//     for(int j=0;j<m-1;j++){
//       Ai = calc_Ai(m,n,phi,J,A1,A2,A3,rho,fp_prmtrs->C,-1,j);
//       for(int i=0;i<n-1;i++){
//         beta_update(&beta_super,L2_res,res,config->M,&fpb_prmtrs,iter);
//         // beta_switch(&beta,beta_sub,beta_super,supersonic????); // pendente: calcular velocidade resultante por meio de u e v

//         // Ai = calc_Ai(m,n,phi,J,A1,A2,A3,rho,fp_prmtrs->C,i-1,j);
//         Aip1 = calc_Ai(m,n,phi,J,A1,A2,A3,rho,fp_prmtrs->C,i,j);

//         if(i == 0)
//           contraU = calc_contraU(uniform_scheme_der1_o2_central_prdc_ksi(m,n,
//                                                                          phi,j),
//                                  uniform_scheme_der1_o2_central(m,n,phi,i,j,2),
//                                  A1[j][i],A2[j][i]);

//         else
//           contraU = calc_contraU(uniform_scheme_der1_o2_central(m,n,phi,i,j,1),
//                                  uniform_scheme_der1_o2_central(m,n,phi,i,j,2),
//                                  A1[j][i],A2[j][i]);

//         if(contraU > 0){
//           ak[i] = -Ai - beta*fabs(contraU);
//           ck[i] = -Aip1;
//         }
//         else{
//           ak[i] = -Ai;
//           ck[i] = -Aip1 - beta*fabs(contraU);
//         }

//         bk[i] = alpha + Ai + Aip1 + beta*fabs(contraU);

//         fk[i] = omega*alpha*L_phi[j][i];

//         Ai = Aip1;
//       }

//       tridiagonal_pmatrix_solver(n-1,ak,bk,ck,fk,uk);

//       for(int i=0;i<n-1;i++)
//         f[j][i] = uk[i]; 
      
//       f[j][n-1] = f[j][0];
//     }

//     print_2d_array_to_file(m,n,f,"mat_f.dat",0);
//     exit(0);

//     // Solve for corrections - step 2 (eta)
//     for(int i=0;i<n-1;i++){
//       bn[0] =
//       cn[0] =

//       fn[0] = f[1][i];
      
//       for(int j=2;j<m-2;j++){
//         an[j-2] = 
//         bn[j-1] = 
//         cn[j-1] = 

//         fn[j-1] = f[j][i];
//       }

//       an[m-4] = 
//       bn[m-3] = 

//       fn[m-3] = f[m-2][i];

//       tridiagonal_matrix_solver(m-2,an,bn,cn,fn,un);

//       for(int j=1;j<m-1;j++)
//         Cij[j][i] = un[j-1];
//     }

//     // operador N
//     // sobre o amortecimento artificial:
//     // contraU > 0 -> backward
//     // contraU < 0 -> forward

//     for(int j=0;j<m-1;j++){
//       for(int i=0;i<n-1;i++)
//         phi[j][i] += Cij[j][i];
//     }

//     // atualizar condições de contorno

//     // forçar periodicidade aqui
//     for(int j=0;j<m-1;j++)
//       phi[j][n-1] = phi[j][0];

//     if(!config->save_last_only)
//       save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,
//                           &str_end_idx);
    
//   }

//   // Save last iteration if it wasn't saved
//   if(iter%config->qtimes != 0 || config->save_last_only){
//     iter--; // To get the correct iteration number
//     sprintf(buffer,"L");
//     save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,&str_end_idx);
//   }

//   free(L_phi);
//   free(Cij);
//   free(rho);
//   free(f);
//   free(ak);
//   free(bk);
//   free(ck);
//   free(uk);
//   free(fk);
//   free(an);
//   free(bn);
//   free(cn);
//   free(un);
//   free(fn);
//   free(filename_save);
//   free(buffer);
//   free(filename_log);
//   fclose(file_log);
//   free(fpb_prmtrs.L2);
//   free(fpb_prmtrs.Linf);
// }
*/

/*
- gigiaero, 15/07/2026, 1241 hours
*/
void save_prmtrs_fullp(char *casename,fullp_prmtrs *fp_prmtrs,
                       char *fname_msh_x,char *fname_msh_y){
  char *filename = malloc(sizeof(char)*200);
  FILE *output;

  sprintf(filename,"%s_prmtrs_fullp.dat",casename);

  output = fopen(filename,"w");

  fprintf(output,"fp.prmtrs.alpha = %f;\n",fp_prmtrs->alpha);
  fprintf(output,"fp.prmtrs.Ma = %f;\n",fp_prmtrs->Ma);
  fprintf(output,"fp.prmtrs.C = %f;\n",fp_prmtrs->C);
  fprintf(output,"fp.prmtrs.beta_sub = %f;\n",fp_prmtrs->beta_sub);
  fprintf(output,"fp.prmtrs.beta_super = %f;\n",fp_prmtrs->beta_super);
  fprintf(output,"fp.prmtrs.lift = %d;\n",fp_prmtrs->lift);
  fprintf(output,"fp.prmtrs.Rg = %f;\n\n",fp_prmtrs->Rg);

  fprintf(output,"char filename_msh_x[] = \"%s\";\n",fname_msh_x);
  fprintf(output,"char filename_msh_y[] = \"%s\";\n",fname_msh_y);

  fclose(output);
  free(filename);
}

/*
FINALLY done
- gigiaero, 12/07/2026, 1926 hours
*/
void solve_af2_2d_rectangular_fullp(int m,int n,double phi[m][n],double J[m][n],
                                    double A1[m][n],double A2[m][n],
                                    double A3[m][n],double J_ih[m][n],
                                    double A1_ih[m][n],double A2_ih[m][n],
                                    double A3_ih[m][n],double J_jh[m][n],
                                    double A1_jh[m][n],double A2_jh[m][n],
                                    double A3_jh[m][n],double *x,double *y,
                                    sim_prmtrs *config,fullp_prmtrs *fp_prmtrs){
  // Solver variables
  double (*L_phi)[n] = calloc(m,sizeof *L_phi);
  double (*Cij)[n] = calloc(m,sizeof *Cij);
  double (*f)[n] = calloc(m,sizeof *f);
  double (*rho)[n] = calloc(m,sizeof *rho);
  double (*rho_ih)[n] = calloc(m,sizeof *rho_ih);
  double (*rho_jh)[n] = calloc(m,sizeof *rho_jh);
  double (*contraU)[n] = calloc(m,sizeof *contraU);
  double (*contraU_ih)[n] = calloc(m,sizeof *contraU_ih);
  double (*contraU_jh)[n] = calloc(m,sizeof *contraU_jh);
  double (*contraV)[n] = calloc(m,sizeof *contraV);
  double (*contraV_ih)[n] = calloc(m,sizeof *contraV_ih);
  double (*contraV_jh)[n] = calloc(m,sizeof *contraV_jh);
  double (*rho_til)[n] = calloc(m,sizeof *rho_til);
  double (*rho_bar)[n] = calloc(m,sizeof *rho_bar);
  double (*nu_ih)[n] = calloc(m,sizeof *nu_ih);
  double (*nu_jh)[n] = calloc(m,sizeof *nu_jh);
  double (*L_phi_terms_ih)[n-1] = calloc(m,sizeof *L_phi_terms_ih);
  double (*L_phi_terms_jh)[n-1] = calloc(m,sizeof *L_phi_terms_jh);
  double (*Ai)[n] = calloc(m,sizeof *Ai);
  double (*Aj)[n] = calloc(m,sizeof *Aj);
  double *ax = malloc(sizeof(double)*(n-1));
  double *bx = malloc(sizeof(double)*(n-1));
  double *cx = malloc(sizeof(double)*(n-1));
  double *fx = malloc(sizeof(double)*(n-1));
  double *ux = malloc(sizeof(double)*(n-1));
  double *bn = malloc(sizeof(double)*(m-1));
  double *cn = malloc(sizeof(double)*(m-2));
  double *fn = malloc(sizeof(double)*(m-1));
  double *un = malloc(sizeof(double)*(m-1));
  double dphi_dksi,dphi_deta;
  double beta_sub = fp_prmtrs->beta_sub;
  double beta_super;
  double beta;
  fp_beta_prmtrs fpb_prmtrs;
  double (*Gamma)[3] = calloc(m,sizeof *Gamma);
  // double Gamma[] = {0.,0.,0.};
  // double Gamma[] = {phi[m-1][1] - phi[m-1][n-2],phi[m-1][1] - phi[m-1][n-2],
  //                   phi[m-1][1] - phi[m-1][n-2]};
  double omega = config->w;
  double alpha;
  int k = 1;
  int iter = 1;
  int idx;
  // Save files
  char *filename_save = malloc(sizeof(char)*200);
  char *buffer = malloc(sizeof(char)*200);
  int str_end_idx;
  // Residuals
  char *filename_log = malloc(sizeof(char)*200);
  double res,L2_res;
  FILE *file_log;

  // Configure log files
  sprintf(filename_log,"%s.log",config->casename);
  file_log = fopen(filename_log,"w");

  // Prepare string to save simulation data  
  sprintf(filename_save,"%s_iter_",config->casename);
  find_str_end(filename_save,&str_end_idx);
  
  beta_initial_config(config,&fpb_prmtrs,fp_prmtrs->beta_super);

  char filename[500];                                             // apagar depois

  for(iter;iter<=config->max_iter;iter++){
    alpha_sequence(&alpha,&k,iter,config);

    calc_contraUV(m,n,contraU_ih,phi,A1_ih,A2_ih,A3_ih,1,1);
    calc_contraUV(m,n,contraU_jh,phi,A1_jh,A2_jh,A3_jh,2,1);
    calc_contraUV(m,n,contraU,phi,A1,A2,A3,3,1);
    calc_contraUV(m,n,contraV_ih,phi,A1_ih,A2_ih,A3_ih,1,2);
    calc_contraUV(m,n,contraV_jh,phi,A1_jh,A2_jh,A3_jh,2,2);
    calc_contraUV(m,n,contraV,phi,A1,A2,A3,3,2);

    calc_rho(m,n,rho_ih,phi,A2_ih,A3_ih,contraU_ih,contraV_ih,1);
    calc_rho(m,n,rho_jh,phi,A2_jh,A3_jh,contraU_jh,contraV_jh,2);
    calc_rho(m,n,rho,phi,A2,A3,contraU,contraV,3);

    nu_switch(m,n,nu_ih,contraU_ih,rho,fp_prmtrs->C,1);
    nu_switch(m,n,nu_jh,contraV_jh,rho,fp_prmtrs->C,2);

    rho_coeffs(m,n,rho_til,rho_ih,contraU_ih,nu_ih,1,1);
    rho_coeffs(m,n,rho_bar,rho_jh,contraV_jh,nu_jh,2,2);

    L_phi_fullp_der_terms_ih(m,n,L_phi_terms_ih,rho_til,contraU_ih,J_ih);
    L_phi_fullp_der_terms_jh(m,n,L_phi_terms_jh,rho_bar,contraV_jh,J_jh);

    L_phi_fullp(m,n,L_phi,L_phi_terms_ih,L_phi_terms_jh);

    calc_Ai(m,n,Ai,rho_til,A1_ih,J_ih);
    calc_Aj(m,n,Aj,rho_bar,A3_jh,J_jh);

    beta_update(&beta_super,L2_res,res,config->M,&fpb_prmtrs,iter);

    res = 0.;
    for(int j=0;j<m;j++){
      for(int i=0;i<n-1;i++){
        if(fabs(L_phi[j][i]) > res)
          res = fabs(L_phi[j][i]);
      }
    }

    L2_res = norm_L2(m,n,L_phi);

    printf("AF2 Iteration %010d | Res %.6e\n",iter,res);

    fprintf(file_log,"%.6e\n",res);

    // Test for convergence
    if(res <= config->eps && iter != 0){
      puts("<< Convergence! >>");
      iter++;
      break;
    }

    if(res >= div_ref){
      puts("- Divergence");
      iter++;
      break;
    }

    // Solve for corrections - step 1 (eta)
    for(int i=0;i<n-1;i++){
      bn[0] = alpha + Aj[0][i];
      cn[0] = -Aj[1][i];

      fn[0] = alpha*omega*L_phi[1][i];

      for(int j=2;j<m-1;j++){
        bn[j-1] = alpha + Aj[j-1][i];
        cn[j-1] = -Aj[j][i];

        fn[j-1] = alpha*omega*L_phi[j][i];
      }

      bn[m-2] = alpha + Aj[m-2][i];

      fn[m-2] = alpha*omega*L_phi[m-1][i];

      bidiagonal_up_matrix_solver(m-1,bn,cn,fn,un);

      for(int j=1;j<m;j++)
        f[j][i] = un[j-1];
    }
    
    // Solve for corrections - step 2 (ksi)
    for(int j=1;j<m;j++){
      for(int i=0;i<n-1;i++){
        dphi_dksi = get_dphi_dksi(m,n,phi,i,j,3);
        dphi_deta = get_dphi_deta(m,n,phi,A2[j][i],A3[j][i],dphi_dksi,i,j,3);

        beta_switch(&beta,beta_sub,beta_super,
                    contraU[j][i]*dphi_dksi + contraV[j][i]*dphi_deta >= 1.);

        if(i == 0)
          idx = n-2;
          
        else
          idx = i-1;

        ax[i] = -Ai[j][idx];
        bx[i] = alpha + Ai[j][idx] + Ai[j][i] + alpha*beta;
        cx[i] = -Ai[j][i];

        fx[i] = f[j][i] + alpha*Cij[j-1][i]; 

        // if(contraU_ih[j][i] >= 0)
        if(contraU[j][i] >= 0)
          ax[i] -= alpha*beta;

        else
          cx[i] -= alpha*beta;

        if(fp_prmtrs->lift && i == 0)
          fx[i] -= Ai[j][idx]*(Gamma[j][0] - Gamma[j][1]);
        
        // if(fp_prmtrs->lift && i == n-1)
        //   fx[i] += Ai[j][i]*(Gamma[j][0] - Gamma[j][1]);
      }

      tridiagonal_pmatrix_solver(n-1,ax,bx,cx,fx,ux);

      for(int i=0;i<n-1;i++)
        Cij[j][i] = ux[i];
    }

    for(int j=0;j<m;j++){
      for(int i=0;i<n-1;i++)
        phi[j][i] += Cij[j][i];

      phi[j][n-1] = phi[j][0];
    }

    if(fp_prmtrs->lift){
      for(int j=0;j<m;j++){
        calc_circulation(Gamma[j],n,phi[j],fp_prmtrs->Rg);
        phi[j][0] = phi[j][n-1] + Gamma[j][0];
        // phi[j][n-1] = phi[j][0] - Gamma[0];
      }

      update_vortex(n,phi[0],x,y,Gamma[0][0],fp_prmtrs);
    }

    // disp(Gamma[0][0]);

    if(!config->save_last_only){
      flip_2d_array(m,n,phi);
      save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,
                          &str_end_idx);
      flip_2d_array(m,n,phi);
    }
  }

  // Save last iteration if it wasn't saved
  if(iter%config->qtimes != 0 || config->save_last_only){
    flip_2d_array(m,n,phi);
    iter--; // To get the correct iteration number
    sprintf(buffer,"L");
    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,&str_end_idx);
    flip_2d_array(m,n,phi);
  }

  free(L_phi);
  free(Cij);
  free(f);
  free(rho);
  free(rho_ih);
  free(rho_jh);
  free(contraU);
  free(contraU_ih);
  free(contraU_jh);
  free(contraV);
  free(contraV_ih);
  free(contraV_jh);
  free(rho_til);
  free(rho_bar);
  free(nu_ih);
  free(nu_jh);
  free(L_phi_terms_ih);
  free(L_phi_terms_jh);
  free(Ai);
  free(Aj);
  free(ax);
  free(bx);
  free(cx);
  free(ux);
  free(fx);
  free(bn);
  free(cn);
  free(un);
  free(fn);
  free(filename_save);
  free(buffer);
  free(filename_log);
  fclose(file_log);
  free(fpb_prmtrs.L2);
  free(fpb_prmtrs.Linf);
  free(Gamma);
}

/*
- gigiaero, 03/06/2026, 1544 hours
*/
void three_point_pol2_extrp(int m,int n,double A[m][n],double x[m][n],
                            double y[m][n],int i,int end){
  double d1,d2,d3;

  if(!end){
    d1 = sqrt(pow(x[1][i] - x[0][i],2.) + pow(y[1][i] - y[0][i],2.));
    d2 = sqrt(pow(x[2][i] - x[0][i],2.) + pow(y[2][i] - y[0][i],2.));
    d3 = sqrt(pow(x[3][i] - x[0][i],2.) + pow(y[3][i] - y[0][i],2.));
    
    A[0][i] = (d1*d1*d2*A[3][i] - d1*d1*d3*A[2][i] - d1*d2*d2*A[3][i] + 
               d1*d3*d3*A[2][i] + d2*d2*d3*A[1][i] - d2*d3*d3*A[1][i])/
              (d1*d1*d2 - d1*d1*d3 - d1*d2*d2  + d1*d3*d3  + d2*d2*d3 - 
               d2*d3*d3);
  }
  else{
    d1 = sqrt(pow(x[m-2][i] - x[m-1][i],2.) + pow(y[m-2][i] - y[m-1][i],2.));
    d2 = sqrt(pow(x[m-3][i] - x[m-1][i],2.) + pow(y[m-3][i] - y[m-1][i],2.));
    d3 = sqrt(pow(x[m-4][i] - x[m-1][i],2.) + pow(y[m-4][i] - y[m-1][i],2.));

    A[m-1][i] = (d1*d1*d2*A[m-4][i] - d1*d1*d3*A[m-3][i] - d1*d2*d2*A[m-4][i] + 
                 d1*d3*d3*A[m-3][i] + d2*d2*d3*A[m-2][i] - d2*d3*d3*A[m-2][i])/
                (d1*d1*d2 - d1*d1*d3 - d1*d2*d2  + d1*d3*d3  + d2*d2*d3 - 
                 d2*d3*d3);
  }
}

/*
- gigiaero, 03/06/2026, 1306 hours
*/
void two_point_linear_extrp(int m,int n,double A[m][n],double x[m][n],
                            double y[m][n],int i,int end){
  double d1,d2;
  
  if(!end){
    d1 = sqrt(pow(x[1][i] - x[0][i],2.) + pow(y[1][i] - y[0][i],2.));
    d2 = sqrt(pow(x[2][i] - x[0][i],2.) + pow(y[2][i] - y[0][i],2.));

    A[0][i] = (d1*A[2][i] - d2*A[1][i])/(d1 - d2);
  }
  else{
    d1 = sqrt(pow(x[m-2][i] - x[m-1][i],2.) + pow(y[m-2][i] - y[m-1][i],2.));
    d2 = sqrt(pow(x[m-3][i] - x[m-1][i],2.) + pow(y[m-3][i] - y[m-1][i],2.));
    // d1 = sqrt(x[m-2][i]*x[m-2][i] + y[m-2][i]*y[m-2][i]);
    // d2 = sqrt(x[m-3][i]*x[m-3][i] + y[m-3][i]*y[m-3][i]);

    A[m-1][i] = (d1*A[m-3][i] - d2*A[m-2][i])/(d1 - d2);
  }
}

/*
- gigiaero, 12/07/2026, 2337 hours
*/
void update_vortex(int n,double *phi,double *x,double *y,double Gamma,
                   fullp_prmtrs *fp_prmtrs){
  double Uinf = freestream_u(fp_prmtrs->Ma);
  double beta = sqrt(1. - fp_prmtrs->Ma*fp_prmtrs->Ma);

  for(int i=0;i<n;i++){
    phi[i] = Uinf*cos(fp_prmtrs->alpha)*x[i] + Uinf*sin(fp_prmtrs->alpha)*y[i]
             - Gamma/2./pi*atan2(beta*y[i],x[i]);
  }
}