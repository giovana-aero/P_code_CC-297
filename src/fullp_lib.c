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
- gigiaero, 29/05/2026, 1446 hours
*/
void calc_A_metrics(int m,int n,double A1[m][n],double A2[m][n],double A3[m][n],
                    double x[m][n],double y[m][n],double J[m][n]){
  double ksix,ksiy;
  double etax,etay;

  for(int j=1;j<m-1;j++){
    for(int i=0;i<n-1;i++){
      metric_terms(&ksix,&ksiy,&etax,&etay,m,n,J,x,y,i,j);

      A1[j][i] = ksix*ksix + ksiy*ksiy;
      A2[j][i] = ksix*etax + ksiy*etay;
      A3[j][i] = etax*etax + etay*etay;
    }
    
    A1[j][n-1] = A1[j][0];
    A2[j][n-1] = A2[j][0];
    A3[j][n-1] = A3[j][0];
  }

  for(int i=0;i<n-1;i++){
    // print_2d_array(m,n,A3);
    three_point_pol2_extrp(m,n,A1,x,y,i,0);
    three_point_pol2_extrp(m,n,A2,x,y,i,0);
    three_point_pol2_extrp(m,n,A3,x,y,i,0);
    three_point_pol2_extrp(m,n,A1,x,y,i,1);
    three_point_pol2_extrp(m,n,A2,x,y,i,1);
    three_point_pol2_extrp(m,n,A3,x,y,i,1);
    // two_point_linear_extrp(m,n,A1,x,y,i,0);
    // two_point_linear_extrp(m,n,A2,x,y,i,0);
    // two_point_linear_extrp(m,n,A3,x,y,i,0);
    // two_point_linear_extrp(m,n,A1,x,y,i,1);
    // two_point_linear_extrp(m,n,A2,x,y,i,1);
    // two_point_linear_extrp(m,n,A3,x,y,i,1);
  }

  A1[0][n-1] = A1[0][0];
  A2[0][n-1] = A2[0][0];
  A3[0][n-1] = A3[0][0];
  A1[m-1][n-1] = A1[m-1][0];
  A2[m-1][n-1] = A2[m-1][0];
  A3[m-1][n-1] = A3[m-1][0];
}

/*
- gigiaero, 03/06/2026, 2105 hours
NOTE: adapted to af2 only for now
*/
double calc_Ai(int m,int n,double phi[m][n],double J[m][n],double A1[m][n],
               double A2[m][n],double A3[m][n],double rho[m][n-1],double C,
               int i,int j,int af2_flip){
  i--;

  if(i < 0)
    i = n-2;

  double dphi_dksi,dphi_deta;
  double contraU,rho_til,nu_x;
  double A1_ih = (A1[j][i+1] + A1[j][i])*.5;
  double A2_ih = (A2[j][i+1] + A2[j][i])*.5;
  
  dphi_dksi = uniform_scheme_der1_o1_forward(m,n,phi,i,j,1);

  if(j == 0 && !af2_flip || j == m-1 && af2_flip)
    dphi_deta = -A2_ih*dphi_dksi/(A3[j][i+1] + A3[j][i])*2.;

  else
    dphi_deta = (uniform_scheme_der1_o2_central(m,n,phi,i,j,2) + 
                 uniform_scheme_der1_o2_central(m,n,phi,i+1,j,2))*.5;

  contraU = calc_contraU(dphi_dksi,dphi_deta,A1_ih,A2_ih);

  nu_x = nu_switch(m,n,rho,C,contraU,i,j,1);

  rho_til = rho_coeffs(m,n,rho,nu_x,rs_idx(contraU),i,j,1);

  return rho_til*A1_ih/((J[j][i+1] + J[j][i])*.5);
}

/*
- gigiaero, 18/06/2026, 2109 hours
NOTE: adapted to af2 only for now
*/
double calc_Aj(int m,int n,double phi[m][n],double J[m][n],double A2[m][n],
               double A3[m][n],double rho[m][n-1],double C,int i,int j){
  j--;

  double dphi_dksi,dphi_deta;
  double contraV,rho_bar,nu_y;
  double A3_ih = (A3[j+1][i] + A3[j][i])*.5;

  if(i == 0)
    dphi_dksi = (uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j) + 
                 uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j+1))*.5;
  else
    dphi_dksi = (uniform_scheme_der1_o2_central(m,n,phi,i,j,1) + 
                 uniform_scheme_der1_o2_central(m,n,phi,i,j+1,1))*.5;

  dphi_deta = uniform_scheme_der1_o1_forward(m,n,phi,i,j,2);

  contraV = calc_contraU(dphi_dksi,dphi_deta,(A2[j+1][i] + A2[j][i])*.5,
                         A3_ih);

  nu_y = nu_switch(m,n,rho,C,contraV,i,j,2);

  rho_bar = rho_coeffs(m,n,rho,nu_y,rs_idx(contraV),i,j,2);

  return rho_bar*A3_ih/((J[j+1][i] + J[j][i])*.5);
}

// /*
// - gigiaero, 01/06/2026,1332 hours
// */
// // NOTA: PENDENTE ADAPTAÇÃO PRO AF2
// double L_phi_fullp_der_terms_jh(int m,int n,double phi[m][n],double J[m][n],
//                                 double A2[m][n],double A3[m][n],
//                                 double rho[m][n-1],double C,int i,int j){
//   double dphi_dksi,dphi_deta;
//   double contraV,rho_bar,nu_y;
//   double A3_jh = (A3[j+1][i] + A3[j][i])*.5;

//   if(i == 0)
//     dphi_dksi = (uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j) + 
//                  uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j+1))*.5;
//   else
//     dphi_dksi = (uniform_scheme_der1_o2_central(m,n,phi,i,j,1) + 
//                  uniform_scheme_der1_o2_central(m,n,phi,i,j+1,1))*.5;
  
//   dphi_deta = uniform_scheme_der1_o1_forward(m,n,phi,i,j,2);

//   contraV = calc_contraV(dphi_dksi,dphi_deta,(A2[j+1][i] + A2[j][i])*.5,
//                          (A3[j+1][i] + A3[j][i])*.5);

//   nu_y = nu_switch(m,n,rho,C,contraV,i,j,2);

//   rho_bar = rho_coeffs(m,n,rho,nu_y,rs_idx(contraV),i,j,2);

//   return rho_bar*contraV/((J[j+1][i] + J[j][i])*.5);
// }

/*
- gigiaero, 31/05/2026, 1143h hours
*/
double calc_contraU(double dphi_dksi,double dphi_deta,double A1_val,
                    double A2_val){
  return A1_val*dphi_dksi + A2_val*dphi_deta;
}

/*
- gigiaero, 31/05/2026, 1321 hours
*/
double calc_contraV(double dphi_dksi,double dphi_deta,double A2_val,
                    double A3_val){
  return A2_val*dphi_dksi + A3_val*dphi_deta;
}

/*
- gigiaero, 29/05/2026, 1435 hours
*/
// a se considerar: calcular nas bordas usando interpolação linear
void calc_J(int m,int n,double J[m][n],double x[m][n],double y[m][n]){
  for(int j=1;j<m-1;j++){
    J[j][0] = 1./(uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j)*
                  uniform_scheme_der1_o2_central(m,n,y,0,j,2) - 
                  uniform_scheme_der1_o2_central(m,n,x,0,j,2)*
                  uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j));

    for(int i=1;i<n-1;i++)
      J[j][i] = 1./(uniform_scheme_der1_o2_central(m,n,x,i,j,1)*
                    uniform_scheme_der1_o2_central(m,n,y,i,j,2) - 
                    uniform_scheme_der1_o2_central(m,n,x,i,j,2)*
                    uniform_scheme_der1_o2_central(m,n,y,i,j,1));
  }

  // for(int i=0;i<n-1;i++){
  //   // two_point_linear_extrp(m,n,J,x,y,i,0);
  //   // two_point_linear_extrp(m,n,J,x,y,i,1);
  //   three_point_pol2_extrp(m,n,J,x,y,i,0);
  //   three_point_pol2_extrp(m,n,J,x,y,i,1);
  // }

  for(int i=1;i<n-1;i++){
    J[0][i] = 1./(uniform_scheme_der1_o2_central(m,n,x,i,0,1)*
                  uniform_scheme_der1_o2_forward(m,n,y,i,0,2) - 
                  uniform_scheme_der1_o2_forward(m,n,x,i,0,2)*
                  uniform_scheme_der1_o2_central(m,n,y,i,0,1));
    J[m-1][i] = 1./(uniform_scheme_der1_o2_central(m,n,x,i,m-1,1)*
                    uniform_scheme_der1_o2_backward(m,n,y,i,m-1,2) - 
                    uniform_scheme_der1_o2_backward(m,n,x,i,m-1,2)*
                    uniform_scheme_der1_o2_central(m,n,y,i,m-1,1));
  }

  J[0][0] = 1./(uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,0)*
                uniform_scheme_der1_o2_forward(m,n,y,0,0,2) - 
                uniform_scheme_der1_o2_forward(m,n,x,0,0,2)*
                uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,0));
  J[m-1][0] = 1./(uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,m-1)*
                  uniform_scheme_der1_o2_backward(m,n,y,0,m-1,2) - 
                  uniform_scheme_der1_o2_backward(m,n,x,0,m-1,2)*
                  uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,m-1));

  for(int j=0;j<m;j++)
    J[j][n-1] = J[j][0];
}

/*
- gigiaero, 01/06/2026, 0956 hours
*/
void calc_rho(int m,int n,double rho[m][n-1],double phi[m][n],double A1[m][n],
              double A2[m][n],double A3[m][n],int af2_flip){
  double dphi_dksi,dphi_deta;
  double contraU,contraV;
  double A1_mean,A2_mean,A3_mean;

  for(int j=0+af2_flip;j<m-1+af2_flip;j++){
    for(int i=0;i<n-1;i++){
      dphi_dksi = uniform_scheme_der1_o1_forward(m,n,phi,i,j,1);

      A1_mean = (A1[j][i+1] + A1[j][i])*.5;
      A2_mean = (A2[j][i+1] + A2[j][i])*.5;
      A3_mean = (A3[j][i+1] + A3[j][i])*.5;

      // if(j == 0 && !af2_flip || j == m-1 && af2_flip)
      if(j == 0 || j == m-1)
        dphi_deta = -A2_mean*dphi_dksi/A3_mean;

      else
        dphi_deta = (uniform_scheme_der1_o2_central(m,n,phi,i,j,2) + 
                     uniform_scheme_der1_o2_central(m,n,phi,i+1,j,2))*.5;

      contraU = calc_contraU(dphi_dksi,dphi_deta,A1_mean,A2_mean);
      contraV = calc_contraV(dphi_dksi,dphi_deta,A2_mean,A3_mean);

      rho[j][i] = pow(1 - ((gamma - 1.)/(gamma + 1.))*(contraU*dphi_dksi + 
                      contraV*dphi_deta),1./(gamma - 1.));
      // disp(A1_mean);
      // disp(A2_mean);
      // disp(A3_mean);
      // disp(contraU);
      // disp(contraV);
      // disp(dphi_dksi);
      // disp(dphi_deta);
      // // disp(contraU*dphi_dksi);
      // // disp(contraV*dphi_deta);
      // disp(contraU*dphi_dksi + contraV*dphi_deta);
      // disp(1 - ((gamma - 1.)/(gamma + 1.))*(contraU*dphi_dksi + contraV*dphi_deta));
      // disp(1./(gamma - 1.));
      // disp(rho[j][i]);
      // putchar('\n');
    }
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
  double (*phi)[n] = calloc(m,sizeof *phi);
  double (*J)[n] = calloc(m,sizeof *J);
  double (*A1)[n] = calloc(m,sizeof *A1);
  double (*A2)[n] = calloc(m,sizeof *A2);
  double (*A3)[n] = calloc(m,sizeof *A3);
  int af2_flip = 0;

  save_prmtrs_sim(config);
  // salvar aqui parâmetros específicos do fullp ////////

  read_2d_array_from_file(m,n,x,fname_msh_x);
  read_2d_array_from_file(m,n,y,fname_msh_y);

  char *filename = malloc(sizeof(char)*200);
  sprintf(filename,"%s_x_mesh.dat",config->casename);
  print_2d_array_to_file(m,n,x,filename,0);
  sprintf(filename,"%s_y_mesh.dat",config->casename);
  print_2d_array_to_file(m,n,y,filename,0);

  if(config->Ntype == 3){
    flip_2d_array(m,n,x);
    flip_2d_array(m,n,y);
    af2_flip = 1;
  }

  print_2d_array_to_file(m,n,x,"mesh_x.dat",0);
  print_2d_array_to_file(m,n,y,"mesh_y.dat",0);

  calc_J(m,n,J,x,y);
  calc_A_metrics(m,n,A1,A2,A3,x,y,J);

  print_2d_array_to_file(m,n,J,"mat_J.dat",0);
  print_2d_array_to_file(m,n,A1,"mat_A1.dat",0);
  print_2d_array_to_file(m,n,A2,"mat_A2.dat",0);
  print_2d_array_to_file(m,n,A3,"mat_A3.dat",0);

  initialize_fullp(m,n,phi,x,y,fp_prmtrs);

  alpha_sequence_aH(m,n,x,y,config,1);

  if(config->save_i_c){
    sprintf(filename,"%s_iter_%010d.dat",config->casename,0);
    print_2d_array_to_file(m,n,phi,filename,0);
  }

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
      // puts("fullp af2 pending");
      puts("AF2, full potential flow over airfoil");
      // solve_af2_2d_rectangular_fullp(m,n,phi,x,y,config);
      flip_2d_array(m,n,phi);
      flip_2d_array(m,n,J);
      break;
      
    default:
      puts("evaluate_delta_form_eom: Invalid Ntype");
      exit(32);
  }

  double (*u)[n] = calloc(m,sizeof *u);
  double (*v)[n] = calloc(m,sizeof *v);
  double (*Ve)[n] = calloc(m,sizeof *Ve);

  // nota: pensar se devo pegar a densidade nas fronteiras

  get_u_v_potential_fullp(m,n,phi,x,y,J,u,v,Ve);

  sprintf(filename,"%s_u.dat",config->casename);
  print_2d_array_to_file(m,n,u,filename,0);
  sprintf(filename,"%s_v.dat",config->casename);
  print_2d_array_to_file(m,n,v,filename,0);
  sprintf(filename,"%s_Ve.dat",config->casename);
  print_2d_array_to_file(m,n,Ve,filename,0);

  

  free(filename);
  free(x);
  free(y);
  free(phi);
  free(J);
  free(A1);
  free(A2);
  free(A3);
  free(u);
  free(v);
  free(Ve);
}

/*
- gigiaero, 01/06/2026, 1555 hours
*/
double freestream_u(double Ma){
  return sqrt((gamma + 1.)/(gamma - 1. + 2./Ma/Ma));
}

/*
- gigiaero, 01/06/2026, 2204 hours
*/
void get_u_v_potential_fullp(int m,int n,double phi[m][n],double x[m][n],
                             double y[m][n],double J[m][n],double u[m][n],
                             double v[m][n],double Ve[m][n]){
  double ksix,ksiy;
  double etax,etay;

  // Domain interior
  for(int j=1;j<m-1;j++){
    metric_terms(&ksix,&ksiy,&etax,&etay,m,n,J,x,y,0,j);
    u[j][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j)*ksix + 
              uniform_scheme_der1_o2_central(m,n,phi,0,j,2)*etax;
    v[j][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j)*ksiy + 
              uniform_scheme_der1_o2_central(m,n,phi,0,j,2)*etay;
    
    for(int i=1;i<n-1;i++){
      metric_terms(&ksix,&ksiy,&etax,&etay,m,n,J,x,y,i,j);
      u[j][i] = uniform_scheme_der1_o2_central(m,n,phi,i,j,1)*ksix + 
                uniform_scheme_der1_o2_central(m,n,phi,i,j,2)*etax;
      v[j][i] = uniform_scheme_der1_o2_central(m,n,phi,i,j,1)*ksiy + 
                uniform_scheme_der1_o2_central(m,n,phi,i,j,2)*etay;
    }
  }

  // External boundary and airfoil surface
  for(int i=1;i<n-1;i++){
    metric_terms(&ksix,&ksiy,&etax,&etay,m,n,J,x,y,i,0);
    u[0][i] = uniform_scheme_der1_o2_central(m,n,phi,i,0,1)*ksix + 
              uniform_scheme_der1_o2_forward(m,n,phi,i,0,2)*etax;
    v[0][i] = uniform_scheme_der1_o2_central(m,n,phi,i,0,1)*ksiy + 
              uniform_scheme_der1_o2_forward(m,n,phi,i,0,2)*etay;

    metric_terms(&ksix,&ksiy,&etax,&etay,m,n,J,x,y,i,m-1);
    u[m-1][i] = uniform_scheme_der1_o2_central(m,n,phi,i,m-1,1)*ksix + 
                uniform_scheme_der1_o2_backward(m,n,phi,i,m-1,2)*etax;
    v[m-1][i] = uniform_scheme_der1_o2_central(m,n,phi,i,m-1,1)*ksiy + 
                uniform_scheme_der1_o2_backward(m,n,phi,i,m-1,2)*etay;
  }

  // Corners
  metric_terms(&ksix,&ksiy,&etax,&etay,m,n,J,x,y,0,0);
  u[0][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,0)*ksix + 
            uniform_scheme_der1_o2_forward(m,n,phi,0,0,2)*etax;
  v[0][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,0)*ksiy + 
            uniform_scheme_der1_o2_forward(m,n,phi,0,0,2)*etay;

  metric_terms(&ksix,&ksiy,&etax,&etay,m,n,J,x,y,0,m-1);
  u[m-1][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,m-1)*ksix + 
              uniform_scheme_der1_o2_backward(m,n,phi,0,m-1,2)*etax;
  v[m-1][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,m-1)*ksiy + 
              uniform_scheme_der1_o2_backward(m,n,phi,0,m-1,2)*etay;

  // Velocity resultant
  for(int j=0;j<m;j++){
    for(int i=0;i<n-1;i++)
      Ve[j][i] = sqrt(pow(u[j][i],2) + pow(v[j][i],2));

    u[j][n-1] = u[j][0];
    v[j][n-1] = v[j][0];
    Ve[j][n-1] = Ve[j][0];
  }
}

/*
- gigiaero, 0/06/2026, 2057 hours
*/
void initialize_fullp(int m,int n,double phi[m][n],double x[m][n],double y[m][n],
                      fullp_prmtrs *fp_prmtrs){
  double Uinf = freestream_u(fp_prmtrs->Ma);

  fp_prmtrs->alpha *= pi/180.;

  // puts("Uinf = "); disp(Uinf);

  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      phi[j][i] = Uinf*cos(fp_prmtrs->alpha)*x[j][i] + 
                  Uinf*sin(fp_prmtrs->alpha)*y[j][i];
  }
}

/*
- gigiaero, 01/06/2026,1314 hours
*/
double L_phi_fullp_der_terms_ih(int m,int n,double phi[m][n],double J[m][n],
                                double A1[m][n],double A2[m][n],double A3[m][n],
                                double rho[m][n-1],double C,int i,int j,
                                int af2_flip){
  double dphi_dksi,dphi_deta;
  double contraU,rho_til,nu_x;
  double A2_ih = (A2[j][i+1] + A2[j][i])*.5;

  // if(i != 0 && j == 0){
  //   dphi_dksi = uniform_scheme_der1_o1_forward(m,n,phi,0,j,1);
  //   dphi_deta = 

  //   contraU = calc_contraU(dphi_dksi,dphi_deta,(A1[j][1] + A1[j][0])*.5,
  //                          (A2[j][1] + A2[j][0])*.5);

  //   nu_x = nu_switch(m,n,rho,C,contraU,0,j,1);

  //   rho_til = rho_coeffs(m,n,rho,nu_x,rs_idx(contraU),0,j,1);

  //   return rho_til*contraU/((J[j][1] + J[j][0])*.5);
  // }

  // if(i == 0){
  //   dphi_dksi = uniform_scheme_der1_o1_forward(m,n,phi,i,j,1);
  //   dphi_deta = (uniform_scheme_der1_o2_central(m,n,phi,i,j,2) + 
  //                uniform_scheme_der1_o2_central(m,n,phi,i+1,j,2))*.5;

  //   contraU = calc_contraU(dphi_dksi,dphi_deta,(A1[j][i+1] + A1[j][i])*.5,
  //                          (A2[j][i+1] + A2[j][i])*.5);

  //   nu_x = nu_switch(m,n,rho,C,contraU,i,j,1);

  //   rho_til = rho_coeffs(m,n,rho,nu_x,rs_idx(contraU),i,j,1);

  //   return rho_til*contraU/((J[j][i+1] + J[j][i])*.5);
  // }
  // else{
  
  dphi_dksi = uniform_scheme_der1_o1_forward(m,n,phi,i,j,1);

  if(j == 0 && !af2_flip || j == m-1 && af2_flip)
    dphi_deta = -A2_ih*dphi_dksi/(A3[j][i+1] + A3[j][i])*2.;

  else
    dphi_deta = (uniform_scheme_der1_o2_central(m,n,phi,i,j,2) + 
                 uniform_scheme_der1_o2_central(m,n,phi,i+1,j,2))*.5;

  contraU = calc_contraU(dphi_dksi,dphi_deta,(A1[j][i+1] + A1[j][i])*.5,
                         A2_ih);

  nu_x = nu_switch(m,n,rho,C,contraU,i,j,1);

  rho_til = rho_coeffs(m,n,rho,nu_x,rs_idx(contraU),i,j,1);

  return rho_til*contraU/((J[j][i+1] + J[j][i])*.5);
  // }
}

/*
- gigiaero, 01/06/2026,1332 hours
*/
// NOTA: PENDENTE ADAPTAÇÃO PRO AF2
double L_phi_fullp_der_terms_jh(int m,int n,double phi[m][n],double J[m][n],
                                double A2[m][n],double A3[m][n],
                                double rho[m][n-1],double C,int i,int j){
  double dphi_dksi,dphi_deta;
  double contraV,rho_bar,nu_y;

  if(i == 0)
    dphi_dksi = (uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j) + 
                 uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j+1))*.5;
  else
    dphi_dksi = (uniform_scheme_der1_o2_central(m,n,phi,i,j,1) + 
                 uniform_scheme_der1_o2_central(m,n,phi,i,j+1,1))*.5;
  
  dphi_deta = uniform_scheme_der1_o1_forward(m,n,phi,i,j,2);

  contraV = calc_contraV(dphi_dksi,dphi_deta,(A2[j+1][i] + A2[j][i])*.5,
                         (A3[j+1][i] + A3[j][i])*.5);

  nu_y = nu_switch(m,n,rho,C,contraV,i,j,2);

  rho_bar = rho_coeffs(m,n,rho,nu_y,rs_idx(contraV),i,j,2);

  return rho_bar*contraV/((J[j+1][i] + J[j][i])*.5);
}

/*
- gigiaero, 01/06/2026,1315 hours
*/
void L_phi_fullp(int m,int n,double L_phi[m][n],double phi[m][n],
                 double J[m][n],double A1[m][n],double A2[m][n],double A3[m][n],
                 double rho[m][n-1],double C,int af2_flip){
  double term_ih,term_ihm1;
  double term_jh,term_jhm1;

  for(int j=0+af2_flip;j<m-1+af2_flip;j++){
    term_ih = L_phi_fullp_der_terms_ih(m,n,phi,J,A1,A2,A3,rho,C,0,j,af2_flip);
    term_ihm1=L_phi_fullp_der_terms_ih(m,n,phi,J,A1,A2,A3,rho,C,n-2,j,af2_flip);

    term_jh = L_phi_fullp_der_terms_jh(m,n,phi,J,A2,A3,rho,C,0,j);
    term_jhm1 = L_phi_fullp_der_terms_jh(m,n,phi,J,A2,A3,rho,C,0,j-1);

    L_phi[j][0] = (term_ih - term_ihm1) + (term_jh - term_jhm1);

    for(int i=1;i<n-1;i++){
      term_ihm1 = term_ih;
      term_ih = L_phi_fullp_der_terms_ih(m,n,phi,J,A1,A2,A3,rho,C,i,j,af2_flip);
      // term_ihm1 = L_phi_fullp_der_terms_ih(m,n,phi,J,A1,A2,rho,C,i-1,j);

      term_jh = L_phi_fullp_der_terms_jh(m,n,phi,J,A2,A3,rho,C,i,j);
      term_jhm1 = L_phi_fullp_der_terms_jh(m,n,phi,J,A2,A3,rho,C,i,j-1);

      L_phi[j][i] = (term_ih - term_ihm1) + (term_jh - term_jhm1);
    }
  }

  // superfície do aerofólio, i=0 aqui -----------------------

  // for(int i=1;i<n-1;i++){
  //   term_ih = L_phi_fullp_der_terms_ih(m,n,phi,J,A1,A2,A3,rho,C,i,0);
  //   term_ihm1 = L_phi_fullp_der_terms_ih(m,n,phi,J,A1,A2,A3,rho,C,i-1,0);

  //   term_jh = L_phi_fullp_der_terms_jh(m,n,phi,J,A2,A3,rho,C,i,0);
  //   term_jhm1  = -term_jh;

  //   L_phi[0][i] = (term_ih - term_ihm1) + (term_jh - term_jhm1);
  // }

  for(int j=0+af2_flip;j<m-1+af2_flip;j++)
    L_phi[j][n-1] = L_phi[j][0];
}

/*
- gigiaero, 31/05/2026, 1326 hours
*/
double max2(double n1,double n2){
  if(n1 >= n2)
    return n1;
  
  else
    return n2;
}

// /*
// - gigiaero, 31/05/2026, 1123 hours
// */
// double mean2(double n1,double n2){
//   return (n1 + n2)/2.;
// }

/*
- gigiaero, 31/05/2026, 1434 hours
*/
double mean4_j(int m,int n,double rho[m][n-1],int i,int j){
  if(i == 0)
    return (rho[j][i] + rho[j][n-2] + rho[j+1][n-2] + rho[j+1][i])*.25;
  else
    return (rho[j][i] + rho[j][i-1] + rho[j+1][i-1] + rho[j+1][i])*.25;
}

/*
- gigiaero, 02/06/2026, 1024 hours
*/
void metric_terms(double *ksix,double *ksiy,double *etax,double *etay,int m,
                  int n,double J[m][n],double x[m][n],double y[m][n],
                  int i,int j){
  double dy_deta,dx_deta;
  double dy_dksi,dx_dksi;

  if(i == 0){
    dy_dksi = uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j);
    dx_dksi = uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j);
  }
  else{
    dy_dksi = uniform_scheme_der1_o2_central(m,n,y,i,j,1);
    dx_dksi = uniform_scheme_der1_o2_central(m,n,x,i,j,1);
  }

  if(j == 0){
    dy_deta = uniform_scheme_der1_o2_forward(m,n,y,i,j,2);
    dx_deta = uniform_scheme_der1_o2_forward(m,n,x,i,j,2);
  }
  else if(j == m-1){
    dy_deta = uniform_scheme_der1_o2_backward(m,n,y,i,j,2);
    dx_deta = uniform_scheme_der1_o2_backward(m,n,x,i,j,2);
  }
  else{
    dy_deta = uniform_scheme_der1_o2_central(m,n,y,i,j,2);
    dx_deta = uniform_scheme_der1_o2_central(m,n,x,i,j,2);
  }

  *ksix = J[j][i]*dy_deta;
  *ksiy = -J[j][i]*dx_deta;
  *etax = -J[j][i]*dy_dksi;
  *etay = J[j][i]*dx_dksi;
}

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

// /*
// - gigiaero, 04/06/2026, 1605 hours
// */
// double norm_Linf(int m,int n,double A[m][n]){
//   double k = 0.;
  
//   for(int j=0;j<m;j++){
//     for(int i=0;i<n-1;i++){
//       if(fabs(A[j][i]) > k)
//         k = fabs(A[j][i]);
//     }
//   }

//   return k;
// }

/*
- gigiaero, 01/06/2026, 1048 hours
*/
// NOTA: PENDENTE ADAPTAÇÃO PRO AF2
double nu_switch(int m,int n,double rho[m][n-1],double C,double contraUV,
                 int i,int j,int axis){
  double C1 = .633939;
  double C2 = 4.9325;
  double rho_integer_ij;
  double nu;

  if(contraUV >= 0){
    if(i == 0)
      rho_integer_ij = (rho[j][i] + rho[j][n-2])*.5;
    else
      rho_integer_ij = (rho[j][i] + rho[j][i-1])*.5;
  }
  else{
    switch(axis){
      case 1:
        if(i == n-2)
          rho_integer_ij = (rho[j][0] + rho[j][i])*.5;
        else
          rho_integer_ij = (rho[j][i+1] + rho[j][i])*.5;

        break;

      case 2:
        if(i == 0)
          rho_integer_ij = (rho[j+1][i] + rho[j+1][n-2])*.5;
        else
          rho_integer_ij = (rho[j+1][i] + rho[j+1][i-1])*.5;

        break;

      default:
        puts("nu_switch: invalid axis");
        exit(15);
    }
  }

  nu = max2(0.,(C1 - rho_integer_ij)*C2*C);

  if(nu > 1.)
    return 1.;
  
  else
    return nu;
}

/*
- gigiaero, 31/05/2026, 1444 hours
*/
double rho_coeffs(int m,int n,double rho[m][n-1],double nu,int rs,int i,int j,
                  int axis){
  switch(axis){
    case 1: // rho_til
      if(i+rs < 0)
        return (1. - nu)*rho[j][i] + nu*rho[j][n-2];
      else if(i+rs > n-2)
        return (1. - nu)*rho[j][i] + nu*rho[j][0];
      else
        return (1. - nu)*rho[j][i] + nu*rho[j][i+rs];

      break;

    case 2: // rho_bar
      return (1. - nu)*mean4_j(m,n,rho,i,j) + nu*mean4_j(m,n,rho,i,j+rs);
      break;

    default:
      puts("rho_coeffs: invalid axis");
      exit(15);
  }
}

/*
- gigiaero, 31/05/2026, 1406 hours
*/
int rs_idx(double contraUV){
  if(contraUV <= 0)
    return +1;
  else
    return -1;
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

void solve_af2_2d_rectangular_fullp(int m,int n,double phi[m][n],double J[m][n],
                                    double A1[m][n],double A2[m][n],
                                    double A3[m][n],sim_prmtrs *config,
                                    fullp_prmtrs *fp_prmtrs){
  // Solver variables
  double (*L_phi)[n] = calloc(m,sizeof *L_phi);
  double (*Cij)[n] = calloc(m,sizeof *Cij);
  double (*rho)[n-1] = calloc(m,sizeof *rho);
  double (*f)[n] = calloc(m,sizeof *f);
  double *ak = malloc(sizeof(double)*(n-1));
  double *bk = malloc(sizeof(double)*(n-1));
  double *ck = malloc(sizeof(double)*(n-1));
  double *fk = malloc(sizeof(double)*(n-1));
  double *uk = malloc(sizeof(double)*(n-1));
  double *an = malloc(sizeof(double)*(m-2));
  double *bn = malloc(sizeof(double)*(m-1));
  double *cn = malloc(sizeof(double)*(m-2));
  double *fn = malloc(sizeof(double)*(m-1));
  double *un = malloc(sizeof(double)*(m-1));
  double Ai,Aip1,Aj,Ajp1;
  double dphi_dksi,dphi_deta;
  double contraU,contraV;
  double beta_sub = fp_prmtrs->beta_sub;
  double beta_super;
  double beta;
  fp_beta_prmtrs fpb_prmtrs;
  // double beta,beta_n;
  double omega = config->w;
  double alpha;
  int k = 1;
  int iter = 1;
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

  // a se pensar: calc_contraU e calc_contraV poderiam ser a mesma função
  
  beta_initial_config(config,&fpb_prmtrs,fp_prmtrs->beta_super);

  for(iter;iter<=config->max_iter;iter++){

    calc_rho(m,n,rho,phi,A1,A2,A3,1);

    L_phi_fullp(m,n,L_phi,phi,J,A1,A2,A3,rho,fp_prmtrs->C,1);

    res = 0.;
    for(int j=0;j<m-1;j++){
      for(int i=0;i<n-1;i++){
        if(fabs(L_phi[j][i]) > res)
          res = fabs(L_phi[j][i]);
      }
    }

    L2_res = norm_L2(m,n,L_phi);

    printf("ADI Iteration %010d | Res %.6e\n",iter,res);

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
      Aj = calc_Aj(m,n,phi,J,A2,A3,rho,fp_prmtrs->C,i,j);
      Ajp1 = calc_Aj(m,n,phi,J,A2,A3,rho,fp_prmtrs->C,i,j+1);

      bn[0] = alpha + Aj;
      cn[0] = -Ajp1;

      fn[0] = alpha*omega*L_phi[1][i];

      for(int j=2;j<m-2;j++){
        Aj = calc_Aj(m,n,phi,J,A2,A3,rho,fp_prmtrs->C,i,j);
        Ajp1 = calc_Aj(m,n,phi,J,A2,A3,rho,fp_prmtrs->C,i,j+1);

        an[j-2] = 0.;
        bn[j-1] = alpha + Aj;
        cn[j-1] = -Ajp1;

        fn[j-1] = alpha*omega*L_phi[j][i];
      }

      Aj = calc_Aj(m,n,phi,J,A2,A3,rho,fp_prmtrs->C,i,j);
      Ajp1 = calc_Aj(m,n,phi,J,A2,A3,rho,fp_prmtrs->C,i,j+1);

      an[m-3] = alpha + Aj;
      bn[m-2] = -Ajp1;

      fn[m-2] = alpha*omega*L_phi[m-1][i];

      tridiagonal_matrix_solver(m-1,an,bn,cn,fn,un);  // NOTA: escrever uma função pra sistema bidiagonal superior e substituir aqui

      for(int j=1;j<m;j++)
        f[j][i] = un[j-1];      
    }
    
    
    // Solve for corrections - step 2 (ksi)


    // LEMBRETE: o phi deve ser espelhado antes de ser salvo em arquivo
  }

  free(L_phi);
  free(Cij);
  free(rho);
  free(f);
  free(ak);
  free(bk);
  free(ck);
  free(uk);
  free(fk);
  free(an);
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
}

double sound_speed(int m,int n,double phi[m][n],int i,int j){

}