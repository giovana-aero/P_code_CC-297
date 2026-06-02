#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"../include/2d_arrays.h"
#include"../include/eom_lib.h"
#include"../include/fullp_lib.h"
#include"../include/num_methods.h"

#define div_ref 1e100
#define gamma 1.4
#define pi 3.14159265358979323846261003

void boundary_condition_freestream(){
  // phi_inf = u_inf*x (simples assim)

  // sobre o amortecimento artificial:
  // contraU > 0 -> backward
  // contraU < 0 -> forward

  // calcular por exemplo no ponto i=0 e depois copiar pros demais?
}

/*
- gigiaero, 29/05/2026, 1446 hours
*/
void calc_A_metrics(int m,int n,double A1[m][n],double A2[m][n],double A3[m][n],
                    double x[m][n],double y[m][n],double J[m][n]){
  double ksix,ksiy;
  double etax,etay;

  for(int j=1;j<m-1;j++){
    // ksix = J[j][0]*uniform_scheme_der1_o2_central(m,n,y,0,j,2);
    // ksiy = -J[j][0]*uniform_scheme_der1_o2_central(m,n,x,0,j,2);
    // etax = -J[j][0]*uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j);
    // etay = J[j][0]*uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j);
    metric_terms(&ksix,&ksiy,&etax,&etay,m,n,J,x,y,0,j);

    A1[j][0] = ksix*ksix + ksiy*ksiy;
    A2[j][0] = ksix*etax + ksiy*etay;
    A3[j][0] = etax*etax + etay*etay;

    for(int i=1;i<n-1;i++){
      // ksix = J[j][i]*uniform_scheme_der1_o2_central(m,n,y,i,j,2);
      // ksiy = -J[j][i]*uniform_scheme_der1_o2_central(m,n,x,i,j,2);
      // etax = -J[j][i]*uniform_scheme_der1_o2_central(m,n,y,i,j,1);
      // etay = J[j][i]*uniform_scheme_der1_o2_central(m,n,x,i,j,1);
      metric_terms(&ksix,&ksiy,&etax,&etay,m,n,J,x,y,i,j);

      A1[j][i] = ksix*ksix + ksiy*ksiy;
      A2[j][i] = ksix*etax + ksiy*etay;
      A3[j][i] = etax*etax + etay*etay;
    }
  }
}

/*
- gigiaero, 29/05/2026, 1435 hours
*/
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
- gigiaero, 01/06/2026, 0956 hours
*/
void calc_rho(int m,int n,double rho[m][n-1],double phi[m][n],double A1[m][n],
              double A2[m][n],double A3[m][n]){
  double dphi_dksi,dphi_deta;
  double contraU,contraV;
  double A2_mean;

  for(int j=1;j<m-1;j++){
    for(int i=0;i<n-1;i++){
      dphi_dksi = uniform_scheme_der1_o1_forward(m,n,phi,i,j,1);
      dphi_deta = (uniform_scheme_der1_o2_central(m,n,phi,i,j,2) + 
                   uniform_scheme_der1_o2_central(m,n,phi,i+1,j,2))*.5;

      A2_mean = (A2[j][i+1] + A2[j][i])*.5;
      contraU = calc_contraU(dphi_dksi,dphi_deta,(A1[j][i+1] + A1[j][i])*.5,
                             A2_mean);
      contraV = calc_contraV(dphi_dksi,dphi_deta,A2_mean,
                             (A3[j][i+1] + A3[j][i])*.5);

      rho[j][i] = pow(1 - ((gamma - 1.)/(gamma + 1.))*(contraU*dphi_dksi + 
                  contraV*dphi_deta),1./(gamma - 1.));
    }
  }

  for(int j=0;j<m;j++)
    rho[j][n-1] = rho[j][0];
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

  save_prmtrs_sim(config);
  // salvar aqui parâmetros específicos do fullp ----------

  read_2d_array_from_file(m,n,x,fname_msh_x);
  read_2d_array_from_file(m,n,y,fname_msh_y);

  char *filename = malloc(sizeof(char)*200);
  sprintf(filename,"%s%s",config->casename,"_x_mesh.dat");
  print_2d_array_to_file(m,n,x,filename,0);
  sprintf(filename,"%s%s",config->casename,"_y_mesh.dat");
  print_2d_array_to_file(m,n,y,filename,0);

  calc_J(m,n,J,x,y);
  calc_A_metrics(m,n,A1,A2,A3,x,y,J);

  initialize_fullp(m,n,phi,x,y,fp_prmtrs);

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
      // solve_adi_2d_rectangular_fullp();
      break;
    
    case 3:
      // puts("fullp af2 pending");
      puts("AF2, full potential flow over airfoil");
      // solve_af2_2d_rectangular_fullp(m,n,phi,x,y,config);
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

  sprintf(filename,"%s%s",config->casename,"_u.dat");
  print_2d_array_to_file(m,n,u,filename,0);
  sprintf(filename,"%s%s",config->casename,"_v.dat");
  print_2d_array_to_file(m,n,v,filename,0);
  sprintf(filename,"%s%s",config->casename,"_Ve.dat");
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
    // ksix = J[j][0]*uniform_scheme_der1_o2_central(m,n,y,0,j,2);
    // ksiy = -J[j][0]*uniform_scheme_der1_o2_central(m,n,x,0,j,2);
    // etax = -J[j][0]*uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j);
    // etay = J[j][0]*uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j);

    metric_terms(&ksix,&ksiy,&etax,&etay,m,n,J,x,y,0,j);
    u[j][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j)*ksix + 
              uniform_scheme_der1_o2_central(m,n,phi,0,j,2)*etax;
    v[j][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j)*ksiy + 
              uniform_scheme_der1_o2_central(m,n,phi,0,j,2)*etay;
    
    for(int i=1;i<n-1;i++){
      // ksix = J[j][i]*uniform_scheme_der1_o2_central(m,n,y,i,j,2);
      // ksiy = -J[j][i]*uniform_scheme_der1_o2_central(m,n,x,i,j,2);
      // etax = -J[j][i]*uniform_scheme_der1_o2_central(m,n,y,i,j,1);
      // etay = J[j][i]*uniform_scheme_der1_o2_central(m,n,x,i,j,1);

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

  /* --- */

  // // Domain interior
  // for(int j=1;j<m-1;j++){
  //   u[j][0] = scheme_der1_o2_central_prdc_ksi(m,n,phi,x,j);
  //   v[j][0] = scheme_der1_o2_central_2dxy(m,n,phi,y,0,j,2);
    
  //   for(int i=1;i<n-1;i++){
  //     u[j][i] = scheme_der1_o2_central_2dxy(m,n,phi,x,i,j,1);
  //     v[j][i] = scheme_der1_o2_central_2dxy(m,n,phi,y,i,j,2);
  //   }
  // }

  // // External boundary and airfoil surface
  // for(int i=1;i<n-1;i++){
  //   u[0][i] = scheme_der1_o2_central_2dxy(m,n,phi,x,i,0,1);
  //   v[0][i] = scheme_der1_o2_forward_2dxy(m,n,phi,y,i,0,2);
  //   u[m-1][i] = scheme_der1_o2_central_2dxy(m,n,phi,x,i,m-1,1);
  //   v[m-1][i] = scheme_der1_o2_backward_2dxy(m,n,phi,y,i,m-1,2);
  // }

  // // Corners
  // u[0][0] = scheme_der1_o2_central_prdc_ksi(m,n,phi,x,0);
  // v[0][0] = scheme_der1_o2_forward_2dxy(m,n,phi,y,0,0,2);
  // u[m-1][0] = scheme_der1_o2_central_prdc_ksi(m,n,phi,x,m-1);
  // v[m-1][0] = scheme_der1_o2_backward_2dxy(m,n,phi,y,0,m-1,2);

  /* --- */

  // // Domain interior
  // for(int j=1;j<m-1;j++){
  //   u[j][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j);
  //   v[j][0] = uniform_scheme_der1_o2_central(m,n,phi,0,j,2);
    
  //   for(int i=1;i<n-1;i++){
  //     u[j][i] = uniform_scheme_der1_o2_central(m,n,phi,i,j,1);
  //     v[j][i] = uniform_scheme_der1_o2_central(m,n,phi,i,j,2);
  //   }
  // }

  // // External boundary and airfoil surface
  // for(int i=1;i<n-1;i++){
  //   u[0][i] = uniform_scheme_der1_o2_central(m,n,phi,i,0,1);
  //   v[0][i] = uniform_scheme_der1_o2_forward(m,n,phi,i,0,2);
  //   u[m-1][i] = uniform_scheme_der1_o2_central(m,n,phi,i,m-1,1);
  //   v[m-1][i] = uniform_scheme_der1_o2_backward(m,n,phi,i,m-1,2);
  // }

  // // Corners
  // u[0][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,0);
  // v[0][0] = uniform_scheme_der1_o2_forward(m,n,phi,0,0,2);
  // u[m-1][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,m-1);
  // v[m-1][0] = uniform_scheme_der1_o2_backward(m,n,phi,0,m-1,2);

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
                                double A1[m][n],double A2[m][n],
                                double rho[m][n-1],double C,int i,int j){
  double dphi_dksi,dphi_deta;
  double contraU,rho_til,nu_x;

  if(i == 0){
    dphi_dksi = uniform_scheme_der1_o1_forward(m,n,phi,0,j,1);
    dphi_deta = (uniform_scheme_der1_o2_central(m,n,phi,0,j,2) + 
                 uniform_scheme_der1_o2_central(m,n,phi,1,j,2))*.5;

    contraU = calc_contraU(dphi_dksi,dphi_deta,(A1[j][1] + A1[j][0])*.5,
                           (A2[j][1] + A2[j][0])*.5);

    nu_x = nu_switch(m,n,rho,C,contraU,0,j,1);

    rho_til = rho_coeffs(m,n,rho,nu_x,rs_idx(contraU),0,j,1);

    return rho_til*contraU/((J[j][1] + J[j][0])*.5);
  }
  else{
    dphi_dksi = uniform_scheme_der1_o1_forward(m,n,phi,i,j,1);
    dphi_deta = (uniform_scheme_der1_o2_central(m,n,phi,i,j,2) + 
                 uniform_scheme_der1_o2_central(m,n,phi,i+1,j,2))*.5;

    contraU = calc_contraU(dphi_dksi,dphi_deta,(A1[j][i+1] + A1[j][i])*.5,
                           (A2[j][i+1] + A2[j][i])*.5);

    nu_x = nu_switch(m,n,rho,C,contraU,i,j,1);

    rho_til = rho_coeffs(m,n,rho,nu_x,rs_idx(contraU),i,j,1);

    return rho_til*contraU/((J[j][i+1] + J[j][i])*.5);
  }
}

/*
- gigiaero, 01/06/2026,1332 hours
*/
double L_phi_fullp_der_terms_jh(int m,int n,double phi[m][n],double J[m][n],
                                double A2[m][n],double A3[m][n],
                                double rho[m][n-1],double C,int i,int j){
  double dphi_dksi,dphi_deta;
  double contraV,rho_bar,nu_y;

  if(i == 0){
    dphi_dksi = (uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j) + 
                  uniform_scheme_der1_o2_central_prdc_ksi(m,n,phi,j+1))*.5;
    dphi_deta = uniform_scheme_der1_o1_forward(m,n,phi,0,j,2);

    contraV = calc_contraV(dphi_dksi,dphi_deta,(A2[j+1][0] + A2[j][0])*.5,
                           (A3[j+1][0] + A3[j][0])*.5);

    nu_y = nu_switch(m,n,rho,C,contraV,0,j,2);

    rho_bar = rho_coeffs(m,n,rho,nu_y,rs_idx(contraV),0,j,2);

    return rho_bar*contraV/((J[j+1][0] + J[j][0])*.5);
  }
  else{
    dphi_dksi = (uniform_scheme_der1_o2_central(m,n,phi,i,j,1) + 
                  uniform_scheme_der1_o2_central(m,n,phi,i,j+1,1))*.5;
    dphi_deta = uniform_scheme_der1_o1_forward(m,n,phi,i,j,2);

    contraV = calc_contraV(dphi_dksi,dphi_deta,(A2[j+1][i] + A2[j][i])*.5,
                            (A3[j+1][i] + A3[j][i])*.5);

    nu_y = nu_switch(m,n,rho,C,contraV,i,j,2);

    rho_bar = rho_coeffs(m,n,rho,nu_y,rs_idx(contraV),i,j,2);

    return rho_bar*contraV/((J[j+1][i] + J[j][i])*.5);
  }
}

/*
- gigiaero, 01/06/2026,1315 hours
*/
void L_phi_fullp(int m,int n,double L_phi[m][n],double phi[m][n],
                 double J[m][n],double A1[m][n],double A2[m][n],double A3[m][n],
                 double rho[m][n-1],double C){
  double term_ih,term_ihm1;
  double term_jh,term_jhm1;

  for(int j=1;j<m-1;j++){
    term_ih = L_phi_fullp_der_terms_ih(m,n,phi,J,A1,A2,rho,C,0,j);
    term_ihm1 = L_phi_fullp_der_terms_ih(m,n,phi,J,A1,A2,rho,C,n-2,j);

    term_jh = L_phi_fullp_der_terms_jh(m,n,phi,J,A2,A3,rho,C,0,j);
    term_jhm1 = L_phi_fullp_der_terms_jh(m,n,phi,J,A2,A3,rho,C,0,j-1);

    L_phi[j][0] = (term_ih - term_ihm1) + (term_jh - term_jhm1);

    for(int i=1;i<n-1;i++){
      term_ihm1 = term_ih;
      term_ih = L_phi_fullp_der_terms_ih(m,n,phi,J,A1,A2,rho,C,i,j);
      // term_ihm1 = L_phi_fullp_der_terms_ih(m,n,phi,J,A1,A2,rho,C,i-1,j);

      term_jh = L_phi_fullp_der_terms_jh(m,n,phi,J,A2,A3,rho,C,i,j);
      term_jhm1 = L_phi_fullp_der_terms_jh(m,n,phi,J,A2,A3,rho,C,i,j-1);

      L_phi[j][i] = (term_ih - term_ihm1) + (term_jh - term_jhm1);
    }
  }

  for(int i=1;i<n-1;i++){
    term_ih = L_phi_fullp_der_terms_ih(m,n,phi,J,A1,A2,rho,C,i,0);
    term_ihm1 = L_phi_fullp_der_terms_ih(m,n,phi,J,A1,A2,rho,C,i-1,0);

    term_jh = L_phi_fullp_der_terms_jh(m,n,phi,J,A2,A3,rho,C,i,0);
    // term_jhm1 = L_phi_fullp_der_terms_jh(m,n,phi,J,A2,A3,rho,C,i,-1); !!!!!

    L_phi[0][i] = (term_ih - term_ihm1) + (term_jh - term_jhm1);
  }

  for(int j=0;j<m-1;j++)
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
    // *etax = -J[j][i]*uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j);
    // *etay = J[j][i]*uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j);
  }
  else{
    dy_dksi = uniform_scheme_der1_o2_central(m,n,y,i,j,1);
    dx_dksi = uniform_scheme_der1_o2_central(m,n,x,i,j,1);
    // *etax = -J[j][i]*uniform_scheme_der1_o2_central(m,n,y,i,j,1);
    // *etay = J[j][i]*uniform_scheme_der1_o2_central(m,n,x,i,j,1);
  }

  if(j == 0){
    dy_deta = uniform_scheme_der1_o2_forward(m,n,y,i,j,2);
    dx_deta = uniform_scheme_der1_o2_forward(m,n,x,i,j,2);
    // *ksix = J[j][i]*uniform_scheme_der1_o2_forward(m,n,y,i,j,2);
    // *ksiy = -J[j][i]*uniform_scheme_der1_o2_forward(m,n,x,i,j,2);
  }
  else if(j == m-1){
    dy_deta = uniform_scheme_der1_o2_backward(m,n,y,i,j,2);
    dx_deta = uniform_scheme_der1_o2_backward(m,n,x,i,j,2);
    // *ksix = J[j][i]*uniform_scheme_der1_o2_nackward(m,n,y,i,j,2);
    // *ksiy = -J[j][i]*uniform_scheme_der1_o2_backward(m,n,x,i,j,2);
  }
  else{
    dy_deta = uniform_scheme_der1_o2_central(m,n,y,i,j,2);
    dx_deta = uniform_scheme_der1_o2_central(m,n,x,i,j,2);
    // *ksix = J[j][i]*uniform_scheme_der1_o2_central(m,n,y,i,j,2);
    // *ksiy = -J[j][i]*uniform_scheme_der1_o2_central(m,n,x,i,j,2);
  }

  *ksix = J[j][i]*dy_deta;
  *ksiy = -J[j][i]*dx_deta;
  *etax = -J[j][i]*dy_dksi;
  *etay = J[j][i]*dx_dksi;

  // if(i == 0 && j != 0 && j != m-1){ // periodicity line
  //   *ksix = J[j][i]*uniform_scheme_der1_o2_central(m,n,y,i,j,2);
  //   *ksiy = -J[j][i]*uniform_scheme_der1_o2_central(m,n,x,i,j,2);
  //   *etax = -J[j][i]*uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j);
  //   *etay = J[j][i]*uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j);
  // }
  // else if(i != 0 && j == 0){ // lower edge
  //   *ksix = J[j][i]*uniform_scheme_der1_o2_forward(m,n,y,i,j,2);
  //   *ksiy = -J[j][i]*uniform_scheme_der1_o2_forward(m,n,x,i,j,2);
  //   *etax = -J[j][i]*uniform_scheme_der1_o2_central(m,n,y,i,j,1);
  //   *etay = J[j][i]*uniform_scheme_der1_o2_central(m,n,x,i,j,1);
  // }
  // else if(i != 0 && j == m-1){ // upper edge
  //   *ksix = J[j][i]*uniform_scheme_der1_o2_backward(m,n,y,i,j,2);
  //   *ksiy = -J[j][i]*uniform_scheme_der1_o2_backward(m,n,x,i,j,2);
  //   *etax = -J[j][i]*uniform_scheme_der1_o2_central(m,n,y,i,j,1);
  //   *etay = J[j][i]*uniform_scheme_der1_o2_central(m,n,x,i,j,1);
  // }
  // else if(i == 0 && j == 0){ // low left corner
  //   *ksix = J[j][i]*uniform_scheme_der1_o2_forward(m,n,y,i,j,2);
  //   *ksiy = -J[j][i]*uniform_scheme_der1_o2_forward(m,n,x,i,j,2);
  //   *etax = -J[j][i]*uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j);
  //   *etay = J[j][i]*uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j);
  // }
  // else if(i == 0 && j == m-1){ // high left corner
  //   *ksix = J[j][i]*uniform_scheme_der1_o2_backward(m,n,y,i,j,2);
  //   *ksiy = -J[j][i]*uniform_scheme_der1_o2_backward(m,n,x,i,j,2);
  //   *etax = -J[j][i]*uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j);
  //   *etay = J[j][i]*uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j);
  // }
  // else{ // domain interior
  //   *ksix = J[j][i]*uniform_scheme_der1_o2_central(m,n,y,i,j,2);
  //   *ksiy = -J[j][i]*uniform_scheme_der1_o2_central(m,n,x,i,j,2);
  //   *etax = -J[j][i]*uniform_scheme_der1_o2_central(m,n,y,i,j,1);
  //   *etay = J[j][i]*uniform_scheme_der1_o2_central(m,n,x,i,j,1);
  // }
}

/*
- gigiaero, 01/06/2026, 1048 hours
*/
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

void solve_adi_2d_rectangular_fullp(int m,int n,double phi[m][n],double J[m][n],
                                    double A1[m][n],double A2[m][n],
                                    double A3[m][n],sim_prmtrs *config,
                                    fullp_prmtrs *fp_prmtrs){
  // Solver variables
  double (*L_phi)[n] = calloc(m,sizeof *L_phi);
  double (*Cij)[n] = calloc(m,sizeof *Cij);
  double (*f)[n] = calloc(m,sizeof *f);
  double (*rho)[n-1] = calloc(m,sizeof *rho);

  int iter = 1;

  // a se pensar: calc_contraU e calc_contraV poderiam ser a mesma função

  for(iter;iter<=config->max_iter;iter++){                                
    calc_rho(m,n,rho,phi,A1,A2,A3);
    L_phi_fullp(m,n,L_phi,phi,J,A1,A2,A3,rho,fp_prmtrs->C);

    // calcular resíduos

    // operador N
    // sobre o amortecimento artificial:
    // contraU > 0 -> backward
    // contraU < 0 -> forward

    // atualizar condições de contorno

    // forçar periodicidade aqui
  }

  free(rho);
}