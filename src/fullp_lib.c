#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"../include/2d_arrays.h"
#include"../include/eom_lib.h"
#include"../include/num_methods.h"

#define gamma 1.4

/*
- gigiaero, 29/05/2026, 1446 hours
*/
void calc_A_metrics(int m,int n,double A1[m][n],double A2[m][n],double A3[m][n],
                    double x[m][n],double y[m][n],double J[m][n]){
  double ksix,ksiy;
  double etax,etay;

  for(int j=1;j<m-1;j++){
    ksix = J[j][0]*uniform_scheme_der1_o2_central(m,n,y,0,j,2);
    ksiy = -J[j][0]*uniform_scheme_der1_o2_central(m,n,x,0,j,2);
    etax = -J[j][0]*uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j);
    etay = J[j][0]*uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j);

    A1[j][0] = ksix*ksix + ksiy*ksiy;
    A2[j][0] = ksix*etax + ksiy*etay;
    A3[j][0] = etax*etax + etay*etay;

    for(int i=1;i<n-1;i++){
      ksix = J[j][i]*uniform_scheme_der1_o2_central(m,n,y,i,j,2);
      ksiy = -J[j][i]*uniform_scheme_der1_o2_central(m,n,x,i,j,2);
      etax = -J[j][i]*uniform_scheme_der1_o2_central(m,n,y,i,j,1);
      etay = J[j][i]*uniform_scheme_der1_o2_central(m,n,x,i,j,1);

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
  // switch(axis){ 
  //   case 1: // i+1/2,j
  return A2_val*dphi_dksi + A3_val*dphi_deta;
  //     break;

  //   case 2: // i,j+1/2
  //     return 
  //     break;

  //   default:
  //     puts("calc_contraV: invalid axis");
  //     exit(15);

  // }
}

// /*
// - gigiaero, 31/05/2026, 1119 hours
// */
// double calc_rho(double contraU,double contraV,double dphi_dksi,
//                 double dphi_deta){
//   return pow(1 - ((gamma - 1.)/(gamma + 1.))*(contraU*dphi_dksi + 
//              contraV*dphi_deta),1./(gamma - 1.));
// }
// double calc_rho(double *contraU,double *contraV,int m,int n,double phi[m][n],
//                 double A1[m][n],double A2[m][n],double A3[m][n],int i,int j){
//   double dphi_dksi,dphi_deta;
//   // double contraU,contraV;

//   dphi_dksi = uniform_scheme_der1_o1_forward(m,n,phi,i,j,1);
//   dphi_deta = (uniform_scheme_der1_o2_central(m,n,phi,i,j,2) + 
//                 uniform_scheme_der1_o2_central(m,n,phi,i+1,j,2))*.5;

//   *contraU = calc_contraU(dphi_dksi,dphi_deta,(A1[j][i+1] + A1[j][i])*.5,
//                           (A2[j][i+1] + A2[j][i])*.5);
//   *contraV = calc_contraV(dphi_dksi,dphi_deta,(A2[j][i+1] + A2[j][i])*.5,
//                           (A3[j][i+1] + A3[j][i])*.5);

//   return pow(1 - ((gamma - 1.)/(gamma + 1.))*((*contraU)*dphi_dksi + 
//              (*contraV)*dphi_deta),1./(gamma - 1.));
// }

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
void evaluate_delta_form_fullp(int m,int n,sim_prmtrs *config,char *fname_msh_x,
                               char *fname_msh_y){
  double (*x)[n] = calloc(m,sizeof *x);
  double (*y)[n] = calloc(m,sizeof *y);
  double (*phi)[n] = calloc(m,sizeof *phi);

  read_2d_array_from_file(m,n,x,fname_msh_x);
  read_2d_array_from_file(m,n,y,fname_msh_y);

  // condição
  // inicial
  // aqui

  if(config->save_i_c){
    char *filename = malloc(sizeof(char)*200);

    sprintf(filename,"%s_iter_%010d.dat",config->casename,0);
    print_2d_array_to_file(m,n,phi,filename,0);

    free(filename);
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
      solve_af2_2d_rectangular_fullp(m,n,phi,x,y,config);
      break;
      
    default:
      puts("evaluate_delta_form_eom: Invalid Ntype");
      exit(32);
  }

  /*
  calcular aqui: velocidades e densidades
  */



  // /* apagar depois */
  // char *filename = malloc(sizeof(char)*200);
  // sprintf(filename,"%s%s",config->casename,"_x_initial.dat");
  // print_2d_array_to_file(m,n,x,filename,0);
  // sprintf(filename,"%s%s",config->casename,"_y_initial.dat");
  // print_2d_array_to_file(m,n,y,filename,0);
  // free(filename);
  // /* apagar depois */

  free(x);
  free(y);
  free(phi);
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

void solve_af2_2d_rectangular_fullp(int m,int n,double phi[m][n],double x[m][n],
                                    double y[m][n],sim_prmtrs *config){
  double (*J)[n] = calloc(m,sizeof *J);
  double (*A1)[n] = calloc(m,sizeof *A1);
  double (*A2)[n] = calloc(m,sizeof *A2);
  double (*A3)[n] = calloc(m,sizeof *A3);
  double (*rho)[n-1] = calloc(m,sizeof *rho);

  int iter = 1;

  // a se pensar: calc_contraU e calc_contraV poderiam ser a mesma função

  for(iter;iter<=config->max_iter;iter++){
    calc_J(m,n,J,x,y);
    calc_A_metrics(m,n,A1,A2,A3,x,y,J);                                        
    calc_rho(m,n,rho,phi,A1,A2,A3);


    // forçar periodicidade aqui

  }

  free(J);
  free(A1);
  free(A2);
  free(A3);
  free(rho);
}