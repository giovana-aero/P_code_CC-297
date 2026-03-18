

/*
Scheme: second derivative, second order, central
- gigiaero, 13/03/2026, 2315 hours
*/
void scheme_der2_o2_central(double phi_ip1,double phi_i,double phi_im1,
                            double x_ip1,double x_i,double x_im1){

  return 2./(x_ip1 - x_im1)*((phi_ip1 - phi_i)/(x_ip1 - x_i) - 
         (phi_i - phi_im1)/(x_i - x_im1));

}