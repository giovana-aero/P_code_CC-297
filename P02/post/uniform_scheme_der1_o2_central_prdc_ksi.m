function val = uniform_scheme_der1_o2_central_prdc_ksi(phi,j)

val = (phi(j,2) - phi(j,end-1))*.5;

end