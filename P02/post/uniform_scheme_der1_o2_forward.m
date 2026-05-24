function val = uniform_scheme_der1_o2_forward(phi,i,j,axis)

if axis == 1
  val = (-3*phi(j,i) + 4*phi(j,i+1) - phi(j,i+2))*.5;
elseif axis == 2
  val = (-3*phi(j,i) + 4*phi(j+1,i) - phi(j+2,i))*.5;
end

end