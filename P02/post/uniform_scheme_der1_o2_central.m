function val = uniform_scheme_der1_o2_central(phi,i,j,axis)

if axis == 1
  val = (phi(j,i+1) - phi(j,i-1))*.5;
elseif axis == 2
  val = (phi(j+1,i) - phi(j-1,i))*.5;
end

end