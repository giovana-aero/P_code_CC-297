clc,clear

address = '../results/';
casename = 'eom';

iter = -1;

if iter < 0
  iter_list_x = dir([address,casename,'_x_iter_*']);
  iter_list_y = dir([address,casename,'_y_iter_*']);
  result = {iter_list_x(end).name,iter_list_y(end).name};
else
  result = {sprintf("%s_x_iter_%010d.dat",casename,iter),...
            sprintf("%s_y_iter_%010d.dat",casename,iter)};
end

initial = {'_x_initial.dat','_y_initial.dat'};
f_orthogonality_eom(address,casename,initial),disp('')
f_orthogonality_eom(address,'',result)