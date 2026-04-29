clc,clear

address = '../results/';
casename = 'eom';

iter = -1;

initial = {'_initial_x.dat','_initial_y.dat'};

if iter < 0
  iter_list_x = dir([address,casename,'_x_iter_*']);
  iter_list_y = dir([address,casename,'_y_iter_*']);
  result = {iter_list_x(end).name,iter_list_y(end).name};
else
  result = {sprintf("%s_x_iter_%010d.dat",casename,iter),...
            sprintf("%s_y_iter_%010d.dat",casename,iter)};
end

f_plot_mesh(1,address,casename,initial)
f_plot_mesh(2,address,'',result)

% mesh_x = dlmread([address,casename,'_initial_x.dat']);
% mesh_y = dlmread([address,casename,'_initial_y.dat']);

