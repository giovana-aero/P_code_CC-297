clc,clear

address = '../results/';
casename = 'eom';

initial = {'_initial_x.dat','_initial_y.dat'};
result = {'_x_iter_0000000002.dat','_y_iter_0000000002.dat'};

f_plot_mesh(1,address,casename,initial)
f_plot_mesh(2,address,casename,result)

% mesh_x = dlmread([address,casename,'_initial_x.dat']);
% mesh_y = dlmread([address,casename,'_initial_y.dat']);

