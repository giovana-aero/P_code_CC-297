clc,clear

address = '../results/';
casename = 'fullp';

save_fig = 0;
base_size = 800;
% fontsize = base_size*12/500;
fontsize = 22;
xlims = [];
ylims = [];

iter = -1;

if iter < 0
  iter_list = dir([address,casename,'_iter_*']);
  result = iter_list(end).name;
else
  result = sprintf("%s_iter_%010d.dat",casename,iter);
end

mesh.x = dlmread([address,casename,'_x_mesh.dat']);
mesh.y = dlmread([address,casename,'_y_mesh.dat']);

figure(1),clf
config.cb_label = '\phi';
config.fontsize = fontsize;
f_plot_contours(mesh,dlmread([address,result]),config)
% contourf(mesh.x,mesh.y,dlmread([address,result]))
% grid on

if save_fig
  set_fontsize_position(fontsize,pos,pba)
end
