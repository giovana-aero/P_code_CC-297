clc,clear

address = '../results/';
casename = 'fullp';

save_fig = 0;
base_size = 800;
% fontsize = base_size*12/500;
fontsize = 22;
xlims = [];
ylims = [];

mesh.x = dlmread([address,casename,'_x_mesh_jh.dat']);
mesh.y = dlmread([address,casename,'_y_mesh_jh.dat']);

figure(1),clf
config.cb_label = '\phi';
config.fontsize = fontsize;
f_plot_contours(mesh,dlmread([address,'fullp_potential_jh.dat']),config)
% contourf(mesh.x,mesh.y,dlmread([address,result]))
% grid on

if save_fig
  set_fontsize_position(fontsize,pos,pba)
end
