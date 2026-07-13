clc,clear

address = '../results/';
casename = 'fullp';

save_fig = 0;
base_size = 800;
% fontsize = base_size*12/500;
fontsize = 22;
xlims = [];
ylims = [];

mesh.x = dlmread([address,casename,'_x_mesh.dat']);
mesh.y = dlmread([address,casename,'_y_mesh.dat']);

figure(1),clf,hold on
config.cb_label = 'p';
config.fontsize = fontsize;
f_plot_contours(mesh,dlmread([address,casename,'_p.dat']),config)
plot(mesh.x(1,:),mesh.y(1,:),'k')
% colormap(purples_discrete())

if save_fig
  set_fontsize_position(fontsize,pos,pba)
end
