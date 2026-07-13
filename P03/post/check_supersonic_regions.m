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

figure(1),clf
config.cb_label = 'q';
config.fontsize = fontsize;
q = dlmread([address,casename,'_q.dat']);
q = double(q >= 1);

if(sum(sum(q)) == 0)
  disp('no supersonic flow in the domain')

else
  % contour(mesh.x,mesh.y,q);
  f_plot_contours(mesh,q,config)

  if save_fig
    set_fontsize_position(fontsize,pos,pba)
  end

end