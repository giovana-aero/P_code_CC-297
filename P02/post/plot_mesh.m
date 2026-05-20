clc,clear

address = '../results/';
casename = 'eom';

save_fig = 0;
base_size = 800;
fontsize = base_size*12/500;
xlims = [-1,2];
% xlims = [];
ylims = [-1,1];

iter = -1;

initial = {'_x_initial.dat','_y_initial.dat'};

if iter < 0
  iter_list_x = dir([address,casename,'_x_iter_*']);
  iter_list_y = dir([address,casename,'_y_iter_*']);
  result = {iter_list_x(end).name,iter_list_y(end).name};
else
  result = {sprintf("%s_x_iter_%010d.dat",casename,iter),...
            sprintf("%s_y_iter_%010d.dat",casename,iter)};
end

figure(1),clf
f_plot_mesh(1,address,casename,initial)
if ~isempty(xlims),xlim(xlims),end
figure(2),clf
f_plot_mesh(2,address,'',result)
if ~isempty(xlims),xlim(xlims),end

% mesh_x = dlmread([address,casename,'_initial_x.dat']);
% mesh_y = dlmread([address,casename,'_initial_y.dat']);

if save_fig
  fig_height = base_size*ylims/xlims;
  fig_width = base_size;
  pba = [fig_width,fig_height]/fig_width;
  pos = [fig_width,fig_height];
  
  figure(1)
  set_rndrr_painters()
  set_fontsize_position(fontsize,pos,pba)
  xylabel_latex('x','y')
  save_fig_eps('mesh_initial.eps');

end