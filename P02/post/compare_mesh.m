clc,clear

address = {'../results/','../results_unmod/'};
% address = {'../results/','../results/'};
initial = {'_x_initial.dat','_y_initial.dat'};
casename = 'eom';

save_fig = 1;
base_size = 800;
fontsize = 24;
xlims = [-.2,1.2];
ylims = [-1,1]*.2;
% xlims = [-.04,.1]; % le zoom
% ylims = [-1,1]*.06; 
% xlims = [.9,1.04]; % te zoom
% ylims = [-1,1]*.06;
% xlims = [];
% ylims = [];

iter = -1;

result = cell(1,2);
for i = 1:2
  if iter < 0
    iter_list_x = dir([address{i},casename,'_x_iter_*']);
    iter_list_y = dir([address{i},casename,'_y_iter_*']);
    result{i} = {iter_list_x(end).name,iter_list_y(end).name};
  else
    result{i} = {sprintf("%s_x_iter_%010d.dat",casename,iter),...
                 sprintf("%s_y_iter_%010d.dat",casename,iter)};
  end
end

% figure(1),clf
% f_plot_mesh(1,address{1},casename,initial,'-')
% f_plot_mesh(1,address{2},'',result{1},'--')

figure(1),clf
f_plot_mesh(1,address{1},'',result{1})
f_plot_mesh(1,address{2},'',result{2},'--')

if ~isempty(xlims),xlim(xlims),end
if ~isempty(ylims),ylim(ylims),end

if save_fig
  % fig_height = base_size*diff(ylims)/diff(xlims);
  fig_height = base_size;
  fig_width = base_size;
  pba = [fig_width,fig_height]/fig_width;
  pos = [fig_width,fig_height];
  
  figure(1)
  set_rndrr_painters()
  set_fontsize_position(fontsize,pos,pba)
  xylabel_latex('x','y')
  save_fig_eps('mesh_comparison.eps');
end