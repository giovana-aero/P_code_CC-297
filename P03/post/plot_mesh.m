clc,clear
addpath('../../P02/post')

casename = 'fullp';

save_fig = 0;
fontsize = 14;

% % results
% address = '../results/';
% pos = [500,500];
% pba = [1,1];
xlims = [];
ylims = [];

% mesh convergence, mesh 4, eom
address = '../mesh_conv/mesh4/';
pos = [1000,500];
pba = [2,1];
% xlims = [-1,1]*4e-3;
% ylims = [-1,1]*4e-3;
filename = '../../../P03/figuras/testes/mesh_conv_mesh4_eom_LE';

% % mesh convergence, mesh 4, pom
% address = '../mesh_conv/mesh4_pom/';
% pos = [1000,500];
% pba = [2,1];
% xlims = [-1,1]*4e-3;
% ylims = [-1,1]*4e-3;
% filename = '../../../P03/figuras/testes/mesh_conv_mesh4_pom_LE';

mesh = {'_x_mesh.dat','_y_mesh.dat'};

figure(1),clf
hold on
f_plot_mesh(1,address,casename,mesh)
if ~isempty(xlims),xlim(xlims),end
if ~isempty(ylims),ylim(ylims),end
xylabel_latex('x','y')

if save_fig
  set_rndrr_painters()
  set_fontsize_position(fontsize,pos,pba)
  save_fig_eps(filename);
end