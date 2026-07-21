clc,clear

save_fig = 1;
fontsize = 16;

casename = 'fullp';

% case 1 results
address = '../cases_results/case1/';
pos = [1000,500];
pba = [2,1];
xlims = [-1.5,2.5];
ylims = [-1,1];
filename = '../../../P03/figuras/resultados/case1_cp_contours';

% % case 2 results
% address = '../cases_results/case2/';
% pos = [1000,500];
% pba = [2,1];
% xlims = [-1.5,2.5];
% ylims = [-1,1];
% filename = '../../../P03/figuras/resultados/case2_cp_contours';

% % naca 0012 transonic
% address = '../transonic_validations/naca0012_M08_mesh4/';
% pos = [1000,500];
% pba = [2,1];
% xlims = [-1.5,2.5];
% ylims = [-1,1];
% filename = '../../../P03/figuras/resultados/naca0012_M08_mesh4_cp_contours';

% % symmetrical supercritical transonic
% address = '../transonic_validations/supercritical_tm-x-1831_M08_mesh3/';
% pos = [1000,500];
% pba = [2,1];
% xlims = [-1.5,2.5];
% ylims = [-1,1];
% filename = '../../../P03/figuras/resultados/supercritical_tm-x-1831_M08_mesh3_cp_contours';

mesh.x = dlmread([address,casename,'_x_mesh.dat']);
mesh.y = dlmread([address,casename,'_y_mesh.dat']);

figure(1),clf,hold on
config.cb_label = 'C_p';
config.fontsize = fontsize;
f_plot_contourfs(mesh,dlmread([address,casename,'_cp.dat']),config)
plot(mesh.x(1,:),mesh.y(1,:),'k','linewidth',.5)
xylabel_latex('x','y')

if ~isempty(xlims),xlim(xlims),end
if ~isempty(ylims),ylim(ylims),end

% figure(2),clf,hold on
% f_plot_contours(mesh,dlmread([address,casename,'_cp.dat']),config)
% plot(mesh.x(1,:),mesh.y(1,:),'k')
% xylabel_latex('x','y')

if save_fig
  figure(1)
  set_fontsize_position(fontsize,pos,pba)
  set_rndrr_painters()
  save_fig_eps(filename);

  % figure(2)
  % set_fontsize_position(fontsize,pos,pba)
  % set_rndrr_painters()
end
