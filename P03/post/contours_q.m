clc,clear

save_fig = 0;
fontsize = 16;

casename = 'fullp';

% results
address = '../results/';
xlims = [];
ylims = [];

% % mesh_conv, mesh4 eom
% address = '../mesh_conv/mesh4/';
% pos = [1000,500];
% pba = [2,1];
% % xlims = [-1,1]*1e-3; % leading edge
% % ylims = [-1,1]*1e-3;
% % filename = '../../../P03/figuras/testes/mesh_conv_mesh4_eom_LE_q_contours';
% xlims = [-1,1]*1e-3+1; % trailing edge
% ylims = [-1,1]*1e-3;
% filename = '../../../P03/figuras/testes/mesh_conv_mesh4_eom_TE_q_contours';

% % case 1 results
% address = '../cases_results/case1/';
% pos = [1000,500];
% pba = [2,1];
% xlims = [-1.5,2.5];
% ylims = [-1,1];
% filename = '../../../P03/figuras/resultados/case1_q_contours';

% % case 2 results
% address = '../cases_results/case2/';
% pos = [1000,500];
% pba = [2,1];
% xlims = [-1.5,2.5];
% ylims = [-1,1];
% filename = '../../../P03/figuras/resultados/case2_q_contours';

% % naca 0012 transonic
% address = '../transonic_validations/naca0012_M08_mesh4/';
% pos = [1000,500];
% pba = [2,1];
% xlims = [-1.5,2.5];
% ylims = [-1,1];
% filename = '../../../P03/figuras/resultados/naca0012_M08_mesh4_q_contours';

% % symmetrical supercritical transonic
% address = '../transonic_validations/supercritical_tm-x-1831_M08_mesh3/';
% pos = [1000,500];
% pba = [2,1];
% xlims = [-1.5,2.5];
% ylims = [-1,1];
% filename = '../../../P03/figuras/resultados/supercritical_tm-x-1831_M08_mesh3_q_contours';

mesh.x = dlmread([address,casename,'_x_mesh.dat']);
mesh.y = dlmread([address,casename,'_y_mesh.dat']);

data = dlmread([address,casename,'_q.dat']);

figure(1),clf,hold on
config.cb_label = 'q';
config.fontsize = fontsize;
f_plot_contourfs(mesh,data,config)
plot(mesh.x(1,:),mesh.y(1,:),'k','linewidth',.5)
xylabel_latex('x','y')
contour(mesh.x,mesh.y,data,[0,1],'color',[1,1,1],'linestyle','--') % supersonic regions

if ~isempty(xlims),xlim(xlims),end
if ~isempty(ylims),ylim(ylims),end

if save_fig
  figure(1)
  set_fontsize_position(fontsize,pos,pba)
  set_rndrr_painters()
  save_fig_eps(filename);

  % figure(2)
  % set_fontsize_position(fontsize,pos,pba)
  % set_rndrr_painters()
end
