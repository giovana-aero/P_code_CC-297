clc,clear

save_fig = 0;
fontsize = 18;

casename = 'fullp';

% % case 1 results
% address = '../cases_results/case1/';
% pos = [1000,500];
% pba = [2,1];
% xlims = [-1.5,2.5];
% ylims = [-1,1]*1.5;
% filename = '../../../P03/figuras/resultados/case1_rho_contours';

% case 2 results
address = '../cases_results/case2/';
pos = [1000,500];
pba = [2,1];
xlims = [-1.5,2.5];
ylims = [-1,1]*1.5;
filename = '../../../P03/figuras/resultados/case2_rho_contours';



mesh.x = dlmread([address,casename,'_x_mesh.dat']);
mesh.y = dlmread([address,casename,'_y_mesh.dat']);

figure(1),clf,hold on
config.cb_label = '\rho';
config.fontsize = fontsize;
f_plot_contourfs(mesh,dlmread([address,casename,'_rho.dat']),config)
plot(mesh.x(1,:),mesh.y(1,:),'k','linewidth',.5)
xylabel_latex('x','y')

if save_fig
  set_fontsize_position(fontsize,pos,pba)
  set_rndrr_painters()
  save_fig_eps(filename);
end
