clc,clear

address = '../results/';
casename = 'fullp';

cmap = slanCM('bupu');

% % results
% address_list = {'../results/'};
% lgds = {'../results/'};
% linestyles = {'-'};
% colors = {[0,0,0]};
% filename = '[]';

% % case 1 and 2 results
% address_list = {'../cases_results/case1/';
%                 '../cases_results/case2/'};
% lgds = {'Caso 1','Caso 2'};
% linestyles = {'-','-','-',':'};
% colors = {cmap(210,:),cmap(100,:)}
% filename = '../../../P03/figuras/resultados/case1_case2_cp_airfoil';

% % mesh convergence
% address_list = {'../mesh_conv/mesh1/';
%                 '../mesh_conv/mesh2/';
%                 '../mesh_conv/mesh3/';
%                 '../mesh_conv/mesh4/'};
% lgds = {'Malha 1','Malha 2','Malha 3','Malha 4'};
% linestyles = {'-','--','-',':'};
% colors = {cmap(210,:),cmap(210,:),[0,0,0],cmap(210,:)};
% filename = '../../../P03/figuras/testes/af2_mesh_conv_cp';

% % mesh convergence, mesh 4 eom/pom comparison
% address_list = {'../mesh_conv/mesh4/';
%                 '../mesh_conv/mesh4_pom/'};
% lgds = {'Malha 4 (E.)','Malha 4 (P.)'};
% linestyles = {':','-'};
% xlims = [-10,1100];
% ylims = [5e-7,2e-2];
% filename = '../../../P03/figuras/testes/af2_mesh_convergence_cp_mesh4_eom_pom';
% % colors = {cmap(210,:),cmap(100,:)};
% colors = {cmap(210,:),[0,0,0]};

% % outer boundary effect
% address_list = {'../outer_boundary_effects/M084_OR1/';
%                 '../outer_boundary_effects/M084_OR2/';
%                 '../outer_boundary_effects/M094_OR1/';
%                 '../outer_boundary_effects/M094_OR2/'};
% lgds = {'Mach .84, R1','Mach .84, R2','Mach .94, R1','Mach .94, R2'};
% linestyles = {'-','--','-.',':'};
% xlims = [];
% ylims = [];
% filename = '../../../P03/figuras/testes/outer_boundary_cp';
% colors = {cmap(210,:),cmap(100,:),cmap(210,:),cmap(100,:)};

% % mach sweep, biair
% address_list = {'../mach_sweep/biair_M010/';
%                 '../mach_sweep/biair_M020/';
%                 '../mach_sweep/biair_M030/';
%                 '../mach_sweep/biair_M040/';
%                 '../mach_sweep/biair_M050/';
%                 '../mach_sweep/biair_M060/';
%                 '../mach_sweep/biair_M070/';
%                 '../mach_sweep/biair_M080/';
%                 '../mach_sweep/biair_M090/';
%                 '../mach_sweep/biair_M100/'};
% lgds = {'Mach .1','Mach .2','Mach .3','Mach .4','Mach .5','Mach .6','Mach .7','Mach .8','Mach .9','Mach 1.0'};
% % linestyles = {'-','-','-','-','-','-','-','-','-','-','-'};
% linestyles = {'-','--',':','-.','-','--',':','-.','-','--',':','-.'};
% xlims = [];
% ylims = [];
% filename = '../../../P03/figuras/testes/mach_sweep_biair_cp_airfoil';
% % colors = {[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]};
% colors = config_colors_cmap(length(address_list),cmap,1);

% mach sweep, biair
address_list = {'../transonic_validations/naca0012_M08_mesh4/';
                '../transonic_validations/supercritical_tm-x-1831_M08_mesh3/'};
lgds = {'NACA 0012 (PC)','PS TM X-1831 (PC)','NACA 0012 (EXP)','PS TM X-1831 (EXP)'};
linestyles = {'-','-'};
xlims = [];
ylims = [];
filename = '../../../P03/figuras/resultados/naca0012_supercritical_foil_M08_comparison';
colors = {cmap(210,:),cmap(100,:)};


save_fig = 1;
fontsize = 14;
lnwdth = 1.5;
pos = [1000,500];
pba = [2,1];

figure(1),clf,grid on,hold on
for i = 1:length(address_list)
  % mesh_x = readmatrix([address_list{i},casename,'_x_mesh.dat'],'filetype','delimitedtext');
  % cp = readmatrix([address_list{i},casename,'_cp.dat'],'filetype','delimitedtext');
  mesh_x = dlmread([address_list{i},casename,'_x_mesh.dat']);
  cp = dlmread([address_list{i},casename,'_cp.dat']);
  
  chord_nump = (size(cp,2) + 1)/2;

  plot(mesh_x(1,1:chord_nump),cp(1,1:chord_nump),'color',colors{i},'linestyle',linestyles{i},'displayname',lgds{i},'linewidth',1.5)
  % plot(mesh_x(1,chord_nump:end),cp(1,chord_nump:end),'color',colors{i},'linestyle',linestyles{i})
end


if contains(address_list{1},'transonic_validations')
  exp_data = dlmread('../../../P03/fontes_validacoes/tm-x-1831_fig13_naca0012.dat');
  plot(exp_data(:,1),exp_data(:,2),'.','color',cmap(210,:),'markersize',20,'displayname',lgds{3})
  % plot(exp_data(:,1),exp_data(:,2),'^k','markersize',20,'displayname',lgds{3})
  exp_data = dlmread('../../../P03/fontes_validacoes/tm-x-1831_fig13_supercritical.dat');
  plot(exp_data(:,1),exp_data(:,2),'.','color',cmap(100,:),'markersize',20,'displayname',lgds{4})
  % plot(exp_data(:,1),exp_data(:,2),'vk','markersize',20,'displayname',lgds{4})
end

xylabel_latex('x/c','C_p')
lgd_latex(legend())
inv_yaxis()


% mesh.x = dlmread([address,casename,'_x_mesh.dat']);
% mesh.y = dlmread([address,casename,'_y_mesh.dat']);

% cp = dlmread([address,casename,'_cp.dat']);

% chord_nump = (size(cp,2) + 1)/2;
% figure(1),clf
% % plot(mesh.x(1,1:chord_nump),mesh.y(1,1:chord_nump)),hold on,grid on
% % plot(mesh.x(1,chord_nump:end),mesh.y(1,chord_nump:end))
% plot(mesh.x(1,1:chord_nump),cp(1,1:chord_nump),'k'),hold on,grid on
% plot(mesh.x(1,chord_nump:end),cp(1,chord_nump:end),'k')
% xylabel_latex('x/c','C_p')
% inv_yaxis();

xlim([-.1,1.1])
ylim([-1.5,1.5])

set_fontsize_position(fontsize,pos,pba);

if save_fig
  set_rndrr_painters()
  save_fig_eps(filename)
end

