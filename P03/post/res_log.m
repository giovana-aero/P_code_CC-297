clc,clear

casename = 'fullp';
ref_val = 1e-6;
save_figs = 1;
cmap = slanCM('bupu');

% % results
% address_list = {'../results/'};
% lgds = {'../results/'};
% linestyles = {'-'};
% xlims = [];
% ylims = [];
% filename = [];
% colors = {[0,0,0]};

% % w
% address_list = {'../optimizations/base/';
%                 '../optimizations/omega/w09/';
%                 '../optimizations/omega/w15/';
%                 '../optimizations/omega/w20/';
%                 '../optimizations/omega/w25/';};
% lgds = {'$\omega=1.0$','$\omega=0.9$','$\omega=1.5$','$\omega=2.0$','$\omega=2.5$'};
% linestyles = {'-','-','--','-','--'};
% xlims = [-100,900];
% ylims = [5e-7,1e-2];
% filename = '../../../P03/figuras/testes/w_optimization_af2';
% colors = {[0,0,0],cmap(210,:),cmap(210,:),cmap(100,:),cmap(100,:)};

% % a 
% address_list = {'../optimizations/base/';
%                 '../optimizations/alpha/a1e-2/';
%                 '../optimizations/alpha/a1e-1/';
%                 '../optimizations/alpha/a1e1/';
%                 '../optimizations/alpha/aseq_suggestion/';
%                 '../optimizations/alpha/aseq_opt/'};
% lgds = {'$\alpha=1$','$\alpha=0.01$','$\alpha=0.1$','$\alpha=10$','Seq. sugerida','$\alpha_L=0.1,\alpha_H=0.2$'};
% linestyles = {'-','-','--',':','-','--'};
% xlims = [-.05,.6]*1e4;
% ylims = [5e-7,100];
% filename = '../../../P03/figuras/testes/a_optimization_af2';
% colors = {[0,0,0],cmap(210,:),cmap(210,:),cmap(210,:),cmap(100,:),cmap(100,:)};

% % optimized
% address_list = {'../optimizations/base/';
%                 '../optimizations/optimized/'};
% lgds = {'Base','Otimizado'};
% linestyles = {'-','-'};
% xlims = [-100,800];
% ylims = [5e-7,1e-2];
% filename = '../../../P03/figuras/testes/af2_optimized';
% colors = {[0,0,0],cmap(210,:)};

% % mesh convergence
% address_list = {'../mesh_conv/mesh1/';
%                 '../mesh_conv/mesh2/';
%                 '../mesh_conv/mesh3/';
%                 '../mesh_conv/mesh4/'};
% lgds = {'Malha 1','Malha 2','Malha 3','Malha 4'};
% linestyles = {'-','--','-',':'};
% xlims = [-10,300];
% ylims = [5e-7,2e-2];
% filename = '../../../P03/figuras/testes/af2_mesh_convergence_res';
% colors = {cmap(210,:),cmap(210,:),[0,0,0],cmap(210,:)};

% % mesh convergence, mesh 4 eom/pom comparison
% address_list = {'../mesh_conv/mesh4/';
%                 '../mesh_conv/mesh4_pom/'};
% lgds = {'Malha 4 (E.)','Malha 4 (P.)'};
% linestyles = {':','-'};
% xlims = [-50,800];
% ylims = [5e-7,2e-2];
% filename = '../../../P03/figuras/testes/af2_mesh_convergence_res_mesh4_eom_pom';
% colors = {cmap(210,:),[0,0,0]};

% mach sweep, biair
address_list = {'../mach_sweep/biair_M010/';
                '../mach_sweep/biair_M020/';
                '../mach_sweep/biair_M030/';
                '../mach_sweep/biair_M040/';
                '../mach_sweep/biair_M050/';
                '../mach_sweep/biair_M060/';
                '../mach_sweep/biair_M070/';
                '../mach_sweep/biair_M080/';
                '../mach_sweep/biair_M090/';
                '../mach_sweep/biair_M100/';
                '../mach_sweep/biair_M110/'};
lgds = {'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','1.1'};
% linestyles = {'-','-','-','-','-','-','-','-','-','-','-'};
linestyles = {'-','--',':','-.','-','--',':','-.','-','--',':','-.'};
xlims = [-50,500];
ylims = [5e-7,5];
filename = '../../../P03/figuras/testes/mach_sweep_biair_res';
% colors = {[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]};
% colors = config_colors_cmap(length(address_list),cmap,1);
% colors = {[0.3020 0 0.2941],[0.3020 0 0.2941],[0.3020 0 0.2941],[0.3020 0 0.2941],[0.5490 0.4533 0.7106],[0.5490 0.4533 0.7106],[0.5490 0.4533 0.7106],[0.5490 0.4533 0.7106],[0.8048 0.8697 0.9256],[0.8048 0.8697 0.9256],[0.8048 0.8697 0.9256]};
colors = {[0,0,0],[0,0,0],[0,0,0],[0,0,0],cmap(210,:),cmap(210,:),cmap(210,:),cmap(210,:),cmap(100,:),cmap(100,:),cmap(100,:),cmap(100,:),};

% % mach sweep, naca0010
% address_list = {'../mach_sweep/naca0010_M010/';
%                 '../mach_sweep/naca0010_M020/';
%                 '../mach_sweep/naca0010_M030/';
%                 '../mach_sweep/naca0010_M040/';
%                 '../mach_sweep/naca0010_M050/';
%                 '../mach_sweep/naca0010_M060/';
%                 '../mach_sweep/naca0010_M070/';
%                 '../mach_sweep/naca0010_M080/';
%                 '../mach_sweep/naca0010_M090/';
%                 '../mach_sweep/naca0010_M100/';
%                 '../mach_sweep/naca0010_M110/'};
% lgds = {'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','1.1'};
% % linestyles = {'-','-','-','-','-','-','-','-','-','-','-'};
% linestyles = {'-','--',':','-.','-','--',':','-.','-','--',':','-.'};
% xlims = [-50,500];
% ylims = [5e-7,5];
% filename = '../../../P03/figuras/testes/mach_sweep_naca0010_res';
% colors = {[0,0,0],[0,0,0],[0,0,0],[0,0,0],cmap(210,:),cmap(210,:),cmap(210,:),cmap(210,:),cmap(100,:),cmap(100,:),cmap(100,:),cmap(100,:),};

fontsize = 14;
lnwdth = 1.5;
pos = [500,500];
pba = [1,1];

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure(1),clf
f_res_log(address_list,casename,colors,linestyles,lgds,ref_val,fontsize,lnwdth,pos,pba,xlims,ylims,save_figs,filename)
