clc,clear


casename = 'eom';
ref_val = 1e-6;
save_figs = 1;
cmap = slanCM('bupu');

% % w - slor
% address_list = {'../optimizations/base/slor/';
%                 '../optimizations/omega/slor/w09/';
%                 '../optimizations/omega/slor/w15/';
%                 '../optimizations/omega/slor/w20/';
%                 '../optimizations/omega/slor/w25/';
%                 '../optimizations/omega/slor/w30/'};
% lgds = {'$\omega=1.0$','$\omega=0.9$','$\omega=1.5$','$\omega=2.0$','$\omega=2.5$','$\omega=3.0$'};
% linestyles = {'-','-','-.','--','-','-.'};
% xlims = [-.2,3.5]*1e4;
% ylims = [5e-7,.1];
% filename = '../../../P02/figuras/testes/w_optimization_slor';
% colors = {[0,0,0],cmap(210,:),cmap(210,:),cmap(210,:),cmap(100,:),cmap(100,:)};

% % w - adi
% address_list = {'../optimizations/base/adi/';
%                 '../optimizations/omega/adi/w09/';
%                 '../optimizations/omega/adi/w15/';
%                 '../optimizations/omega/adi/w20/';
%                 '../optimizations/omega/adi/w25/';
%                 '../optimizations/omega/adi/w30/'};
% lgds = {'$\omega=1.0$','$\omega=0.9$','$\omega=1.5$','$\omega=2.0$','$\omega=2.5$','$\omega=3.0$'};
% linestyles = {'-','-','-.','--','-','-.'};
% xlims = [-.2,3.]*1e4;
% ylims = [5e-7,.1];
% filename = '../../../P02/figuras/testes/w_optimization_adi';
% colors = {[0,0,0],cmap(210,:),cmap(210,:),cmap(210,:),cmap(100,:),cmap(100,:)};

% % w - af2
% address_list = {'../optimizations/base/af2/';
%                 '../optimizations/omega/af2/w09/';
%                 '../optimizations/omega/af2/w15/';
%                 '../optimizations/omega/af2/w20/';
%                 '../optimizations/omega/af2/w25/'};
% lgds = {'$\omega=1.0$','$\omega=0.9$','$\omega=1.5$','$\omega=2.0$','$\omega=2.5$'};
% linestyles = {'-','-','--','-','--'};
% xlims = [-.1,1.1]*1e4;
% ylims = [5e-7,.1];
% filename = '../../../P02/figuras/testes/w_optimization_af2';
% colors = {[0,0,0],cmap(210,:),cmap(210,:),cmap(100,:),cmap(100,:)};

% a - adi
address_list = {'../optimizations/base/adi/';
                '../optimizations/alpha/adi/a15/';
                '../optimizations/alpha/adi/a1e-1/';
                '../optimizations/alpha/adi/a1e-2/';
                '../optimizations/alpha/adi/a1e-3/';
                '../optimizations/alpha/adi/a1e-4/';
                '../optimizations/alpha/adi/a1e-5/';
                '../optimizations/alpha/adi/aseq_suggestion/';
                '../optimizations/alpha/adi/aseq_opt/'};
lgds = {'$\alpha=1.0$','$\alpha=1.5$','$\alpha=1\mathrm{E}-01$','$\alpha=1\mathrm{E}-02$','$\alpha=1\mathrm{E}-03$','$\alpha=1\mathrm{E}-04$','$\alpha=1\mathrm{E}-05$','Seq. sugerida','$\alpha_L=1\mathrm{E}-04,\alpha_H=1\mathrm{E}-01$'};
linestyles = {'-','-','-.','--','-','-.','--',':','--'};
xlims = [-.05,1]*1e4;
ylims = [5e-7,.1];
filename = '../../../P02/figuras/testes/a_optimization_adi';
colors = {[0,0,0],cmap(210,:),cmap(210,:),cmap(210,:),cmap(100,:),cmap(100,:),cmap(100,:),'k','k'};

% % a - af2
% address_list = {'../optimizations/base/af2/';
%                 '../optimizations/alpha/af2/a01/';
%                 '../optimizations/alpha/af2/a20/';
%                 '../optimizations/alpha/af2/a50/';
%                 '../optimizations/alpha/af2/a1e2/';
%                 '../optimizations/alpha/af2/aseq_suggestion/';
%                 '../optimizations/alpha/af2/aseq_opt/'};
% lgds = {'$\alpha=1.0$','$\alpha=0.1$','$\alpha=20$','$\alpha=50$','$\alpha=1\mathrm{E}02$','Seq. sugerida','$\alpha_L=1,\alpha_H=1\mathrm{E}04$'};
% linestyles = {'-','-','-.','--','-','-.','--',':','--'};
% % xlims = [-.05,1]*1e4;
% % ylims = [5e-7,.1];
% % filename = '../../../P02/figuras/testes/a_optimization_af2';
% xlims = [8100,8150];
% ylims = [1.04,1.18]*1e-6;
% filename = '../../../P02/figuras/testes/a_optimization_af2_detail';
% colors = {[0,0,0],cmap(210,:),cmap(210,:),cmap(100,:),cmap(100,:),cmap(100,:),'k'};

% % r - slor
% address_list = {'../optimizations/base/slor/';
%                 '../optimizations/r/slor/r09/';
%                 '../optimizations/r/slor/r15/';
%                 '../optimizations/r/slor/r20/';
%                 '../optimizations/r/slor/r25/';
%                 '../optimizations/r/slor/r30/';
%                 '../optimizations/r/slor/r40/'};
% lgds = {'$r=1.0$','$r=0.9$','$r=1.5$','$r=2.0$','$r=2.5$','$r=3.0$','$r=4.0$'};
% linestyles = {'-','-','-.','--','-','-.','--'};
% xlims = [-.2,5]*1e4;
% ylims = [5e-7,.1];
% filename = '../../../P02/figuras/testes/r_optimization_slor';
% colors = {[0,0,0],cmap(210,:),cmap(210,:),cmap(210,:),cmap(100,:),cmap(100,:),cmap(100,:)};

% % optimized - slor
% address_list = {'../optimizations/base/slor/';
%                 '../optimizations/optimized/slor/'};
% lgds = {'Base','Otimizado'};
% linestyles = {'-','-'};
% xlims = [-.1,3.05]*1e4;
% ylims = [5e-7,.1];
% filename = '../../../P02/figuras/testes/optimized_slor';
% colors = {[0,0,0],cmap(210,:)};

% % optimized - adi
% address_list = {'../optimizations/base/adi/';
%                 '../optimizations/optimized/adi/'};
% lgds = {'Base','Otimizado'};
% linestyles = {'-','-'};
% xlims = [-.1,2.5]*1e4;
% ylims = [5e-7,.1];
% filename = '../../../P02/figuras/testes/optimized_adi';
% colors = {[0,0,0],cmap(210,:)};

% % optimized - af2
% address_list = {'../optimizations/base/af2/';
%                 '../optimizations/optimized/af2/'};
% lgds = {'Base','Otimizado'};
% linestyles = {'-','-'};
% xlims = [-.05,1]*1e4;
% ylims = [5e-7,.1];
% filename = '../../../P02/figuras/testes/optimized_af2';
% colors = {[0,0,0],cmap(210,:)};

% % ecm cases
% address_list = {'../final_results/ecm/';
%                 '../final_results/ecm_control/'};
% lgds = {'Sem controle','Com controle'};
% linestyles = {'-','-'};
% xlims = [-100,1500];
% ylims = [1e-8,1];
% filename = '../../../P02/figuras/resultados/ecm_res';
% colors = {cmap(210,:),cmap(100,:)};

fontsize = 14;
lnwdth = 1.5;
pos = [500,500];
pba = [1,1];


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure(1),clf
f_res_log('x',address_list,casename,colors,linestyles,lgds,ref_val,fontsize,lnwdth,pos,pba,xlims,ylims,save_figs,[filename,'_x.eps'])
figure(2),clf
f_res_log('y',address_list,casename,colors,linestyles,lgds,ref_val,fontsize,lnwdth,pos,pba,xlims,ylims,save_figs,[filename,'_y.eps'])
