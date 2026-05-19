clc,clear


casename = 'eom';
ref_val = 1e-6;
save_figs = 0;
cmap = slanCM('bupu');


address_list = {'../results/'};
lgds = {''};
linestyles = {'-'};
xlims = [];
ylims = [];
filename = 'r_optimization_sor';

colors = config_colors_cmap(length(address_list),cmap,1);

fontsize = 14;
lnwdth = 1.5;
pos = [1000,500];
pba = [2,1];

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure(1)
f_res_log('x',address_list,casename,colors,linestyles,lgds,ref_val,fontsize,lnwdth,pos,pba,xlims,ylims,save_figs)
figure(2)
f_res_log('y',address_list,casename,colors,linestyles,lgds,ref_val,fontsize,lnwdth,pos,pba,xlims,ylims,save_figs)
