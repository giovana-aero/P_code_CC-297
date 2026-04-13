clc,clear

casename = 'bi_air';
cmap = slanCM('bupu');
save_figs = 1;


% % mesh convergence, pj
% address_list = {'../mesh_conv/pj/mtype1/',
%                 '../mesh_conv/pj/mtype2/',
%                 '../mesh_conv/pj/mtype3/',
%                 '../mesh_conv/pj/mtype4/',
%                 '../mesh_conv/pj/mtype5/'};
% lgds = {'Malha 1','Malha 2','Malha 3','Malha 4','Malha 5'};
% linestyles = {'-','--','-','-.','-'};
% filename = 'mesh_conv_pj';

% mesh convergence, gs
address_list = {'../mesh_conv/pj/mtype1/',
                '../mesh_conv/pj/mtype2/',
                '../mesh_conv/pj/mtype3/',
                '../mesh_conv/pj/mtype4/',
                '../mesh_conv/pj/mtype5/'};
lgds = {'Malha 1','Malha 2','Malha 3','Malha 4','Malha 5'};
linestyles = {'-','--','-','-.','-'};
filename = 'mesh_conv_gs';

% % mesh convergence, gs
% address_list = {'../mesh_conv/sor/mtype1/',
%                 '../mesh_conv/sor/mtype2/',
%                 '../mesh_conv/sor/mtype3/',
%                 '../mesh_conv/sor/mtype4/',
%                 '../mesh_conv/sor/mtype5/'};
% lgds = {'Malha 1','Malha 2','Malha 3','Malha 4','Malha 5'};
% linestyles = {'-','--','-','-.','-'};
% filename = 'mesh_conv_sor';

% % mesh convergence, gs
% address_list = {'../mesh_conv/lgs/mtype1/',
%                 '../mesh_conv/lgs/mtype2/',
%                 '../mesh_conv/lgs/mtype3/',
%                 '../mesh_conv/lgs/mtype4/',
%                 '../mesh_conv/lgs/mtype5/'};
% lgds = {'Malha 1','Malha 2','Malha 3','Malha 4','Malha 5'};
% linestyles = {'-','--','-','-.','-'};
% filename = 'mesh_conv_lgs';

% % mesh convergence, slor
% address_list = {'../mesh_conv/slor/mtype1/',
%                 '../mesh_conv/slor/mtype2/',
%                 '../mesh_conv/slor/mtype3/',
%                 '../mesh_conv/slor/mtype4/',
%                 '../mesh_conv/slor/mtype5/'};
% lgds = {'Malha 1','Malha 2','Malha 3','Malha 4','Malha 5'};
% linestyles = {'-','--','-','-.','-'};
% filename = 'mesh_conv_slor';


ILE = [11,21,41,81,6];
ITE = [31,61,121,241,16];

fontsize = 14;
lnwdth = 1.5;
pos = [1000,500];
pba = [2,1];

colors = config_colors_cmap(length(address_list),cmap,0);


figure(1),clf,grid on,hold on
for i = 1:length(address_list)
  mesh_x = readmatrix([address_list{i},casename,'_mesh_x.msh'],'filetype','delimitedtext');
  cp = readmatrix([address_list{i},casename,'_cp_chord.dat'],'filetype','delimitedtext');
  chord = mesh_x(ILE(i):ITE(i));

  plot(chord,cp,'color',colors{i},'linestyle',linestyles{i},'displayname',lgds{i},'linewidth',lnwdth)

  % plot(0:(length(res)-1),res,'color',colors{i},'linestyle',linestyles{i},...
  %         'displayname',lgds{i},'linewidth',lnwdth)
end

xlim([-.1,1.1])
ylim([-.2,.3])

xylabel_latex('x/c','Cp');
set_fontsize_position(fontsize,pos,pba);
inv_yaxis();
lgd = legend();
lgd_latex(lgd)

set_rndrr_painters();

if save_figs
  save_fig_eps(filename);
end


