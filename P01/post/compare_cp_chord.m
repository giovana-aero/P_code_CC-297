clc,clear

casename = 'bi_air';
cmap = slanCM('bupu');
save_figs = 1;

% bi_air, t005
address_list = {'../final_results_slor/t005_mtype1/',...
                '../final_results_slor/t005_mtype3/'};
address_ref = '../original_cp_chord_005.dat';
filename = 'compare_cp_chord_t005';
lgds = {'Malha 1','Malha 3'};
linestyles = {'-','-'};
colors = {cmap(100,:),cmap(200,:)};

% % naca0005
% address_list = {'../final_results_slor/naca0005_mtype1/',...
%                 '../final_results_slor/naca0005_mtype3/'};
% address_ref = '../naca0005_cp_chord.dat';
% filename = 'compare_cp_chord_naca0005';
% lgds = {'Malha 1','Malha 3'};
% linestyles = {'-','-'};
% colors = {cmap(100,:),cmap(200,:)};

% % naca0010
% address_list = {'../final_results_slor/naca0010_mtype1/',...
%                 '../final_results_slor/naca0010_mtype3/'};
% address_ref = '../naca0010_cp_chord.dat';
% filename = 'compare_cp_chord_naca0010';
% lgds = {'Malha 1','Malha 3'};
% linestyles = {'-','-'};
% colors = {cmap(100,:),cmap(200,:)};

% i = 3;
% ILE = [11,21,41,81,6];
% ITE = [31,61,121,241,16];
ILE = [11,41];
ITE = [31,121];

% lgds = {'Malha 1','Malha 2','Malha 3','Malha 4','Malha 5'};
% linestyles = {'-','--','-','-.','-'};

fontsize = 14;
lnwdth = 1.5;
pos = [1000,500];
pba = [2,1];

% colors = config_colors_cmap(length(address_list),cmap,0);

% figure(1),clf,grid on,hold on
% mesh_x = readmatrix([address,casename,'_mesh_x.msh'],'filetype','delimitedtext');
% cp = readmatrix([address,casename,'_cp_chord.dat'],'filetype','delimitedtext');
% chord = mesh_x(ILE(i):ITE(i));
% plot(chord,cp,'color',cmap(30,:),'linestyle','-','displayname','Simulação','linewidth',lnwdth)


figure(1),clf,grid on,hold on
for i = 1:length(address_list)
  mesh_x = readmatrix([address_list{i},casename,'_mesh_x.msh'],'filetype','delimitedtext');
  cp = readmatrix([address_list{i},casename,'_cp_chord.dat'],'filetype','delimitedtext');
  chord = mesh_x(ILE(i):ITE(i));

  plot(chord,cp,'color',colors{i},'linestyle',linestyles{i},'displayname',lgds{i},'linewidth',lnwdth)

  % plot(0:(length(res)-1),res,'color',colors{i},'linestyle',linestyles{i},...
  %         'displayname',lgds{i},'linewidth',lnwdth)
end

cp_ref = readmatrix(address_ref,'filetype','delimitedtext');
if contains(address_ref,'naca')
  plot(cp_ref(:,1),cp_ref(:,2),'k--','displayname','Referencia','linewidth',lnwdth)
else
  plot(cp_ref(:,1),-cp_ref(:,2),'k--','displayname','Referencia','linewidth',lnwdth)
end

if contains(address_ref,'naca0005')
  xlim([-.1,1.1])
  ylim([-.3,.3])
elseif contains(address_ref,'naca0010')
  xlim([-.1,1.1])
  ylim([-.7,.5])
else
  xlim([-.1,1.1])
  ylim([-.2,.3])
end

xylabel_latex('x/c','Cp');
set_fontsize_position(fontsize,pos,pba);
inv_yaxis();
lgd = legend();
lgd_latex(lgd)

set_rndrr_painters();

if save_figs
  save_fig_eps(filename);
end

