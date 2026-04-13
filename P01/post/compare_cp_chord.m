clc,clear

casename = 'bi_air';
cmap = slanCM('bupu');
save_figs = 1;

% bi_air, t005
address_list = {'../final_results_slor/t005_mtype1/',...
                '../final_results_slor/t005_mtype4/'};
address_ref = {'../original_cp_chord_005.dat',...
               '../../../P01/ref_data/t005.dat'};
lgds = {'Malha 1','Malha 4'};
lgds_ref = {'Referencia (corda)','Referencia (geometria completa)'};
linestyles = {'-','-'};
linestyles_ref = {'--',':'};
colors = {cmap(100,:),cmap(150,:)};
colors_ref = {'k','k'};
ILE = [11,81];
ITE = [31,241];
filename = 'compare_cp_chord_t005';

% % bi_air, t010
% address_list = {'../final_results_slor/t010_mtype1/',...
%                 '../final_results_slor/t010_mtype4/'};
% address_ref = {'../../../P01/ref_data/t010.dat'};
% lgds = {'Malha 1','Malha 4'};
% lgds_ref = {'Referencia (geometria completa)'};
% linestyles = {'-','-'};
% linestyles_ref = {'--',':'};
% colors = {cmap(100,:),cmap(150,:)};
% colors_ref = {'k','k'};
% ILE = [11,81];
% ITE = [31,241];
% filename = 'compare_cp_chord_t010';

% % naca0005
% address_list = {'../final_results_slor/naca0005_mtype1/',...
%                 '../final_results_slor/naca0005_mtype3/'};
% address_ref = {'../naca0005_cp_chord.dat'};
% lgds = {'Malha 1','Malha 3'};
% lgds_ref = {'XFOIL'};
% linestyles = {'-','-'};
% linestyles_ref = {'--'};
% colors = {cmap(100,:),cmap(200,:)};
% colors_ref = {'k'};
% ILE = [11,41];
% ITE = [31,121];
% filename = 'compare_cp_chord_naca0005';

% % naca0010
% address_list = {'../final_results_slor/naca0010_mtype1/',...
%                 '../final_results_slor/naca0010_mtype3/'};
% address_ref = {'../naca0010_cp_chord.dat'};
% lgds = {'Malha 1','Malha 3'};
% lgds_ref = {'XFOIL'};
% linestyles = {'-','-'};
% linestyles_ref = {'--'};
% colors = {cmap(100,:),cmap(200,:)};
% colors_ref = {'k'};
% ILE = [11,41];
% ITE = [31,121];
% filename = 'compare_cp_chord_naca0010';

fontsize = 14;
lnwdth = 1.5;
pos = [1000,500];
pba = [2,1];

figure(1),clf,grid on,hold on
for i = 1:length(address_list)
  mesh_x = readmatrix([address_list{i},casename,'_mesh_x.msh'],'filetype','delimitedtext');
  cp = readmatrix([address_list{i},casename,'_cp_chord.dat'],'filetype','delimitedtext');
  chord = mesh_x(ILE(i):ITE(i));

  plot(chord,cp,'color',colors{i},'linestyle',linestyles{i},'displayname',lgds{i},'linewidth',lnwdth)

  % plot(0:(length(res)-1),res,'color',colors{i},'linestyle',linestyles{i},...
  %         'displayname',lgds{i},'linewidth',lnwdth)
end

for i = 1:length(address_ref)
  cp_ref = readmatrix(address_ref{i},'filetype','delimitedtext');

  if contains(address_ref{i},'naca')
    plot(cp_ref(:,1),cp_ref(:,2),'color',colors_ref{i},'linestyle',linestyles_ref{i},'displayname',lgds_ref{i},'linewidth',lnwdth)
  elseif contains(address_ref{i},'/P01/ref_data/')
    plot(cp_ref(:,1),cp_ref(:,3),'color',colors_ref{i},'linestyle',linestyles_ref{i},'displayname',lgds_ref{i},'linewidth',lnwdth)
  else
    plot(cp_ref(:,1),-cp_ref(:,2),'color',colors_ref{i},'linestyle',linestyles_ref{i},'displayname',lgds_ref{i},'linewidth',lnwdth)
  end

end

xlim([-.1,1.1])
if contains(address_list{i},'t005')
  ylim([-.2,.3])
elseif contains(address_list{i},'t010')
  ylim([-.35,.5])
elseif contains(address_ref{1},'naca0005')
  ylim([-.3,.3])
elseif contains(address_ref{1},'naca0010')
  ylim([-.7,.5])
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

