clc,clear

address = '../reverse_cst/';

addpath(address);

% filename = 'bi_air_10.txt';
% filename = 'naca0010.txt';
% filename = 'naca2412.txt';
% filename = 'whitcomb.txt';
% filename = 'fx61-163.txt';
% filename = 'naca64A005.92.txt';
% filename = 's1223.txt';
filename = 'supercritical_TM-X-1831.txt';

fontsize = 12;
base_size = 500;

save_figs = 1;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

xlims = [-.1,1.1];
ylims = [-1,1]*.2;

coo = dlmread([address,filename]);
coo = converter_function_xfoil(coo);

figure(1),clf
plot(coo(:,1),coo(:,2),'k','linewidth',1.5),grid on,axis equal
xlim(xlims)
ylim(ylims)
xylabel_latex('x','y')

fig_height = base_size*diff(ylims)/diff(xlims);
% fig_height = base_size;
fig_width = base_size;
pba = [fig_width,fig_height]/fig_width;
pos = [fig_width,fig_height];

set_rndrr_painters();
set_fontsize_position(fontsize,pos,pba)

if save_figs
  filename = [filename(1:end-3),'eps'];
  save_fig_eps(filename);
end
