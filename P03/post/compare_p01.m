clc,clear

casename = 'fullp';

cmap = slanCM('bupu');

% fullp
address = '../mach_sweep/biair_M010/';

% biair
address_p01 = '../../P01/final_results_slor/t010_mtype4/bi_air_cp_chord.dat';

% biair, reference
address_ref = '../../../P01/ref_data/t010.dat';

lgds = {'Referencia (P01)','Pequenas perturbacoes (P01)','Potencial completo'};
linestyles = {'-','--','--'};
colors = {[0,0,0],cmap(100,:),cmap(210,:)};
filename = '../../../P03/figuras/testes/comparison_p01_cp_airfoil';

save_fig = 1;
fontsize = 14;
lnwdth = 1.5;
pos = [1000,500];
pba = [2,1];

figure(1),clf,grid on,hold on
cp_ref = dlmread(address_ref);
xc_ref = cp_ref(:,1);
cp_ref = cp_ref(:,3);
plot(xc_ref,cp_ref,'color',colors{1},'linestyle',linestyles{1},'displayname',lgds{1},'linewidth',1.5)

cp_p01 = dlmread(address_p01);
xc_p01 = linspace(0,1,length(cp_p01));
plot(xc_p01,cp_p01,'color',colors{2},'linestyle',linestyles{2},'displayname',lgds{2},'linewidth',1.5)

mesh_x = readmatrix([address,casename,'_x_mesh.dat'],'filetype','delimitedtext');
cp = readmatrix([address,casename,'_cp.dat'],'filetype','delimitedtext');
  
chord_nump = (size(cp,2) + 1)/2;

plot(mesh_x(1,1:chord_nump),cp(1,1:chord_nump),'color',colors{3},'linestyle',linestyles{3},'displayname',lgds{3},'linewidth',1.5)


xylabel_latex('x/c','C_p')
lgd_latex(legend())
inv_yaxis();


xlim([-.1,1.1])
ylim([-.5,.6])

set_fontsize_position(fontsize,pos,pba);

if save_fig
  set_rndrr_painters()
  save_fig_eps(filename)
end