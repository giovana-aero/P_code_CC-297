function f_res_log(xy,address_list,casename,colors,linestyles,lgds,ref_val,fontsize,lnwdth,pos,pba,xlims,ylims,save_figs,filename)

% figure(1),clf,grid on,hold on
hold on,grid on
for i = 1:length(address_list)
  % res =readmatrix([address_list{i},casename,'_',xy,'.log'],'filetype','delimitedtext');
  res = dlmread([address_list{i},casename,'_',xy,'.log']);
  disp([address_list{i},casename,'_',xy,'.log'])

  plot(0:(length(res)-1),res,'color',colors{i},'linestyle',linestyles{i},...
          'displayname',lgds{i},'linewidth',lnwdth)

  fprintf("%s: iter %d\n",lgds{i},find(res<ref_val,1));

end

set(gca,'yscale','log')

xylabel_latex('n',['|L',xy,'|_{max}']);
set_fontsize_position(fontsize,pos,pba);
% lgd = legend();
lgd_latex(legend())
if ~isempty(xlims)
  xlim(xlims)
end
if ~isempty(ylims)
  ylim(ylims)
end

set_rndrr_painters();

if save_figs
  save_fig_eps(filename);
end

end