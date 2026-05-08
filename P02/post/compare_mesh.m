clc,clear

address = {'../results/','../results_unmod/'};
casename = 'eom';

iter = -1;

result = cell(1,2);
for i = 1:2
  if iter < 0
    iter_list_x = dir([address{i},casename,'_x_iter_*']);
    iter_list_y = dir([address{i},casename,'_y_iter_*']);
    result{i} = {iter_list_x(end).name,iter_list_y(end).name};
  else
    result{i} = {sprintf("%s_x_iter_%010d.dat",casename,iter),...
              sprintf("%s_y_iter_%010d.dat",casename,iter)};
  end
end

figure(1),clf
f_plot_mesh(1,address{1},'',result{1})
f_plot_mesh(1,address{2},'',result{2},'--')
