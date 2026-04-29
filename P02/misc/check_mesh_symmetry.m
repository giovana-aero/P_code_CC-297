clc,clear

address = '../results/';
casename = 'eom';

iter = 10;

if iter >= 0
  printf('%010d\n',iter);
  result = {sprintf("%s_x_iter_%010d.dat",casename,iter),...
            sprintf("%s_y_iter_%010d.dat",casename,iter)};
else
  disp('initial condition');
  result = {[casename,'_x_initial.dat'],...
            [casename,'_y_initial.dat']};
end

data_x = dlmread([address,result{1}]);
data_y = dlmread([address,result{2}]);

l_ksi = (size(data_x,2)-1)/2 + 1;
disp(isequal(data_x(:,1:l_ksi),flip(data_x(:,l_ksi:end),2)))
disp(isequal(data_y(:,1:l_ksi),-flip(data_y(:,l_ksi:end),2)))
