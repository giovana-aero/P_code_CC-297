clc,clear

clc,clear

address = '../results/';
casename = 'fullp';

mesh.x = dlmread([address,casename,'_x_mesh.dat']);
mesh.y = dlmread([address,casename,'_y_mesh.dat']);
u = dlmread([address,casename,'_u.dat']);
v = dlmread([address,casename,'_v.dat']);

% mid = round(size(u,2)/2);

% [x_up,x_down] = rearrange(mesh.x);
% [y_up,y_down] = rearrange(mesh.y);
% x = rearrange(mesh.x);
% y = rearrange(mesh.y);
x = mesh.x;
y = mesh.y;


% [startx,starty] = meshgrid(-1,0);
% [x,y] = meshgrid(0:0.1:1,0:0.1:1);
% [X,Y] = meshgrid(1:92,1:15);

figure(1),clf,hold on
% scatter(x,y,'k')
plot(x(1,:),y(1,:))
% streamline(x,y,u,v)
% scatter(x_up,y_up,'k')
% scatter(x_down,y_down,'k')
quiver(mesh.x(:,1:end-1),mesh.y(:,1:end-1),u(:,1:end-1),v(:,1:end-1),.1);
% streamline(mesh.x,mesh.y,u,v,-2,1);
% streamline(mesh.x(:,1:end-1),mesh.y(:,1:end-1),u(:,1:end-1),v(:,1:end-1),-2,1);
% streamline(X,Y,u(:,1:end-1),v(:,1:end-1),1,1);
grid on,axis equal


function A = rearrange(A)

mid = round(size(A,2)/2);

half_up = A(:,mid:end);
half_down = A(:,1:mid);

half_down = flip(half_down,2);
half_down = flip(half_down,1);

A = [half_up;half_down];

end

% ideia: implementar manualmente streamlines
% - em um dado ponto (x,y), checar quais são as coordenadas mais próximas na
%   malha
% - obter u e v a partir de uma média ponderada das velocidades na malha
% - encontrar um novo cur_x e cur_y por meio de um euler explícito
% https://www.mathworks.com/matlabcentral/answers/521369-plotting-streamlines-in-3d-from-a-non-rectangular-grid
% function f_streamline(x,y,u,v,startx,starty,euler_prmtrs)

% euler_prmtrs.dt
% euler_prmtrs.n (número total de passos)

% cur_x = startx;
% cur_y = starty;

% dx = min(min(abs(x-cur_x)));
% dy = min(min(abs(y-cur_y)));
% cur_u = 

% end