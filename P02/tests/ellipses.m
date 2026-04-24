clc,clear

% a = 1;
% b = .5;

% phi = linspace(0,2*pi,101);
% phi = phi(1:end-1);

% x = a*cos(phi);
% y = b*sin(phi);

% figure(1),clf
% plot(x,y)
% grid on,axis equal

function x = cosspace_half(startP,endP,N)
    
  % Versão modificada da função cosspace. Gera apenas metade da distribuição.
  
  
  x = zeros(1,N); x(N) = endP;
  angleInc = pi/(N-1);
      
  curAngle=angleInc;
  for i = 2:N-1
      x(i)=endP*(1-cos(curAngle));
      curAngle=curAngle+angleInc;
  end
       
end

% function s=lnspc(xi,xf,n)

% s = zeros(1,n);
% s(1) = xi;
% step = (xf-xi)/(n-1);

% for i = 2:n
%   s(i) = s(i-1) + step;
% end

% end

% disp(linspace(1,5,5))
% disp(lnspc(1,5,5))

n = 50;
x = cosspace_half(0,1,n);
figure(1),clf
plot(x,zeros(1,n),'o'),grid on