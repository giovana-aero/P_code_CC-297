clc,clear

x = linspace(0,1,2000);

% % a generic function
%f = @(x) x.^5;
%dx_f = @(x) 5*x.^4;
%dx2_f = @(x) 4*5*x.^3;
%dx3_f = @(x) 3*4*5*x.^2;
%dx4_f = @(x) 2*3*4*5*x;
%dx5_f = 2*3*4*5;

%y2 = f(a) + dx_f(a)*(x - a) + dx2_f(a)*(x - a).^2/factorial(2) + dx3_f(a)*(x - a).^3/factorial(3); + dx4_f(a)*(x - a).^4/factorial(4) + dx5_f*(x - a).^5/factorial(5);

%y1 = f(x);

a = .01;

% naca4
t = 0.12;
a0 = 0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 = 0.2843;
a4 = -0.1015;

f = @(x) 5*t*(a0*x.^0.5 + a1*x + a2*x.^2 + a3*x.^3 + a4*x.^4);
df = @(x) 5*t*(0.5*a0*x.^(0.5-1) + 2*a2*x + 3*a3*x.^2 + 4*a4*x.^3);
d2f = @(x) 5*t*(-0.25*a0*x.^(0.5-2) + 6*a3*x + 12*a4*x.^2);
d3f = @(x) 5*t*(0.3750*a0*x.^(0.5-2) + 24*a4*x);
d4f = @(x) 5*t*(-0.5625*a0*x.^(0.5-3));
d5f = @(x) 5*t*(1.40625*a0*x.^(0.5-4));
d6f = @(x) 5*t*(-4.921875*a0*x.^(0.5-5));
d7f = @(x) 5*t*(22.1484375*a0*x.^(0.5-6));

y1 = f(x);

y2 = f(a) + df(a)*(x - a) + d2f(a)*(x - a).^2/factorial(2) + d3f(a)*(x - a).^3/factorial(3); + d4f(a)*(x - a).^4/factorial(4) + d5f(a)*(x - a).^5/factorial(5) + d6f(a)*(x - a).^6/factorial(6) + d7f(a)*(x - a).^7/factorial(7);

figure(1),clf
plot(x,y1,'k'),hold on,grid on,axis equal
plot(x,y2,'m--')
xlim([0,1])
ylim([-.5,.5])


y1 = df(x);
y2 = df(a) + d2f(a)*(x - a) + d3f(a)*(x - a).^2/factorial(2) + d4f(a)*(x - a).^3/factorial(3); + d5f(a)*(x - a).^4/factorial(4) + d6f(a)*(x - a).^5/factorial(5) + d7f(a)*(x - a).^6/factorial(6);

figure(2),clf
plot(x,y1,'k'),hold on,grid on,axis equal
plot(x,y2,'m--')
xlim([0,1])
ylim([-.5,.5])


% conclusion: none of interest except that I have no idea how to deal with that
% derivative at x=0

