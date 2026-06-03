clc,clear

syms x0 x1 x2
syms y0 y1 y2

A = [x0^2,x0,1;
     x1^2,x1,1;
     x2^2,x2,1];
b = [y0;y1;y2];

disp(A\b)