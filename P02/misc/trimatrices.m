clc,clear

disp('')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

syms x1 x2 x3 x4 x5
syms w1 w2 w3 w4 w5 w6
syms k1 k2 k3 k4 k5

s1 = sym(1);
s0 = sym(0);

U = [s1,x1,s0,s0,s0,s0;
     s0,s1,x2,s0,s0,s0;
     s0,s0,s1,x3,s0,s0;
     s0,s0,s0,s1,x4,s0;
     s0,s0,s0,s0,s1,x5;
     s0,s0,s0,s0,s0,s1];
disp(U)
disp('')

L = [w1,s0,s0,s0,s0,s0;
     k1,w2,s0,s0,s0,s0;
     s0,k2,w3,s0,s0,s0;
     s0,s0,k3,w4,s0,s0;
     s0,s0,s0,k4,w5,s0;
     s0,s0,s0,s0,k5,w6];
disp(L)
disp('')

%M = U*L;
M = L*U;
disp(M)
disp('')


syms a1 a2 a3 a4 a5
syms b1 b2 b3 b4 b5 b6
syms c1 c2 c3 c4 c5

A = [b1,c1,s0,s0,s0,s0;
     a1,b2,c2,s0,s0,s0;
     s0,a2,b3,c3,s0,s0;
     s0,s0,a3,b4,c4,s0;
     s0,s0,s0,a4,b5,c5;
     s0,s0,s0,s0,a5,b6];
disp(A)
disp('')


syms y1 y2 y3 y4 y5 y6

y = [y1;y2;y3;y4;y5;y6];
disp(L*y)
disp('')


syms u1 u2 u3 u4 u5 u6

f = [u1;u2;u3;u4;u5;u6];
disp(U*f)
