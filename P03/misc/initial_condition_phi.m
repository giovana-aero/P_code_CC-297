clc,clear

syms A1 A2 A3 dphi_dksi dphi_deta U V

A = [A1,A2;
     A2,A3];
b = [U;V];

sol = A\b;
disp(sol);