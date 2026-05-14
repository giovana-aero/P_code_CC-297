clc,clear

s = cosspace(0,1,20);
s = s*10;

Ds = zeros(1,length(s));

for i = 2:length(Ds)-1
  Ds(i) = (s(i) - s(i-1))/(s(end) - s(i-1));
end

disp(Ds);

figure(1),clf
plot(Ds,Ds,'o')