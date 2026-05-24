function x = cosspace_half(startP,endP,N)
    
% Vers�o modificada da fun��o cosspace. Gera apenas metade da distribui��o.


x = zeros(1,N); x(N) = endP;
angleInc = pi/(N-1)/2;
    
curAngle=angleInc;
for i = 2:N-1
    x(i)=endP*(1-cos(curAngle));
    curAngle=curAngle+angleInc;
end
     
end