% (b)
% (i)
Nsol = 100;
% (ii)
a0sol = [1:(N-1)];
a1sol = [1:(N-2)];
a_1sol = [2:N-1];
Bsol = diag(1-a0sol-2*a0sol.^2,0) + diag(a_1sol.^2+a_1sol,-1) + diag(a1sol.^2,1);
xsol = [(1:(Nsol-1))/Nsol]';
bsol = [exp(xsol).*(xsol.^2 + xsol.^3)];
bsol(1) = bsol(1)-2*0;
bsol(end) = bsol(end) - (Nsol-1)^2*exp(1);
% (iii)
Usol = Bsol\bsol;
% (iv)
xasol = [(1:Nsol-1)/Nsol]';
Uerrsol = norm(Usol-xasol.*exp(xasol));

