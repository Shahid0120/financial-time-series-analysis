
% (a)
fsol = @(x) exp(-x).*sin(10*pi*x);
% (b)
IQsol = quad(fsol, 0, 10);
% (c)
% (i)
NIsol = 11;
% (ii)
xIsol = linspace(0,10,NIsol);
xIsol = xIsol(:);
% (iii)
wIsol = ones(NIsol,1);
wIsol(1) = 1/2; wIsol(end) = 1/2;
wIsol = 10*wIsol/(NIsol-1);
ITsol = fsol(xIsol)'*wIsol;
% (d)
fprintf('|IQsol-ITsol| = %e \n', abs(IQsol-ITsol));
% (e)
N2sol = 555;  % This is the smallest number such that |IQ-IT2| < 10^(-3)
x2sol = linspace(0,10,N2sol);
% (iii)
w2sol = ones(N2sol,1);
w2sol(1) = 1/2; w2sol(end) = 1/2;
w2sol = 10*w2sol/(N2sol-1);
IT2sol = fsol(x2sol)*w2sol;
abs(IT2sol-IQsol)

