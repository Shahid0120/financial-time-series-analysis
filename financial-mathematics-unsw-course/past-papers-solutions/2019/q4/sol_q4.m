fsol = @(x) abs(cos(pi*x));

% (b)
%IQ = quad(f, 0, 2);
IQsol = integral(fsol,0,2);

% (c) 
% (i)
NIsol = 1001;
% (ii)
xIsol = linspace(0,2,NIsol);
% (iii)
hIsol = 2/(NIsol-1);
% (iv)
wIsol = ones(1,NIsol);
wIsol(1) = 1/2; 
wIsol(end) = 1/2;

ITsol = sum(fsol(xIsol).*wIsol)*hIsol;
