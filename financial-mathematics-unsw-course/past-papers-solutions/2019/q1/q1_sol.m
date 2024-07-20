nsol = 10;

% (b)
Asol = diag(-2*ones(nsol,1), 0) + diag(ones(nsol-2,1), 2) + diag(3*ones(nsol-4,1), -4);

% (c)
[Lsol, Usol, psol] = lu(Asol, 'vector');

% (d)
bsol = [20, 18, 16, 14, 12, 10, 8, 6, 4, 2]';

% (e)
ysol = Lsol\bsol(psol);
xssol = Usol\ysol;
disp(xssol)
% (f)
rsol = Asol*xssol - bsol;
disp(rsol)
% (g)
r1normsol = norm(rsol,1);

disp(r1normsol)
