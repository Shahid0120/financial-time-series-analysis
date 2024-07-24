n1 = 11;
rnorm = 1.305622276959184e-13;
xsol = [18.97814248; 15.81188301; 23.05101644; 18.41194001; 19.72457047; 14.73766308; 13.82359504; 9.12508583;  8.22294476 ; 4.05003433 ; 3.4891779 ];
rsol = [ 4.26325641e-14 ; 3.73034936e-14 ;  8.88178420e-15 ;  3.55271368e-15 ;  8.88178420e-16  ; -3.55271368e-15 ;  -2.13162821e-14 ; 4.44089210e-15 ; -3.10862447e-15 ; 8.88178420e-16 ; 3.99680289e-15];
% (a)
n1sol = 11;
% (b)
A1sol = diag(5*ones(11,1),0) + diag(-2*ones(9,1), -2) + diag(-3*ones(9,1),2) + diag(-ones(6,1),5);
% (c)
[Q1sol, R1sol] = qr(A1sol);
% (d)
Q1checksol = norm(Q1sol'*Q1sol - eye(n1sol),2);
% (e)
bsol = [n1sol:-1:1];
bsol=bsol(:);
% (f)
xsolsol = R1sol \ (Q1sol'*bsol);
% (g)
rsolsol = A1sol*xsolsol - bsol;
% (h)
rnormsol = norm(rsolsol,1);
% Get the size of rsol
rsol_size = size(rsol);
rsolsol_size = size(rsolsol);
% Print the size of rsol
fprintf('The size of rsol is: %d by %d\n', rsol_size(1), rsol_size(2));
fprintf('The size of rsolsol is: %d by %d\n', rsolsol_size(1), rsolsol_size(2));

fprintf('size', size(rsolsol))

if exist('rsol', 'var') == 1
    if prod(size(rsol) == size(rsolsol))
        fprintf('norm(rsol-rsolsol,inf) = %e\n',norm(rsol-rsolsol,Inf));
        if norm(rsol-rsolsol,Inf) < 1e-4
          fprintf('q1)g) rsol is within the specified tolerance\n\n');
        else
          fprintf('q1)g) rsol is not within the specified tolerance\n\n');
        end;
    else
        fprintf('q1)g) rsol is not of correct size\n\n');
    end;
else
    fprintf('q1)g) rsol is not defined\n\n');
end;

if exist('rnorm', 'var') == 1
    if prod(size(rnorm) == size(rnormsol))
        fprintf('norm(rnorm-rnormsol,inf) = %e\n',norm(rnorm-rnormsol,Inf));
        if norm(rnorm-rnormsol,Inf) < 1e-4
          fprintf('q1)h) rnorm is within the specified tolerance\n\n');
        else
          fprintf('q1)h) rnorm is not within the specified tolerance\n\n');
        end;
    else
        fprintf('q1)h) rnorm is not of correct size\n\n');
    end;
else
    fprintf('q1)h) rnorm is not defined\n\n');
end;
if exist('xsol', 'var') == 1
    if prod(size(xsol) == size(xsolsol))
        fprintf('norm(xsol-xsolsol,inf) = %e\n',norm(xsol-xsolsol,Inf));
        if norm(xsol-xsolsol,Inf) < 1e-4
          fprintf('q1)f) xsol is within the specified tolerance\n\n');
        else
          fprintf('q1)f) xsol is not within the specified tolerance\n\n');
        end;
    else
        fprintf('q1)f) xsol is not of correct size\n\n');
    end;
else
    fprintf('q1)f) xsol is not defined\n\n');
end;