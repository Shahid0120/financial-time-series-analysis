% (a)
Ssol = 100;
Ksol = 90;
rsol = 0.05;
Tsol = 1;
csol = 30;

% (b)
gsol = @(x) blackscholes(Ssol,Ksol,rsol,x, Tsol)-csol;

% (c)
xdsol = linspace(0,1,21);
xdsol = xdsol';

% (d)
ydsol = gsol(xdsol);

% (e)
Bsol = [ones(21,1), xdsol, xdsol.^2, xdsol.^3, xdsol.^4, xdsol.^5];
alssol = Bsol \ ydsol;
% f)
ppsol = spline(xdsol, ydsol);
% g)
xplsol = linspace(0,1,101);
% h)
ylssol = alssol(1) + alssol(2)*xplsol + alssol(3)*xplsol.^2 + alssol(4)*xplsol.^3 + alssol(5)*xplsol.^4 + alssol(6)*xplsol.^5;
yspsol = ppval(ppsol, xplsol);
% i)
%plot(xdsol,ydsol, '*', xplsol,ylssol, xplsol, yspsol);
%legend('data', 'least squares', 'spline')

% Solutions 
S = 100 ;
K = 90;
r = 0.05;
T = 1 ;
c = 30;

xd = [0; 0.05; 0.1 ; 0.15; 0.2 ; 0.25; 0.3 ; 0.35; 0.4; 0.45;0.5 ;0.55;0.6 ; 0.65;0.7 ; 0.75; 0.8;  0.85; 0.9 ; 0.95];
yd = [-15.61064821 ;-15.60944404 ;-15.37116238 ;-14.53284094 ;-13.30055159;-11.85923705 ;-10.30255791  ;-8.67888064 ; -7.01521094  ;-5.32783836; -3.62728512 ; -1.92076745  ;-0.21350029  ; 1.49057152  ; 3.18834342;4.87728811  ; 6.55528873 ;  8.22052999 ;  9.87142519 ; 11.50656615];
als = [-15.56400512 ; -10.44213353 ; 152.0882548 ; -253.24762122 ; 205.85135869; -65.67768568];

marks = 0;

if exist('Ssol', 'var') && exist('Ksol', 'var') == 1
    if Ssol - S == 0
        marks = marks + 1;
        fprintf('q3)a) Ssol, Ksol,... exists and defined correctly\n');
    else
        fprintf('q3)a) Ssol not defined correctly\n');
    end
else
    fprintf('q3)a) Ssol, Ksol,... not defined\n');
end

fprintf('q3)b) Is the function g defined correctly? (+1 mark)\n');

if exist('xd', 'var') == 1
        if size(xd) == size(xdsol)
            if norm(xd-xdsol) < 10^(-8)
                fprintf('q3)c) xd exists and is correct')
                marks = marks + 1;
            else
                fprintf('q3)c) xd exists but is not correct')
                marks = marks + 0;
            end
        else
            fprintf('q3)c) xd is not of the correct size')
        end
end
       

if exist('yd', 'var') == 1
        if size(yd) == size(ydsol)
            if norm(yd-ydsol) < 10^(-8)
                fprintf('q3)d) yd exists and is correct')
                marks = marks + 1;
            else
                fprintf('q3)d) yd exists but is not correct')
                marks = marks + 0;
            end
        else
            fprintf('q3)d) yd is not of the correct size')
        end
end
        

if exist('als', 'var') == 1
    if prod(size(als) == size(alssol))
        fprintf('norm(als-alssol) = %e\n',norm(als-alssol,Inf));
        if norm(als-alssol,Inf) < 1e-10
          marks = marks + 2;
          fprintf('q3)e) als is within the specified tolerance\n');
        else
          fprintf('q3)e) als is not within the specified tolerance\n');
        end;
    else
        fprintf('q3)e) als is not of correct size\n');
    end;
else
    fprintf('q3)e) als is not defined\n');
end;
fprintf('\n als(1) = %1.4f, als(2) = %1.4f, als(3) = %1.4f \n\n', alssol(1), alssol(2), alssol(3));


fprintf('q3)f) Is the cubic spline computed correctly? (+ 1 mark)\n')


if exist('xpl', 'var') == 1
    if prod(size(xpl) == size(xpl))
        fprintf('norm(xpl-xplsol) = %e\n',norm(xpl-xplsol,Inf));
        if norm(xpl-xplsol,Inf) < 1e-10
          marks = marks + 1;
          fprintf('q3)g) xpl is within the specified tolerance\n');
        else
          fprintf('q3)g) xpl is not within the specified tolerance\n');
        end;
    else
        fprintf('q3)g) xpl is not of correct size\n');
    end;
else
    fprintf('q3)g) xpl is not defined\n');
end;



if exist('yls', 'var') == 1
    if prod(size(yls) == size(ylssol))
        fprintf('norm(yls-ylssol) = %e\n',norm(yls-ylssol,Inf));
        if norm(yls-ylssol,Inf) < 1e-10
          marks = marks + 1;
          fprintf('q3)h) yls is within the specified tolerance\n');
        else
          fprintf('q3)h) yls is not within the specified tolerance\n');
        end;
    else
        fprintf('q3)h) yls is not of correct size\n');
    end;
else
    fprintf('q3)h) yls is not defined\n');
end;


if exist('ysp', 'var') == 1
    if prod(size(ysp) == size(yspsol))
        fprintf('norm(ysp-yspsol) = %e\n',norm(ysp-yspsol,Inf));
        if norm(ysp-yspsol,Inf) < 1e-10
          marks = marks + 1;
          fprintf('q3)h) ysp is within the specified tolerance\n');
        else
          fprintf('q3)h) ysp is not within the specified tolerance\n');
        end;
    else
        fprintf('q3)h) ysp is not of correct size\n');
    end;
else
    fprintf('q3)h) ysp is not defined\n');
end;

% fprintf('q3)g) Is the plot correct? (+ 1 mark)\n')


fprintf('\n\n');

fprintf('q3.m suggested marks: %d out of 10 marks.\n', marks);
fprintf('The final marks will depend on how you obtained the solution\n');
fprintf('and may vary from the suggested marks.\n');


marks3 = marks;
