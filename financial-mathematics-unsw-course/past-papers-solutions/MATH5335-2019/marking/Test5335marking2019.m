% Labtest 2018 Math5335 marking program


function Test5335marking

diary off
clear all
close all
format compact
format long e

if exist('Test5335.txt', 'file') == 2
    delete Test5335.txt
end;


diary on

fprintf('q1) ... out of 10\n\n');
fprintf('q2) ... out of 10\n\n');
fprintf('q3) ... out of 10\n\n');
fprintf('q4) ... out of 10\n\n');
fprintf('q5a)... out of 5 \n\n');
fprintf('q5b)... out of 10\n\n');



fprintf('________________________________________\n\n')
fprintf('Total .... out of 55   \n\n\n')



%%
%
% Question 1
%
clear
%
errsol=0;
% List student submission
if exist('q1.m', 'file')==2
    fprintf('\n========== Submitted q1.m ========\n')
    type q1.m
    try
        fprintf('\n\n\n---------- Running submitted q1.m --------\n')
        q1
    catch MSG
        fprintf('ERROR running q1.m:\n');
        MSG.message
        errsol=1;   % q1 is not marked if there is an error
        %quit
    end
else
    fprintf('\n ======== NO file q1.m submitted =======\n');
    %quit
end;

marks = 0;


if exist('q1.m', 'file')==2 && errsol == 0

fprintf('\n\n\n');
fprintf('================ q1 marking program ============\n\n\n');


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

if exist('n1', 'var') == 1
    fprintf('abs(n1-n1sol) = %e\n',norm(n1-n1sol,Inf));
    if abs(n1-n1sol) < 1e-15
        marks = marks + 1;
        fprintf('q1)a) n1 is within the specified tolerance\n\n');
    else
        fprintf('q1)a) n1 is not within the specified tolerance\n\n');
    end;
else
    fprintf('q1)a) n1 is not defined\n\n');
end;

if exist('A1', 'var') == 1
    if prod(size(A1) == size(A1sol))
        fprintf('norm(A1-A1sol,inf) = %e\n',norm(A1-A1sol,Inf));
        if norm(A1-A1sol,Inf) < 1e-6
          marks = marks + 2;
          fprintf('q1)b) A1 is within the specified tolerance\n\n');
        else
          fprintf('q1)b) A1 is not within the specified tolerance\n\n');
        end;
    else
        fprintf('q1)b) A1 is not of correct size\n\n');
    end;
else
    fprintf('q1)b) A1 is not defined\n\n');
end;

if exist('Q1', 'var') == 1
    if prod(size(Q1) == size(Q1sol))
        fprintf('norm(Q1-Q1sol,inf) = %e\n',norm(Q1-Q1sol,Inf));
        if norm(Q1-Q1sol,Inf) < 1e-6
          marks = marks + 1;
          fprintf('q1)c) Q1 is within the specified tolerance\n\n');
        else
          fprintf('q1)c) Q1 is not within the specified tolerance\n\n');
        end;
    else
        fprintf('q1)c) Q1 is not of correct size\n\n');
    end;
else
    fprintf('q1)c) Q1 is not defined\n\n');
end;



if exist('Q1check', 'var') == 1
    if size(Q1check) == size(Q1checksol)
        fprintf('abs(Q1check-Q1checksol) = %e\n',norm(Q1check-Q1checksol,Inf));
        if norm(Q1check-Q1checksol,Inf) < 1e-6
          marks = marks + 1;
          fprintf('q1)d) Q1check is within the specified tolerance\n\n');
        else
          fprintf('q1)d) Q1check is not within the specified tolerance\n\n');
        end;
    else
        fprintf('q1)d) Q1check is not of correct size\n\n');
    end;
else
    fprintf('q1)d) Q1check is not defined\n\n');
end;

if exist('b', 'var') == 1
    if prod(size(b) == size(bsol))
        fprintf('norm(b-bsol,inf) = %e\n',norm(b-bsol,Inf));
        if norm(b-bsol,Inf) < 1e-6
          marks = marks + 1;
          fprintf('q1)e) b is within the specified tolerance\n\n');
        else
          fprintf('q1)e) b is not within the specified tolerance\n\n');
        end;
    else
        fprintf('q1)e) b is not of correct size\n\n');
    end;
else
    fprintf('q1)e) b is not defined\n\n');
end;

if exist('xsol', 'var') == 1
    if prod(size(xsol) == size(xsolsol))
        fprintf('norm(xsol-xsolsol,inf) = %e\n',norm(xsol-xsolsol,Inf));
        if norm(xsol-xsolsol,Inf) < 1e-4
          marks = marks + 2;
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


if exist('rsol', 'var') == 1
    if prod(size(rsol) == size(rsolsol))
        fprintf('norm(rsol-rsolsol,inf) = %e\n',norm(rsol-rsolsol,Inf));
        if norm(rsol-rsolsol,Inf) < 1e-4
          marks = marks + 1;
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
          marks = marks + 1;
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


fprintf('\n\n');

fprintf('q1.m suggested marks: %d out of 10 marks.\n\n', marks);
fprintf('The final marks will depend on how you obtained the solution\n');
fprintf('and may vary from the suggested marks.\n');
end;

marks1 = marks;

%%
% Question 2
%
%
clear all
close all
%
errsol=0;
% List student submission
if exist('q2.m', 'file')==2
    fprintf('\n========== Submitted q2.m ========\n')
    type q2.m
    try
        fprintf('\n\n\n---------- Running submitted q2.m --------\n')
        q2
    catch MSG
        fprintf('ERROR running q2.m:\n');
        MSG.message
        errsol=1;
        %quit
    end
else
    fprintf('\n ======== NO file q2.m submitted =======\n');
    %quit
end;

marks = 0;

if exist('q2.m', 'file')==2 && errsol==0

fprintf('\n\n\n');
fprintf('================ q2 marking program ============\n\n\n');


% (a)
x0sol = 1;
% (b)
fsol = @(x) x.^2 - exp(-x.^2/4);
dfsol = @(x) 2*x + x/2 .* exp(-x.^2/4);
% (c)
x1sol = x0sol - fsol(x0sol)/dfsol(x0sol);
% (d)
%fprintf('x1sol = %e \n',x1);
% (f)
%fprintf('f1(x1sol) - f2(x1) = %e \n',f(x1));
% (e)
xnewsol = x1sol;
xoldsol = x0sol;
while abs(xnewsol - xoldsol) > 2^(-52)
    xoldsol = xnewsol;
    xnewsol = xnewsol - fsol(xnewsol)/dfsol(xnewsol);
end
xappsol = xnewsol;


fprintf('q2)a) Has a plot been created to find x0? (+1 mark) \n');

fprintf('\n');

if exist('x0', 'var') == 1
    fprintf('q2)a) x0 exists and x0 = %e\n', x0);
 %   if abs(x0-x0sol) < 10
        marks = marks + 1;
%        fprintf('q2)a) x0 is within the specified tolerance\n');
%    else
%        fprintf('q2)a) x0 is not within the specified tolerance\n');
%    end;
else
    fprintf('q2)a) x0 is not defined\n');
end;
fprintf('\n');

fprintf('q2)b) Has an autonomous function f been defined? (+1 mark) \n');
fprintf('\n');

if exist('x1', 'var') == 1
        fprintf('abs(x0-xsol) - abs(x1-xsol) = %e\n', abs(x0-xappsol) - abs(x1-xappsol));
        if abs(x0-xappsol) > abs(x1-xappsol)
          marks = marks + 1;
          fprintf('q2)c) x1 is closer to x* than x0 \n');
        else
          fprintf('q2)c) x1 is not closer to x* than x0 \n');
        end;
else
    fprintf('q2)c) x1 is not defined\n');
end;

fprintf('\n');

fprintf('q2)d) Is x1 printed using 3 significant figures? (+1 mark) \n')
fprintf('\n');
fprintf('q2)e) Is |f1(x1) - f2(x1)| printed? (+1 mark) \n')
fprintf('\n');
if exist('xapp', 'var') == 1
    if prod(size(xapp) == size(xappsol))
        fprintf('abs(xapp-xappsol) = %e\n',norm(xapp-xappsol,Inf));
        if norm(xapp-xappsol,Inf) < 1e-8
          marks = marks + 4;
          fprintf('q2)f) xapp is within the specified tolerance\n');
        else
          fprintf('q2)f) xapp is not within the specified tolerance\n');
        end;
    else
        fprintf('q2)f) xapp is not of correct size\n');
    end;
else
    fprintf('q2)f) xapp is not defined\n');
end;


fprintf('\n\n');

fprintf('q2.m suggested marks: %d out of 6 marks marked by this program.\n\n', marks);
fprintf('The final marks will depend on how you obtained the solution\n');
fprintf('and may vary from the suggested marks.\n');


end;

marks2 = marks;

%%
% Question 3
%
clear all
close all
%
errsol=0;
% List student submission
if exist('q3.m', 'file')==2
    fprintf('\n========== Submitted q3.m ========\n')
    type q3.m
    try
        fprintf('\n\n\n---------- Running submitted q3.m --------\n')
        q3
    catch MSG
        fprintf('ERROR running q3.m:\n');
        MSG.message
        errsol=1;
        %quit
    end
else
    fprintf('\n ======== NO file q3.m submitted =======\n');
    %quit
end;

marks = 0;

if exist('q3.m', 'file')==2 && errsol==0

fprintf('\n\n\n');
fprintf('================ q3 marking program ============\n\n\n');


% (a)
load 5data19.mat
tdatasol = tdata(:);
ydatasol = ydata(:);
% (b)
ndatasol = length(tdata);
% (c)
Asol = [cos(2*pi*tdatasol) cos(4*pi*tdatasol) tdatasol.^2];
xlssol = Asol\ydatasol;
% (d)
ppsol = spline(tdatasol,ydatasol);
xplsol = linspace(0,1,201);
yplsol = ppval(ppsol, xplsol);
% (e)
%figure(1)
%plot(tdata,ydata,'*', xpl, ypl, xpl, xls(1)*cos(2*pi*xpl)+xls(2)*cos(4*pi*xpl)+xls(3)*xpl.^2)
%legend('Data', 'Spline', 'Least squares app');

%q3

marks = 0;


if exist('ndata', 'var') == 1
            if norm(ndata-ndatasol) < 10^(-8)
                fprintf('q3)b) ndata exists and is correct\n')
                marks = marks + 1;
            else
                fprintf('q3)b) ndata exists but is not correct\n')
                marks = marks + 0;
            end
else
    fprintf('q3)b) ndata is not defined\n')
end

fprintf('\n');

if exist('xls', 'var') == 1
        if size(xls) == size(xlssol)
            if norm(xls-xlssol) < 10^(-3)
                fprintf('q3)c) xls exists and is correct\n')
                marks = marks + 3;
            else
                fprintf('q3)c) xls exists but is not correct\n')
                marks = marks + 0;
            end
        else
            fprintf('q3)d) xls is not of the correct size\n')
        end
end
        
fprintf('\n');
fprintf('q3)d) Has the spline been calculated using the not-a-knot condition? (+2 marks) \n')
fprintf('\n');
fprintf('q3)e) Has a plot been created which includes the data, the least squares approximation, the spline and a legend? (+4 marks) \n')


fprintf('\n\n');


fprintf('q3.m suggested marks: %d out of 4 marks checked by this program.\n\n', marks);
fprintf('The final marks will depend on how you obtained the solution\n');
fprintf('and may vary from the suggested marks.\n');


end;

marks3 = marks;


%% Question 4
%
clear all
close all
%
errsol=0;
% List student submission
if exist('q4.m', 'file')==2
    fprintf('\n========== Submitted q4.m ========\n')
    type q4.m
    try
        fprintf('\n\n\n---------- Running submitted q4.m --------\n')
        q4
    catch MSG
        fprintf('ERROR running q4.m:\n');
        MSG.message
        errsol=1;
        %quit
    end
else
    fprintf('\n ======== NO file q4.m submitted =======\n');
    %quit
end;


marks = 0;

if exist('q4.m', 'file')==2 && errsol==0

fprintf('\n\n\n');
fprintf('================ q4 marking program ============\n\n\n');

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





fprintf('q4)a) Is the function f defined correctly? (+ 1 mark)\n')
fprintf('\n');
if exist('IQ', 'var') == 1
    fprintf('abs(IQ-IQsol) = %e\n',norm(IQ-IQsol,Inf));
    if abs(IQ-IQsol) < 1e-6
        marks = marks + 1;
        fprintf('q4)b) IQ is within the specified tolerance\n');
    else
        fprintf('q4)b) IQ is not within the specified tolerance\n');
    end
else
    fprintf('q4)b) IQ is not defined\n');
end
fprintf('\n');
if exist('NI', 'var') == 1
        fprintf('abs(NI-NIsol) = %d\n',norm(NI-NIsol,Inf));
        if norm(NI-NIsol,Inf) < 1e-10
          marks = marks + 1;
          fprintf('q4)c)i) NI is within the specified tolerance\n');
        else
          fprintf('q4)c)i) NI is not within the specified tolerance\n');
        end
else
    fprintf('q4)c)i) NI is not defined\n');
end
fprintf('\n');
xI = xI(:);

if exist('xI', 'var') == 1
    %if prod(size(xI) == size(xIsol))
        fprintf('abs(xI-xIsol) = %e\n',norm(xI-xIsol,Inf));
        if norm(xI-xIsol,Inf) < 1e-7
          marks = marks + 1;
          fprintf('q4)c)ii) xI is within the specified tolerance\n');
        else
          fprintf('q4)c)ii) xI is not within the specified tolerance\n');
        end
    %else
    %    fprintf('q4)c)ii) xI is not of correct size\n');
    %end;
else
    fprintf('q4)c)ii) xI is not defined\n');
end;
fprintf('\n');

if exist('IT', 'var') == 1
        fprintf('abs(IT-ITsol) = %d\n',norm(IT-ITsol,Inf));
        if norm(IT-ITsol,Inf) < 1e-4
          marks = marks + 2;
          fprintf('q4)c)iii) IT is within the specified tolerance\n');
        else
          fprintf('q4)c)iii) IT is not within the specified tolerance\n');
        end;
else
    fprintf('q4)c)iii) IT is not defined\n');
end;
fprintf('\n');

fprintf('q4)d) Is the value of |IQ-IT| printed? (+ 1 mark)\n');
fprintf('q4)d) My value: |IQ - IT| = %e\n',IQsol-ITsol);

fprintf('\n');

if exist('N2', 'var') == 1
        fprintf('N2-N2sol = %d\n', N2-N2sol);
        if N2 > N2sol-1
          marks = marks + 1;
          fprintf('q4)c)iii) N2 is within the specified tolerance\n');
        else
          fprintf('q4)c)iii) N2 is not within the specified tolerance\n');
        end;
else
    fprintf('q4)c)iii) N2 is not defined\n');
end;

fprintf('\n');

if exist('IT2', 'var') == 1
        fprintf('abs(IQsol-IT2) = %d\n',norm(IT2-IQsol,Inf));
        fprintf('abs(IQ-IT2) = %d\n',norm(IQ-IT2,Inf));
        if norm(IT2-IQsol,Inf) < 0.001
          marks = marks + 1;
          fprintf('q4)e) IT2 is within the specified tolerance\n');
        else
          fprintf('q4)e) IT2 is not within the specified tolerance\n');
        end;
else
    fprintf('q4)e) IT2 is not defined\n');
end;

fprintf('\n');

fprintf('q4)f) Is there an explanation for the large number of nodes needed? (+1 mark) \n')



fprintf('\n\n');

fprintf('q4.m suggested marks: %d out of 7 marks checked by this program.\n\n', marks);
fprintf('The final marks will depend on how you obtained the solution\n');
fprintf('and may vary from the suggested marks.\n');


end

marks4 = marks;

%%
% Question 5 a)
%
clear
%
errsol = 0;
% List student submission
if exist('F519.m', 'file')==2
    fprintf('\n========== Submitted F519.m ========\n')
    type F519.m
    try
        fprintf('\n\n\n---------- Running submitted F519.m --------\n')
        tsol = F519([1 2 3 4]);
    catch MSG
        fprintf('ERROR running F519.m:\n');
        MSG.message
        errsol = 1;
        %quit
    end
else
    fprintf('\n ======== NO file F519.m submitted =======\n');
    %quit
end;

marks = 0;

if exist('F519.m', 'file')==2 && errsol == 0

fprintf('\n\n\n');
fprintf('================ F519 marking program ============\n\n\n');

fprintf('q5)a)i) Is the function F519.m defined correctly? (+2 marks)\n');

fprintf('\n');

tsol = F519([1,2,3,4]);
qsol = size(tsol);
if qsol(1) ~= 4
    fprintf('q5)a)ii) Output is not a column vector.\n');
else
    fprintf('q5)a)ii) F519([1,2,3,4]) - (4,8,14,22)^T = %e\n', norm(tsol-[4;8;14;22]));
    if norm(tsol-[4;8;14;22]) < 10^(-7)
        marks = marks + 2;
        fprintf('q5)a)ii) Output seems to be correct\n');
    end;
end;

fprintf('\n');

fprintf('q5)a)iii) Does help F519 produce explanation? (+1 mark) \n');
help F519


fprintf('\n\n');

fprintf('F519.m suggested marks: %d out of 2 marks checked by the marking program.\n\n', marks);
fprintf('The final marks will depend on how you obtained the solution\n');
fprintf('and may vary from the suggested marks.\n');


end;



marks5a = marks;

%%
% Question 5b)
%
clear all
close all
%
errsol=0;
% List student submission
if exist('q5.m', 'file')==2
    fprintf('\n========== Submitted q5.m ========\n')
    type q5.m
    try
        fprintf('\n\n\n---------- Running submitted q5.m --------\n')
        q5
    catch MSG
        fprintf('ERROR running q5.m:\n');
        MSG.message
        errsol=1;
        %quit
    end
else
    fprintf('\n ======== NO file q5.m submitted =======\n');
    %quit
end;

marks = 0;

if exist('q5.m', 'file')==2 && errsol == 0

fprintf('\n\n\n');
fprintf('================ q5 marking program ============\n\n\n');


% (a)
Nsol = 10^2;
% (b)
Usol = zeros(Nsol+1,1);
Usol(2) = 2/Nsol;
for nsol = 2:Nsol
    Usol(nsol+1) = -Usol(nsol-1)+Usol(nsol)*(2-Nsol^(-2)) + (2 + nsol/Nsol + nsol^2/Nsol^2)/Nsol^2;
end
% Solution
xsol = linspace(0,1,Nsol+1);
usol = sin(xsol) + xsol + xsol.^2;

% Error
Uerrsol = norm(usol-Usol');

%
Nlargesol = 10^6;
Ulargesol = zeros(Nlargesol+1,1);
Ulargesol(2) = 2/Nlargesol;
for nsol = 2:Nlargesol
    Ulargesol(nsol+1) = -Ulargesol(nsol-1)+Ulargesol(nsol)*(2-Nlargesol^(-2)) + (2 + nsol/Nlargesol + nsol^2/Nlargesol^2)/Nlargesol^2;
end
Ulargesol=Ulargesol(:);

fprintf('\n');

if exist('N', 'var') == 1
    fprintf('abs(N-Nsol) = %e\n',norm(N-Nsol,Inf));
    if abs(N-Nsol) < 1e-10
        marks = marks + 1;
        fprintf('q5)b)i) N is within the specified tolerance\n');
    else
        fprintf('q5)b)i) N is not within the specified tolerance\n');
    end;
else
    fprintf('q5)b)i) N is not defined\n');
end;

fprintf('\n');

if exist('U', 'var') == 1
    if prod(size(U) == size(Usol))
        fprintf('norm(U-Usol,Inf) = %e\n',norm(U-Usol,Inf));
        if norm(U-Usol,Inf) < 1e-1
          marks = marks + 4;
          fprintf('q5)b)ii) U is within the specified tolerance\n');
        else
          fprintf('q5)b)ii) U is not within the specified tolerance\n');
        end;
    else
        fprintf('q5)b)ii) U is not of correct size\n');
    end;
else
    fprintf('q5)b)ii) U is not defined\n');
end;
       fprintf('\n');

if exist('Uerr', 'var') == 1
    fprintf('Uerr = %e\n', Uerr);
    if Uerr < 1e-1
        marks = marks + 1;
        fprintf('q5)b)iii) Uerr is within the specified tolerance\n');
    else
        fprintf('q5)b)iii) Uerr is not within the specified tolerance\n');
    end;
else
    fprintf('q5)b)iii) Uerr is not defined\n');
end;   
       fprintf('\n');   
if exist('Nlarge', 'var') == 1
    fprintf('abs(Nlarge-Nlargesol) = %e\n',norm(Nlarge-Nlargesol,Inf));
    if abs(Nlarge-Nlargesol) < 1e-10
        marks = marks + 1;
        fprintf('q5)b)iv) Nlarge is within the specified tolerance\n');
    else
        fprintf('q5)b)iv) Nlarge is not within the specified tolerance\n');
    end;
else
    fprintf('q5)b)iv) Nlarge is not defined\n');
end;

fprintf('\n');
if exist('Ularge', 'var') == 1
    if prod(size(Ularge) == size(Ulargesol))
        fprintf('norm(Ularge-Ulargesol,Inf) = %e\n',norm(Ularge-Ulargesol,Inf));
        if norm(Ularge-Ulargesol,Inf) < 1e-1
          marks = marks + 3;
          fprintf('q5)b)v) Ularge is within the specified tolerance\n');
        else
          fprintf('q5)b)v) Ularge is not within the specified tolerance\n');
        end;
    else
        fprintf('q5)b)v) Ularge is not of correct size\n');
    end;
else
    fprintf('q5)b)v) Ularge is not defined\n');
end;




fprintf('\n\n');

fprintf('q5.m suggested marks: %d out of 10 marks.\n\n', marks);

fprintf('The final marks will depend on how you obtained the solution\n');
fprintf('and may vary from the suggested marks.\n');


end;

marks5b = marks;


%fprintf('Summary: \n')
%fprintf('Question 1:  %d\n', marks1)
%fprintf('Question 2:  %d\n', marks2)
%fprintf('Question 3:  %d\n', marks3)
%fprintf('Question 4:  %d\n', marks4)
%fprintf('Question 5a: %d\n', marks5a)
%fprintf('Question 5b: %d\n', marks5b)

%totalmarks = marks1 + marks2 + marks3 + marks4 + marks5a + marks5b
%fprintf('Total marks: %d\n', totalmarks)



diary Test5335.txt

diary off

clear