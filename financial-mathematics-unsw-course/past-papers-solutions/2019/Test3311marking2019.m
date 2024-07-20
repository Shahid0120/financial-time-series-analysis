% Labtest 2019 Math3311 marking program


function Test5335marking

diary off
clear all
close all
format compact
format long e

if exist('Test3311.txt', 'file') == 2
    delete Test3311.txt
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
nsol = 10;

% (b)
Asol = diag(-2*ones(nsol,1), 0) + diag(ones(nsol-2,1), 2) + diag(3*ones(nsol-4,1), -4);

% (c)
[Lsol, Usol, psol] = lu(Asol, 'vector');

% (d)
bsol = [2*n:-2:2]';

% (e)
ysol = Lsol \ bsol(psol);
xssol = Usol \ ysol;

% (f)
rsol = Asol*xssol - bsol;

% (g)
r1normsol = norm(rsol,1);


if exist('n', 'var') == 1
    fprintf('abs(n-nsol) = %e\n',norm(n-nsol,Inf));
    if abs(n-nsol) < 1e-15
        marks = marks + 1;
        fprintf('q1)a) n is within the specified tolerance\n\n');
    else
        fprintf('q1)a) n is not within the specified tolerance\n\n');
    end;
else
    fprintf('q1)a) n is not defined\n\n');
end;

if exist('A', 'var') == 1
    if prod(size(A) == size(Asol))
        fprintf('norm(A-Asol,inf) = %e\n',norm(A-Asol,Inf));
        if norm(A-Asol,Inf) < 1e-6
          marks = marks + 2;
          fprintf('q1)b) A is within the specified tolerance\n\n');
        else
          fprintf('q1)b) A is not within the specified tolerance\n\n');
        end;
    else
        fprintf('q1)b) A is not of correct size\n\n');
    end;
else
    fprintf('q1)b) A is not defined\n\n');
end;

if exist('L', 'var') == 1
    if prod(size(L) == size(Lsol))
        fprintf('norm(L-Lsol,inf) = %e\n',norm(L-Lsol,Inf));
        if norm(L-Lsol,Inf) < 1e-6
          marks = marks + 1;
          fprintf('q1)c) L is within the specified tolerance\n\n');
        else
          fprintf('q1)c) L is not within the specified tolerance\n\n');
        end;
    else
        fprintf('q1)c) L is not of correct size\n\n');
    end;
else
    fprintf('q1)c) L is not defined\n\n');
end;


if exist('p', 'var') == 1
    if prod(size(p) == size(psol))
        fprintf('norm(p-psol,inf) = %e\n',norm(p-psol,Inf));
        if norm(p-psol,Inf) < 1e-6
          marks = marks + 1;
          fprintf('q1)c) p is within the specified tolerance\n\n');
        else
          fprintf('q1)c) p is not within the specified tolerance\n\n');
        end;
    else
        fprintf('q1)c) p is not of correct size\n\n');
    end;
else
    fprintf('q1)c) p is not defined\n\n');
end;


if exist('b', 'var') == 1
    if prod(size(b) == size(bsol))
        fprintf('norm(b-bsol,inf) = %e\n',norm(b-bsol,Inf));
        if norm(b-bsol,Inf) < 1e-6
          marks = marks + 1;
          fprintf('q1)d) b is within the specified tolerance\n\n');
        else
          fprintf('q1)d) b is not within the specified tolerance\n\n');
        end;
    else
        fprintf('q1)d) b is not of correct size\n\n');
    end;
else
    fprintf('q1)d) b is not defined\n\n');
end;

if exist('xs', 'var') == 1
    if prod(size(xs) == size(xssol))
        fprintf('norm(xs-xssol,inf) = %e\n',norm(xs-xssol,Inf));
        if norm(xs-xssol,Inf) < 1e-4
          marks = marks + 2;
          fprintf('q1)e) xs is within the specified tolerance\n\n');
        else
          fprintf('q1)e) xs is not within the specified tolerance\n\n');
        end;
    else
        fprintf('q1)e) xs is not of correct size\n\n');
    end;
else
    fprintf('q1)f) xs is not defined\n\n');
end;


if exist('r', 'var') == 1
    if prod(size(r) == size(rsol))
        fprintf('norm(r-rsol,inf) = %e\n',norm(r-rsol,Inf));
        if norm(r-rsol,Inf) < 1e-4
          marks = marks + 1;
          fprintf('q1)f) r is within the specified tolerance\n\n');
        else
          fprintf('q1)f) r is not within the specified tolerance\n\n');
        end;
    else
        fprintf('q1)f) r is not of correct size\n\n');
    end;
else
    fprintf('q1)f) r is not defined\n\n');
end;

if exist('r1norm', 'var') == 1
    if prod(size(r1norm) == size(r1normsol))
        fprintf('norm(r1norm-r1normsol,inf) = %e\n',norm(r1norm-r1normsol,Inf));
        if norm(r1norm-r1normsol,Inf) < 1e-4
          marks = marks + 1;
          fprintf('q1)g) r1norm is within the specified tolerance\n\n');
        else
          fprintf('q1)g) r1norm is not within the specified tolerance\n\n');
        end;
    else
        fprintf('q1)g) r1norm is not of correct size\n\n');
    end;
else
    fprintf('q1)g) r1norm is not defined\n\n');
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
x0sol = 3;
x1sol = 2.5;

% (b)
gsol = @(x) (x-1).^3 + exp(x) - 2;
x2sol = x1sol - gsol(x1sol) * (x1sol-x0sol ) / (gsol(x1sol)-gsol(x0sol));

% (c)
errsol = 1;
xoldsol = x1sol;
xoldoldsol = x0sol;
while errsol > eps
    xnewsol = xoldsol - gsol(xoldsol)* (xoldsol - xoldoldsol) / (gsol(xoldsol) - gsol(xoldoldsol));
    errsol = abs(xnewsol - xoldsol);
    xoldoldsol = xoldsol;
    xoldsol = xnewsol;
end;
xesol = xnewsol;



if exist('x0', 'var') == 1
    fprintf('q2)a) x0 exists and x0 = %e\n', x0);
    if abs(x0-x0sol) < 1e-6
        marks = marks + 1;
        fprintf('q2)a) x0 is within the specified tolerance\n');
    else
        fprintf('q2)a) x0 is not within the specified tolerance\n');
    end;
else
    fprintf('q2)a) x0 is not defined\n');
end;
fprintf('\n');

if exist('x1', 'var') == 1
    fprintf('q2)a) x1 exists and x1 = %e\n', x1);
    if abs(x1-x1sol) < 1e-6
        marks = marks + 1;
        fprintf('q2)a) x1 is within the specified tolerance\n');
    else
        fprintf('q2)a) x1 is not within the specified tolerance\n');
    end;
else
    fprintf('q2)a) x1 is not defined\n');
end;
fprintf('\n');

if exist('x2', 'var') == 1
        fprintf('abs(x2-x2sol) = %e\n', abs(x2-x2sol));
        if abs(x2-x2sol) < 1e-3
          marks = marks + 3;
          fprintf('q2)b) x2 is within the specified tolerance \n');
        else
          fprintf('q2)b) x2 is not within the specified tolerance \n');
        end;
else
    fprintf('q2)b) x2 is not defined\n');
end;

fprintf('\n');

if exist('xe', 'var') == 1
        fprintf('abs(xe-xesol) = %e\n', abs(xe-xesol));
        if abs(xe-xesol) < 1e-3
          marks = marks + 4;
          fprintf('q2)c) xe is within the specified tolerance \n');
        else
          fprintf('q2)c) xe is not within the specified tolerance \n');
        end;
else
    fprintf('q2)c) xe is not defined\n');
end;

fprintf('\n');


fprintf('q2)d) Is x2 and/or xe printed? (+1 mark) \n')


fprintf('\n');


% (d)
fprintf('My solution x2sol = %f, xesol = %f \n', x2sol, xesol)




fprintf('\n\n');

fprintf('q2.m suggested marks: %d out of 9 marks marked by this program.\n\n', marks);
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

if exist('x', 'var') == 1
    xsub = x;
end
if exist('y', 'var') == 1
    ysub = y;
end

clear x, y;


marks = 0;

% (a)
load q319data.mat


if exist('xsub', 'var') == 1
        if size(xsub) == size(x)
            if norm(xsub-x) < 10^(-3)
                fprintf('q3)a) x exists and is correct\n')
                marks = marks + 1;
            else
                fprintf('q3)a) x exists but is not correct\n')
                marks = marks + 0;
            end
        else
            fprintf('q3)a) x is not of the correct size\n')
        end
end

fprintf('\n')

% (a)
load q319data.mat

% (b)
%figure(1)
%plot(x,y,'*')

% (c)
Asol  = [ones(size(x)) x x.^2 x.^3 x.^4];
zsol = Asol\y;

% (d)
xplsol = linspace(0,1,401);
yplsol = zsol(1)*ones(size(xplsol)) + zsol(2)*xplsol + zsol(3)*xplsol.^2 + zsol(4)*xplsol.^3+zsol(5)*xplsol.^4;
%hold on
%plot(xpl, ypl)
%hold off

% (e)
lserrsol = norm(Asol*zsol-y);

% (f)
%fprintf('lserrsol = %f\n\n', lserrsol)




fprintf('q3)b) Is there a plot of the data? (+1 marks) \n')

fprintf('\n')

if exist('z', 'var') == 1
        if size(z) == size(zsol)
            if norm(z-zsol) < 10^(-3)
                fprintf('q3)c) z exists and is correct\n')
                marks = marks + 4;
            else
                fprintf('q3)c) z exists but is not correct\n')
                marks = marks + 0;
            end
        else
            fprintf('q3)c) z is not of the correct size\n')
        end
end
        
fprintf('\n');
fprintf('q3)d) Is a plot saved as Figure1.pdf which includes the data points and least squares? (+2 marks) \n')
fprintf('\n');




if exist('lserr', 'var') == 1
        fprintf('abs(lserr-lserrsol) = %e\n', abs(lserr-lserrsol));
        if abs(lserr-lserrsol) < 1e-3
          marks = marks + 1;
          fprintf('q3)e) lserr is within the specified tolerance \n');
        else
          fprintf('q3)e) lserr is not within the specified tolerance \n');
        end;
else
    fprintf('q3)e) lserr is not defined\n');
end;

fprintf('\n');

fprintf('q3)f) Has the value of lserr bin printed? (+1 mark) \n')

% (f)
fprintf('My value lserrsol = %f\n\n', lserrsol)


fprintf('\n\n');


fprintf('q3.m suggested marks: %d out of 6 marks checked by this program.\n\n', marks);
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
fsol = @(x) abs(cos(pi*x));

% (b)
%IQ = quad(f, 0, 2);
IQsol = integral(f,0,2);

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
xIsol = xIsol(:);

if exist('xI', 'var') == 1
    if prod(size(xI) == size(xIsol))
        fprintf('abs(xI-xIsol) = %e\n',norm(xI-xIsol,Inf));
        if norm(xI-xIsol,Inf) < 1e-7
          marks = marks + 1;
          fprintf('q4)c)ii) xI is within the specified tolerance\n');
        else
          fprintf('q4)c)ii) xI is not within the specified tolerance\n');
        end
    else
        fprintf('q4)c)ii) xI is not of correct size\n');
    end;
else
    fprintf('q4)c)ii) xI is not defined\n');
end;
fprintf('\n');


if exist('hI', 'var') == 1
        fprintf('abs(hI-hIsol) = %d\n',norm(hI-hIsol,Inf));
        if norm(hI-hIsol,Inf) < 1e-6
          marks = marks + 1;
          fprintf('q4)c)iii) hI is within the specified tolerance\n');
        else
          fprintf('q4)c)iii) hI is not within the specified tolerance\n');
        end
else
    fprintf('q4)c)iii) hI is not defined\n');
end
fprintf('\n');

if exist('IT', 'var') == 1
        fprintf('abs(IT-ITsol) = %d\n',norm(IT-ITsol,Inf));
        if norm(IT-ITsol,Inf) < 1e-4
          marks = marks + 3;
          fprintf('q4)c)iv) IT is within the specified tolerance\n');
        else
          fprintf('q4)c)iv) IT is not within the specified tolerance\n');
        end;
else
    fprintf('q4)c)iv) IT is not defined\n');
end;
fprintf('\n');

fprintf('q4)d) (+ 2 marks)\n');
fprintf('q4)d) Using 1001 points means that there is one point at x = 1.\n')
fprintf('Thus in this case we are using two trapezoidal rules, one for the\n')
fprintf('interval [0,1] with 501 points and another one for the interval [1,2],\n')
fprintf('also with 501 points. Note that the point at 1 has weight h/2 for each rule.\n')
fprintf('Thus we avoid the singularity at x = 1;')




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
if exist('F319.m', 'file')==2
    fprintf('\n========== Submitted F319.m ========\n')
    type F319.m
    try
        fprintf('\n\n\n---------- Running submitted F319.m --------\n')
        tsol = F319([1 2 3 4]);
    catch MSG
        fprintf('ERROR running F319.m:\n');
        MSG.message
        errsol = 1;
        %quit
    end
else
    fprintf('\n ======== NO file F319.m submitted =======\n');
    %quit
end;

marks = 0;

if exist('F319.m', 'file')==2 && errsol == 0

fprintf('\n\n\n');
fprintf('================ F319 marking program ============\n\n\n');

fprintf('q5)a)i) Is the function F319.m defined correctly? (+2 marks)\n');

fprintf('\n');

tsol = F319([1,2,3,4]);
qsol = size(tsol);
if qsol(1) ~= 4
    fprintf('q5)a)ii) Output is not a column vector.\n');
else
    fprintf('q5)a)ii) F319([1,2,3,4]) - (2*exp(1),12*exp(2),36*exp(3),80*exp(4))^T = %e\n', norm(tsol-[2*exp(1);12*exp(2);36*exp(3);80*exp(4)]));
    if norm(tsol-[2*exp(1);12*exp(2);36*exp(3);80*exp(4)]) < 10^(-4)
        marks = marks + 2;
        fprintf('q5)a)ii) Output seems to be correct\n');
    end;
end;

fprintf('\n');

fprintf('q5)a)iii) Does help F319 produce explanation? (+1 mark) \n');
help F319


fprintf('\n\n');

fprintf('F319.m suggested marks: %d out of 2 marks checked by the marking program.\n\n', marks);
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



if exist('B', 'var') == 1
    if prod(size(B) == size(Bsol))
        fprintf('abs(B-Bsol) = %e\n',norm(B-Bsol,Inf));
        if norm(B-Bsol,Inf) < 1e-4
          marks = marks + 3;
          fprintf('q5)b)iii) B is within the specified tolerance\n');
        else
          fprintf('q5)b)iii) B is not within the specified tolerance\n');
        end
    else
        fprintf('q5)b)iii) B is not of correct size\n');
    end;
else
    fprintf('q5)b)iii) B is not defined\n');
end;
fprintf('\n');



if exist('b', 'var') == 1
    if prod(size(b) == size(bsol))
        fprintf('abs(b-bsol) = %e\n',norm(b-bsol,Inf));
        if norm(b-bsol,Inf) < 1e-4
          marks = marks + 3;
          fprintf('q5)b)iii) b is within the specified tolerance\n');
        else
          fprintf('q5)b)iii) b is not within the specified tolerance\n');
        end
    else
        fprintf('q5)b)iii) b is not of correct size\n');
    end;
else
    fprintf('q5)b)iii) b is not defined\n');
end;
fprintf('\n');


if exist('U', 'var') == 1
    if prod(size(U) == size(Usol))
        fprintf('norm(U-Usol,Inf) = %e\n',norm(U-Usol,Inf));
        if norm(U-Usol,Inf) < 1e-1
          marks = marks + 1;
          fprintf('q5)b)iii) U is within the specified tolerance\n');
        else
          fprintf('q5)b)iii) U is not within the specified tolerance\n');
        end;
    else
        fprintf('q5)b)iii) U is not of correct size\n');
    end;
else
    fprintf('q5)b)iii) U is not defined\n');
end;
       fprintf('\n');

if exist('Uerr', 'var') == 1
    fprintf('Uerr = %e\n', Uerr);
    if Uerr < 2e-1
        marks = marks + 2;
        fprintf('q5)b)iv) Uerr is within the specified tolerance\n');
    else
        fprintf('q5)b)iv) Uerr is not within the specified tolerance\n');
    end;
else
    fprintf('q5)b)iv) Uerr is not defined\n');
end;   







fprintf('\n\n');

fprintf('q5.m suggested marks: %d out of 10 marks.\n\n', marks);

fprintf('The final marks will depend on how you obtained the solution\n');
fprintf('and may vary from the suggested marks.\n');


end;

marks5b = marks;



diary Test3311.txt

diary off

clear