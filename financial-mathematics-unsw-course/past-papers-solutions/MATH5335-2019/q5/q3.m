%(a) Load q3val.mat
load q3val.mat;

tdatasol = tdata;
ydatasol = ydata;

clear tdata;
clear ydata;

%(b) Find ndata;
ndatasol = length(tdatasol);

%(c) Calculate least squares approximation
tdatasol = tdatasol(:); % Make sure tdata is a column vector;
ydatasol = ydatasol(:);
Asol = [ones(ndatasol,1) -tdatasol.^2  log(tdatasol)];
xlssol = Asol \ log(ydatasol);
xlssol(1) = exp(xlssol(1));

%(d) Print coefficients
%fprintf('\n x(1) = %1.4f, x(2) = %1.4f, x(3) = %1.4f \n\n', xlssol(1), xlssol(2), xlssol(3));

%(e) Calculate cubic spline
ppsol = spline(tdatasol,ydatasol);

%(f) Plot
tpsol = linspace(1,2,1001);

ylssol = xlssol(1)*exp(-xlssol(2)*tpsol.^2).* tpsol.^xlssol(3);
yysol = ppval(ppsol,tpsol);

%figure(2)
%plot(tdatasol,ydatasol,'*',tpsol,ylssol,tpsol,yysol);
%legend('Data','Least squares fit','Cubic spline','Location','NorthWest');
