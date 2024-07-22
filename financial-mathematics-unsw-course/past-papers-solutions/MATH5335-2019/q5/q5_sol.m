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