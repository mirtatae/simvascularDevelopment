function v = velEccCylinders(x,y,rv,rc,mu,q,c,alpha,beta,ecc)
% This function calculates the velocity profile between two eccentric
% cylinders based on the exact solution provided in:
% DOI: 10.1002/aic.690110319
%
% inputs:
%   x,y     cartesian cordinates of the points to calculate velocity at
%   rv      radius of the blood vessel
%   rc      radius of the catheter
%   mu      fluid (e.g. blood) viscosity
%   q       fluid (e.g. blood) flowrate
%   c       a constant
%   alpha   line of constant etta in the bipolar system which represents a
%           circle in the cartesian system
%   beta    line of constant etta in the bipolar system which represents a
%           circle in the cartesian system
%   ecc     eccentricity, i.e. distance between the center of catheter and
%           center of blood vessel
%
% outputs:
%   v       velocity in z direction at the points specified with x,y
%
%
% Author: Amirtaha Taebi
% University of California Davis
% Summer 2020
%
% Reference
% Please cite the following manuscript:
%
%

% Body
% parameter definition


etta = 0.5 * log((y.^2+(x+c).^2)./(y.^2+(x-c).^2));
xi = atan(2*y*c./(x.^2+y.^2-c^2));


F = (alpha*coth(beta)-beta*coth(alpha))/(2*(alpha-beta));
E = (coth(alpha)-coth(beta))/(2*(alpha-beta));

syms n
s1 = double(symsum((((coth(alpha)-coth(beta))/(exp(2*n*alpha)-exp(2*n*beta))) * exp(n*etta) + ...
    (((exp(2*n*alpha)*coth(beta)-exp(2*n*beta)*coth(alpha))/(exp(2*n*alpha)-exp(2*n*beta))) - ...
    coth(etta)) .* exp(-n*etta)) .* ...
    cos(n*xi), n, 1, inf));

% non-dimensional velocity
u = F + E*etta - 0.5*coth(etta) + s1;


f = (rv^2-rc^2+ecc^2)/(2*ecc);
M = sqrt(f^2-rv^2);

s2 = double(symsum((n*exp(-n*(alpha+beta)))/(sinh(n*alpha-n*beta)), ...
    n, 1, inf));

% pressure gradient
delP = (8*mu*q/pi) / (rv^4-rc^4-((4*ecc^2*M^2)/(alpha-beta))-8*ecc^2*M^2*s2);

% dimensional velocity
v = u' * c^2 * delP / mu;
end