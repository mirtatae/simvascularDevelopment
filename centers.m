function [cVes,cCat,c,alpha,beta] = centers(rc,rv,ecc)
% This function calculates the center of the vessel and the catheter for
% correct transformation in the bipolar system.
%
% inputs:
%   rc      radius of catheter
%   rv      radius of blood vessel
%   ecc     eccentricity, i.e. distance between the center of catheter and
%           center of blood vessel
%
% outputs:
%   cVes    center of the blood vessel
%   cCat    center of the catheter
%   c       a constant
%   alpha   line of constant etta in the bipolar system which represents a
%           circle in the cartesian system
%   beta    line of constant etta in the bipolar system which represents a
%           circle in the cartesian system
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
gama = rc/rv;
phi = ecc/(rv-rc);

alpha = acosh((gama*(1+phi^2)+(1-phi^2))/(2*phi*gama));
beta = acosh((gama*(1-phi^2)+(1+phi^2))/(2*phi));

c = rc * sinh(alpha);

cCat = c * coth(alpha);
cVes = c * coth(beta);
end