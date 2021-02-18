function h = catheter(x,y,r)
% This function plots a circle (e.g. the catheter) with a radius of r
% centered at (x,y).
%
% inputs:
%   x,y     center of the circle
%   r       radius of the circle

% Body
hold on
theta = 0:pi/50:2*pi;
xunit = r * cos(theta) + x;
yunit = r * sin(theta) + y;
h = plot(xunit, yunit,'k','LineWidth',1.5);
hold off
end