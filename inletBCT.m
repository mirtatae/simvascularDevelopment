% Code Description
%
% Important parameters:
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
%

%%
clear
close all

%% switches
plotOn = 1;

%% parameter definition
mu = 0.004;             % fluid viscosity [use consistent units with other 
                        % parameters]
catR = 1.5;             % catheter radius
ecc = 0.2;              % catheter eccentricity
nl = 201;               % number of time points (entered as "Point Number"
                        % in the "Set Inlet/Outlet BCs>BC Type: Prescribed 
                        % Velocities" in the SimVascular software)
period = 1;             % one period duration (entered as "Period" in the
                        % "Set Inlet/Outlet BCs>BC Type: Prescribed 
                        % Velocities" in the SimVascular software)

%% rotation matrices
syms ang integer

% counter-clock wise rotation around z and y
ccTz(ang) = [cosd(ang) -sind(ang) 0; ...
    sind(ang) cosd(ang) 0; ...
    0 0 1];
ccTy(ang) = [cosd(ang) 0 sind(ang); ...
    0 1 0; ...
    -sind(ang) 0 cosd(ang)];

% clock wise rotation around z and y
cwTz(ang) = [cosd(ang) sind(ang) 0; ...
    -sind(ang) cosd(ang) 0; ...
    0 0 1];
cwTy(ang) = [cosd(ang) 0 -sind(ang); ...
    0 1 0; ...
    sind(ang) 0 cosd(ang)];

%% setting the directory
directory = 'E:\Project A\SimvascularDevelopment\simvascularDevelopment\example\';

%% reading flowrate data
filename = 'flowrate.csv';
flowData = xlsread([directory,filename]);

time = 0:period/(nl-1):period;  % time vector

% interpolated flowrate
flowrate = interp1(flowData(:,1),flowData(:,2),time);

%% reading the inlet mesh geometry data
filename = 'inlet_coordinates.csv';
data = xlsread([directory,filename]);

% finding vessel wall nodes and inlet nodes
ki = find(data(:,2));
inlet = data(ki,:);

kw = find(~data(:,2));
wall = data(kw,:);

if plotOn == 1
    figure(1)
    scatter3(wall(:,5),wall(:,6),wall(:,7))
    xlabel('x')
    ylabel('y')
    zlabel('z')
end

%% normal of the inlet plane
% finding the normal vector of the inlet cross-section

% finding the distance between wall nodes
dist = sqrt(sum((wall(1,5:7) - wall(:,5:7)).^2,2));
[b,i] = sort(dist);

% approximating a center point for the inlet cross-section
ctr = (wall(1,5:7) + wall(i(end),5:7))/2;

% normal vector to the inlet cross-section
n = cross(wall(1,5:7)-wall(i(end),5:7),wall(1,5:7)-wall(i(end-1),5:7));
n = n/norm(n);

gama = atand(n(2)/n(1));
beta = atand(sqrt(n(1)^2+n(2)^2)/n(3));
d = sqrt(sum(ctr.^2));

if plotOn == 1
    figure(1)
    hold on
    scatter3(wall(1,5),wall(1,6),wall(1,7),'filled','MarkerFaceColor','g')
    scatter3(wall(i(end),5),wall(i(end),6),wall(i(end),7),'filled','MarkerFaceColor','r')
    line([wall(1,5),wall(i(end),5)],[wall(1,6),wall(i(end),6)],[wall(1,7),wall(i(end),7)])
    scatter3(ctr(1),ctr(2),ctr(3),'xk')
    quiver3(ctr(1),ctr(2),ctr(3),n(1),n(2),n(3));
    scatter3(inlet(:,5),inlet(:,6),inlet(:,7),'*k')
end


%%
% newWall = (wall(:,5:7)'-ctr');
newWall = double(cwTy(beta))*double(cwTz(gama))*(wall(:,5:7)'-ctr');
% newInlet = (inlet(:,5:7)'-ctr');
newInlet = double(cwTy(beta))*double(cwTz(gama))*(inlet(:,5:7)'-ctr');

if plotOn == 1
    figure(2)
    scatter3(newWall(1,:),newWall(2,:),newWall(3,:))
    hold on
    scatter3(newInlet(1,:),newInlet(2,:),newInlet(3,:),'*k')
end

% force-projecting all the z values to zero
newWall(3,:) = 0;
newInlet(3,:) = 0;

if plotOn == 1
    figure(3)
    subplot(1,2,1)
    scatter3(newWall(1,:),newWall(2,:),newWall(3,:))
    hold on
    scatter3(newInlet(1,:),newInlet(2,:),newInlet(3,:),'*k')
    view(2)
    axis equal
end

newWallR = sqrt(newWall(1,:).^2 + newWall(2,:).^2);
vesR = max(newWallR);
% vesR = mean(newWallR);

if vesR - catR < ecc
    error('Eccentricity must be smaller than the difference between vessel radius and catheter radius!')
end
%%

theta = atan2d(newWall(2,:),newWall(1,:));
for i = 1:length(newWall)
    if newWall(2,i) < 0
        theta(i) = theta(i) + 360;
    end
end

% finding the catheter and vessel center for the bipolar transformation
[vesCtr,catCtr,c,alfa,betta] = centers(catR,vesR,ecc);

% translating the points using the vessel center point
newWall = newWall + [vesCtr;0;0];
newInlet = newInlet + [vesCtr;0;0];

% sorting the wall points based on their angle
[thetaSort,iW] = sort(theta);
newWallSort = newWall(:,iW);

% testP = newInlet(1:2,1);
testP = [8.786 1.058];

% finding the inlet points inside and outside the catheter
k = find(abs(sqrt((newInlet(1,:)-catCtr).^2+(newInlet(2,:)).^2))>catR);
inletOutCat = newInlet(:,k);
inletInCat = newInlet;
inletInCat(:,k) = [];

if plotOn == 1
    figure(3)
    subplot(1,2,2)
    line([newWallSort(1,:),newWallSort(1,1)],[newWallSort(2,:),newWallSort(2,1)],'Color','k','LineWidth',1)
    catheter(catCtr,0,catR);
    axis equal
    hold on
    scatter(vesCtr,0,'xk')
    scatter(catCtr,0,'xb')
    scatter3(inletOutCat(1,:),inletOutCat(2,:),inletOutCat(3,:),'*k')
    scatter3(inletInCat(1,:),inletInCat(2,:),inletInCat(3,:),'.g')
    scatter(testP(1),testP(2),'or')
    catheter(vesCtr,0,vesR)
end

%% velocity profile
v = velEccCylinders(testP(1),testP(2),vesR,catR,mu,flowrate(100),c,alfa,betta,ecc);
disp(['velocity = ',sprintf('%0.2f',v),' mm/s']);

%% saving bct.dat file for simvascular simulation


%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = catheter(x,y,r)
% This function plots the catheter with a radius of r centered at (x,y).

% Body
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
end

function [cVes,cCat,c,alpha,beta] = centers(rc,rv,ecc)
% function description
%
% inputs:
%
% outputs:
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
gama = rc/rv;
phi = ecc/(rv-rc);

alpha = acosh((gama*(1+phi^2)+(1-phi^2))/(2*phi*gama));
beta = acosh((gama*(1-phi^2)+(1+phi^2))/(2*phi));

c = rc * sinh(alpha);

cCat = c * coth(alpha);
cVes = c * coth(beta);
end

function v = velEccCylinders(x,y,rv,rc,mu,q,c,alpha,beta,ecc)
% This function calculates the velocity profile between two eccentric
% cylinders based on the exact solution provided in:
% DOI: 10.1002/aic.690110319
%
% inputs:
%
% outputs:
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


etta = 0.5 * log((y^2+(x+c)^2)/(y^2+(x-c)^2));
xi = atan(2*y*c/(x^2+y^2-c^2));


F = (alpha*coth(beta)-beta*coth(alpha))/(2*(alpha-beta));
E = (coth(alpha)-coth(beta))/(2*(alpha-beta));

syms n
s1 = double(symsum((((coth(alpha)-coth(beta))/(exp(2*n*alpha)-exp(2*n*beta))) * exp(n*etta) + ...
    (((exp(2*n*alpha)*coth(beta)-exp(2*n*beta)*coth(alpha))/(exp(2*n*alpha)-exp(2*n*beta))) - ...
    coth(etta)) * exp(-n*etta)) * ...
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
v = u * c^2 * delP / mu;

end