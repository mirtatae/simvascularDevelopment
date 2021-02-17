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

%% developments
% TODO: add catheter velocity profile
% TODO: add catheter wall thickness

%%
clear
close all

%% switches
plotOn = 1;
vCat = 1;               % defines the type of the velocity profile inside
                        % the catheter:
                        % 0: zero velocity
                        % 1: parabolic profile
outputFormat = 0;       % output file format:
                        % 0: .dat
                        % 1: .csv
%% parameter definition
mu = 0.004;             % fluid viscosity [use consistent units with other 
                        % parameters]
catR = 0.57/2;             % catheter radius
catT = 0.23/2;             % catheter thickness
ecc = 0.05;              % catheter eccentricity
nl = 201;               % number of time points (entered as "Point Number"
                        % in the "Set Inlet/Outlet BCs>BC Type: Prescribed 
                        % Velocities" in the SimVascular software)
period = 1;             % one period duration (entered as "Period" in the
                        % "Set Inlet/Outlet BCs>BC Type: Prescribed 
                        % Velocities" in the SimVascular software)
catFlow = -330;         % flowrate inside the catheter

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
directory = 'example\';

%% reading flowrate data
filename = 'flowrate.csv';
flowData = xlsread([directory,filename]);

time = 0:period/(nl-1):period;  % time vector

% interpolated flowrate
flowrate = interp1(flowData(:,1),flowData(:,2),time);

%% reading the inlet mesh geometry data
filename = 'Patient8_inlet_coordinates_mesh2.csv';
data = xlsread([directory,filename]);

% finding vessel wall nodes and inlet nodes
ki = find(data(:,2));
inlet = data(ki,:);
inletNodeID = inlet(:,1);

kw = find(~data(:,2));
wall = data(kw,:);
wallNodeID = wall(:,1);

if plotOn == 1
    figure(1)
    subplot(1,2,1)
    scatter3(wall(:,5),wall(:,6),wall(:,7))
    xlabel('X [mm]')
    ylabel('Y [mm]')
    zlabel('Z [mm]')
end

%% normal of the inlet plane
% finding the normal vector of the inlet cross-section

% finding the distance between wall nodes
dist = sqrt(sum((wall(1,5:7) - wall(:,5:7)).^2,2));
[~,i] = sort(dist);

% approximating a center point for the inlet cross-section
ctr = (wall(1,5:7) + wall(i(end),5:7))/2;

% normal vector to the inlet cross-section (in the original coordinate
% system)
n = cross(wall(1,5:7)-wall(i(end),5:7),wall(1,5:7)-wall(i(end-1),5:7));
n = n/norm(n);

gama = atand(n(2)/n(1));
beta = abs(atand(sqrt(n(1)^2+n(2)^2)/n(3)));
% TODO: check the sign of beta for different input files.

if plotOn == 1
    figure(1)
    subplot(1,2,1)
    hold on
    scatter3(wall(1,5),wall(1,6),wall(1,7),'filled','MarkerFaceColor','g')
    scatter3(wall(i(end),5),wall(i(end),6),wall(i(end),7),'filled','MarkerFaceColor','r')
    line([wall(1,5),wall(i(end),5)],[wall(1,6),wall(i(end),6)],[wall(1,7),wall(i(end),7)])
    scatter3(ctr(1),ctr(2),ctr(3),'xk')
    quiver3(ctr(1),ctr(2),ctr(3),n(1),n(2),n(3));
    scatter3(inlet(:,5),inlet(:,6),inlet(:,7),'*k')
    axis equal
end


%%
newWall = double(cwTy(beta))*double(cwTz(gama))*(wall(:,5:7)'-ctr');
newInlet = double(cwTy(beta))*double(cwTz(gama))*(inlet(:,5:7)'-ctr');
% In this new coordinate system the inlet coss-section is in the XY plane
% and the axial velocity is in the Z direction.

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
    subplot(2,2,1)
    scatter3(newWall(1,:),newWall(2,:),newWall(3,:))
    hold on
    scatter3(newInlet(1,:),newInlet(2,:),newInlet(3,:),'*k')
    view(2)
    axis equal
    xlabel('x [mm]')
    ylabel('y [mm]')
end

newWallR = sqrt(newWall(1,:).^2 + newWall(2,:).^2);
vesR = max(newWallR);

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
[vesCtr,catCtr,c,alfa,betta] = centers(catR+catT,vesR,ecc);

% translating the points using the vessel center point
newWall = newWall + [vesCtr;0;0];
newInlet = newInlet + [vesCtr;0;0];

% sorting the wall points based on their angle (just for plotting)
[thetaSort,iW] = sort(theta);
newWallSort = newWall(:,iW);


% finding the inlet points inside and outside the catheter
k = find(abs(sqrt((newInlet(1,:)-catCtr).^2+(newInlet(2,:)).^2))>catR+catT);
inletOutCat = newInlet(:,k);
inletInCat = newInlet;
inletInCat(:,k) = [];

inletOutCatID = inletNodeID(k);
inletInCatID = inletNodeID;
inletInCatID(k) = [];

if plotOn == 1
    figure(3)
    subplot(2,2,2)
    line([newWallSort(1,:),newWallSort(1,1)],[newWallSort(2,:),newWallSort(2,1)],'Color','k','LineWidth',1)
    catheter(catCtr,0,catR);
    catheter(catCtr,0,catR+catT);
    axis equal
    hold on
    scatter(vesCtr,0,'xk')
    scatter(catCtr,0,'xb')
    scatter3(inletOutCat(1,:),inletOutCat(2,:),inletOutCat(3,:),'*k')
    scatter3(inletInCat(1,:),inletInCat(2,:),inletInCat(3,:),'.g')
    catheter(vesCtr,0,vesR);
    title(['Vessel center = ',sprintf('%0.2f',vesCtr)])
    xlabel('x [mm]')
end

%% velocity profile
% velocity profile outside of the catheter
v = velEccCylinders(inletOutCat(1,:),inletOutCat(2,:),vesR,catR+catT,mu,flowrate,c,alfa,betta,ecc);

vv = zeros(size(inletOutCat,1),size(inletOutCat,2),nl);
% Now, the velocity profile should be rotated to be alighned to the normal
% vector (n) of the inlet cross-section in the original coordinate system.
for i = 1:nl
    vv(:,:,i) = double(ccTz(gama))*double(ccTy(beta))* ...
        ([inletOutCat(1,:) - vesCtr;inletOutCat(2,:);v(:,i)']);
end

% zero wall velocity
vWall = zeros(length(newWall),nl);

% velocity profile inside the catheter
rCat = sqrt((inletInCat(1,:)-catCtr).^2+inletInCat(2,:).^2);
if vCat == 0
    vInletInCat = zeros(length(inletInCat),1);
elseif vCat == 1
    vInletInCat = transpose(2*catFlow*(1-(rCat/catR).^2)/(pi*catR^2));
    
    % set the velocity on the catheter wall to zero
    vInletInCat(vInletInCat > 0) = 0;
end
vInletInCat = repelem(vInletInCat,1,nl);

vvInletInCat = zeros(size(inletInCat,1),size(inletInCat,2),nl);
% Now, the velocity profile should be rotated to be alighned to the normal
% vector (n) of the inlet cross-section in the original coordinate system.
for i = 1:nl
    vvInletInCat(:,:,i) = double(ccTz(gama))*double(ccTy(beta))* ...
        ([inletInCat(1,:) - vesCtr;inletInCat(2,:);vInletInCat(:,i)']);
end

if plotOn == 1
    figure(3)
    subplot(2,2,3)
    quiver3([inletOutCat(1,:),inletInCat(1,:)],...
        [inletOutCat(2,:),inletInCat(2,:)],...
        [inletOutCat(3,:),inletInCat(3,:)],...
        zeros(1,length(inletOutCat(1,:))+length(inletInCat(1,:))),...
        zeros(1,length(inletOutCat(1,:))+length(inletInCat(1,:))),...
        [v(:,50)',vInletInCat(:,50)'])
    axis equal
    hold on
    scatter3(inletOutCat(1,:),inletOutCat(2,:),inletOutCat(3,:),'*k')
    scatter3(inletInCat(1,:),inletInCat(2,:),inletInCat(3,:),'.g')
    hold off
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('v [m/s]')
    title('Velocity profile')
    
    subplot(2,2,4)
    [xi,yi] = meshgrid(min(newWall(1,:)):0.01:max(newWall(1,:)),...
        min(newWall(2,:)):0.01:max(newWall(2,:)));
    zi = griddata([inletOutCat(1,:),newWall(1,:),inletInCat(1,:)],...
        [inletOutCat(2,:),newWall(2,:),inletInCat(2,:)],...
        [-v(:,50)',vWall(:,50)',-vInletInCat(:,50)'],xi,yi);
    contourf(xi,yi,zi,'k','ShowText','on','LabelSpacing',400)
    axis equal
    colormap jet
    xlabel('x [mm]')
    ylabel('y [mm]')
    
    %{
    figure
    surf(xi,yi,zi,'EdgeColor','none')
    colormap jet
    caxis([0 700])
    c = colorbar;
    c.Label.String = 'Velocity [mm/s]';
    view(2)
    axis equal
    grid off
    xlabel('x [mm]')
    ylabel('y [mm]')
    %}
    
    figure(1)
    subplot(1,2,2)
    scatter3(inlet(k,5),inlet(k,6),inlet(k,7),'*k')
    hold on
    quiver3(inlet(k,5),inlet(k,6),inlet(k,7),vv(1,:,50)',vv(2,:,50)',vv(3,:,50)')
    scatter3(wall(:,5),wall(:,6),wall(:,7))
    xlabel('X [mm]')
    ylabel('Y [mm]')
    zlabel('Z [mm]')
    axis equal
    
    figure(4)
    line([newWallSort(1,:),newWallSort(1,1)],[newWallSort(2,:),newWallSort(2,1)],'Color','k','LineWidth',1.5)
    catheter(catCtr,0,catR);
    catheter(catCtr,0,catR+catT);
    axis equal
    hold on
    scatter(vesCtr,0,'xk')
    scatter(catCtr,0,'xb')
    scatter3(inletOutCat(1,:),inletOutCat(2,:),inletOutCat(3,:),'k','filled')
    scatter3(inletInCat(1,:),inletInCat(2,:),inletInCat(3,:),'g','filled')
    catheter(vesCtr,0,vesR);
    set(gca,'Visible','off')
    
    
    figure(5)
    [xi,yi] = meshgrid(min(newWall(1,:)):0.01:max(newWall(1,:)),...
        min(newWall(2,:)):0.01:max(newWall(2,:)));
    zi = griddata([inletOutCat(1,:),newWall(1,:),inletInCat(1,:)],...
        [inletOutCat(2,:),newWall(2,:),inletInCat(2,:)],...
        [-v(:,50)',vWall(:,50)',-vInletInCat(:,50)'],xi,yi);
    contourf(xi,yi,zi,[0 100 200 300 400 500 600],'k')
    axis equal
    colormap(jet(5))
    set(gca,'Visible','off')
    colorbar
    caxis([0 500])
    hold on
    scatter3(inletOutCat(1,:),inletOutCat(2,:),inletOutCat(3,:),'k','filled')
    scatter3(inletInCat(1,:),inletInCat(2,:),inletInCat(3,:),'g','filled')
    catheter(catCtr,0,catR);
    catheter(catCtr,0,catR+catT);
    
    fig4name = ['cord_fr',num2str(10*catR),'t',num2str(10*catT),'e',num2str(10*ecc)];
    fig5name = ['velo_catV_fr',num2str(10*catR),'t',num2str(10*catT),'e',num2str(10*ecc)];
%     print('-f4',fig4name,'-djpeg','-r300')
%     print('-f5',fig5name,'-djpeg','-r300')
end


%% saving bct.dat file for simvascular simulation

if outputFormat == 0
    dlmwrite('bct.dat',[length(wall)+length(inlet),nl],' ');
elseif outputFormat == 1
    xlswrite('bct.csv',[length(wall)+length(inlet),nl],'A1:B1');
end

% velocity of the nodes inside the catheter
catCoords = inlet(:,5:7);
catCoords(k,:) = [];
temp2 = zeros(nl,4);
for i = 1:length(catCoords)
    temp1(1,:) = [catCoords(i,1),catCoords(i,2),catCoords(i,3),...
        nl, inletInCatID(i)];
    for j = 1:nl
        temp2(j,:) = [vvInletInCat(1,i,j),vvInletInCat(2,i,j),vvInletInCat(3,i,j),...
        time(j)];
    end
    
    if outputFormat == 0
        dlmwrite('bct.dat',temp1,'delimiter',' ','-append')
        dlmwrite('bct.dat',temp2,'delimiter',' ','-append')
    elseif outputFormat == 1
        row1 = i+(i-1)*nl+1;
        dataRange1 = ['A',num2str(row1),':E',num2str(row1)];
        dataRange2 = ['A',num2str(row1+1),':D',num2str(row1+nl)];
        xlswrite('bct.csv',temp1(1,:),dataRange1);
        xlswrite('bct.csv',temp2(1:nl,:),dataRange2);
    end
end

% velocity of the nodes outside of the catheter
outCatCoords = inlet(k,5:7);
for i = 1:length(outCatCoords)
    temp1(1,:) = [outCatCoords(i,1),outCatCoords(i,2),outCatCoords(i,3),...
        nl, inletOutCatID(i)];
    for j = 1:nl
        temp2(j,:) = [vv(1,i,j),vv(2,i,j),vv(3,i,j),...
        time(j)];
    end
    
    if outputFormat == 0
        dlmwrite('bct.dat',temp1,'delimiter',' ','-append')
        dlmwrite('bct.dat',temp2,'delimiter',' ','-append')
    elseif outputFormat == 1
        row1 = i+(i-1)*nl+(nl+1)*length(catCoords)+1;
        dataRange1 = ['A',num2str(row1),':E',num2str(row1)];
        dataRange2 = ['A',num2str(row1+1),':D',num2str(row1+nl+1)];
        xlswrite('bct.csv',temp1(1,:),dataRange1);
        xlswrite('bct.csv',temp2(1:nl,:),dataRange2);
    end
end

% velocity of the wall nodes
for i = 1:length(wall)
    temp1(1,:) = [wall(i,5),wall(i,6),wall(i,7),...
        nl, wall(i,1)];
    temp2 = zeros(nl,3);
    temp2(:,4) = time;
    
    if outputFormat == 0
        dlmwrite('bct.dat',temp1,'delimiter',' ','-append')
        dlmwrite('bct.dat',temp2,'delimiter',' ','-append')
    elseif outputFormat == 1
        row1 = i+(i-1)*nl+(nl+1)*length(inlet)+1;
        dataRange1 = ['A',num2str(row1),':E',num2str(row1)];
        dataRange2 = ['A',num2str(row1+1),':D',num2str(row1+nl+1)];
        xlswrite('bct.csv',temp1(1,:),dataRange1);
        xlswrite('bct.csv',temp2(1:nl,:),dataRange2);
    end
end

%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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