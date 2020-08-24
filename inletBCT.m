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

%% reading the inlet mesh geometry data
directory = 'E:\Project A\SimvascularDevelopment\example\';
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

theta = atand(n(2)/n(1));
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
newWall = double(cwTy(beta))*double(cwTz(theta))*(wall(:,5:7)'-ctr');
% newInlet = (inlet(:,5:7)'-ctr');
newInlet = double(cwTy(beta))*double(cwTz(theta))*(inlet(:,5:7)'-ctr');

figure(2)
scatter3(newWall(1,:),newWall(2,:),newWall(3,:))
hold on
scatter3(newInlet(1,:),newInlet(2,:),newInlet(3,:),'*k')


% force-projecting all the z values to zero
newWall(3,:) = 0;
newInlet(3,:) = 0;

figure(3)
scatter3(newWall(1,:),newWall(2,:),newWall(3,:))
hold on
scatter3(newInlet(1,:),newInlet(2,:),newInlet(3,:),'*k')
view(2)
