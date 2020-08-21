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
%% reading the inlet mesh geometry data
directory = 'E:\Project A\SimvascularDevelopment\example\';
filename = 'inlet_coordinates.csv';

data = xlsread([directory,filename]);

% finding vessel wall nodes and inlet nodes
ki = find(data(:,2));
inlet = data(ki,:);

kw = find(~data(:,2));
wall = data(kw,:);

figure(1)
scatter3(wall(:,5),wall(:,6),wall(:,7))
xlabel('x')
ylabel('y')
zlabel('z')


%% normal of the inlet plane
% finding the normal vector of the inlet cross-section

% finding the distance between wall nodes
dist = sqrt(sum((wall(1,5:7) - wall(:,5:7)).^2,2));
[b,i] = sort(dist);

figure(1)
hold on
scatter3(wall(1,5),wall(1,6),wall(1,7),'filled','MarkerFaceColor','g')
scatter3(wall(i(end),5),wall(i(end),6),wall(i(end),7),'filled','MarkerFaceColor','r')
scatter3(wall(i(end-1),5),wall(i(end-1),6),wall(i(end-1),7),'filled','MarkerFaceColor','r')
line([wall(1,5),wall(i(end),5)],[wall(1,6),wall(i(end),6)],[wall(1,7),wall(i(end),7)])

% approximating a center point for the inlet cross-section
ctr = (wall(1,5:7) + wall(i(end),5:7))/2;

figure(1)
hold on
scatter3(ctr(1),ctr(2),ctr(3),'xk')


c = cross(wall(1,5:7)-wall(i(end),5:7),wall(1,5:7)-wall(i(end-1),5:7));
c = c/norm(c);

figure(1)
hold on
quiver3(ctr(1),ctr(2),ctr(3),c(1),c(2),c(3));