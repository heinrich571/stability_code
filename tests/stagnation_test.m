%% Fresh start

% close all
clearvars
clc


%% Generate coordinates

Ny = 50;
Nx = 2*Ny;

y_limit = 10;
x_limit = 2*y_limit;


%% Compute a simple stagnation point flow

a = 2;
x = linspace(-x_limit,x_limit,Nx+1);
y = linspace(0,y_limit,Ny+1);
[x, y] = meshgrid(x, y);
u =  a*x;
v = -a*y;


%% Plot streamlines

figure('Name','Stagnation point streamlines','NumberTitle','off')
streamline(x,y,u,v,x(1,:),10*ones(size(x(1,:))))
xlim(x_limit*[-1 1])
ylim(y_limit*[0 1])
