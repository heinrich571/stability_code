%% Fresh start

close all
clearvars
clc


%% Define problem parameters

% Reynolds numbers
Re_delta_vec = logspace(1, 3, 3);

% Grid parameters
X_Limit  = 1;
Y_Limit  = 10;
Y_Median = 0.3;
Nx       = 100;
Ny       = Nx/2;

% Definitions
Definitions.initguess            = [1.22 1.24];
Definitions.maxIterations        = 1e2;
Definitions.convergenceTolerance = 1e-6;
Definitions.interval             = [0 Y_Limit];

Base_Flow = get_base_flow(Definitions);


%% Re-normalize the length coordinate

eta  = Base_Flow.eta;
phi  = Base_Flow.phi;
dphi = Base_Flow.dphi;

y_vec = (1./Re_delta_vec).*eta;
y_vec(end+1,:) = y_vec(end,1)+eps;
x_vec = linspace(-X_Limit, X_Limit, Nx + 1);

u = zeros([length(Re_delta_vec) length(x_vec) length(y_vec)]);
v = zeros([length(Re_delta_vec) length(x_vec) length(y_vec)]);
for i = 1:length(Re_delta_vec)
    for j = 1:length(x_vec)
        x = x_vec(j);
        for k = 1:length(y_vec)
            y = y_vec(k,i);
            ind = k;
            if k == length(y_vec)
                ind = k - 1;
            end
            u(i,j,k) = x*dphi(ind);
            v(i,j,k) = -phi(ind);
            V(i,j,k) = sqrt(u(i,j,k)^2 + v(i,j,k)^2);
        end
    end
end


%% Plot results

i = 1;
[X,Y] = meshgrid(x_vec, y_vec(:,i));
[startX, startY] = meshgrid(-X_Limit:0.1:X_Limit, Y_Limit);
streamslice(X, Y, 2000*reshape(u(i,:,:), [length(x_vec) length(y_vec)])', 2000*reshape(v(i,:,:), [length(x_vec) length(y_vec)])',2);

% u velocity component
for i = 1:length(Re_delta_vec)
    figure('Name', 'u map', 'NumberTitle', 'off')
    contourf(x_vec, y_vec(:,i), abs(reshape(u(i,:,:), [length(x_vec) length(y_vec)])'), 50, 'EdgeColor', 'none')
    xlabel('x/\delta')
    ylabel('y/\delta')
    drawnow
end

% v velocity component
for i = 1:length(Re_delta_vec)
    figure('Name', 'u map', 'NumberTitle', 'off')
    contourf(x_vec, y_vec(:,i), abs(reshape(v(i,:,:), [length(x_vec) length(y_vec)])'), 50, 'EdgeColor', 'none')
    xlabel('x/\delta')
    ylabel('y/\delta')
    drawnow
end

% velocity magnitude
for i = 1:length(Re_delta_vec)
    figure('Name', 'u map', 'NumberTitle', 'off')
    contourf(x_vec, y_vec(:,i), abs(reshape(V(i,:,:), [length(x_vec) length(y_vec)])'), 50, 'EdgeColor', 'none')
    xlabel('x/\delta')
    ylabel('y/\delta')
    drawnow
end

















