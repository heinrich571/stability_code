function test_solution_symmetry(Base_Flow, Solution, si, evi, Options)

% Handle special inputs
if nargin < 4
    Options.Solution_Index = length(Solution);
    Options.X_Limit = Solution(Options.Solution_Index).Domain.vec_X(1);
    Options.Y_Limit = 20;
end
x_limit = Options.X_Limit;
y_limit = Options.Y_Limit;

% Parse domain information
Nx = length(Solution(si).Domain.vec_X);
Ny = length(Solution(si).Domain.vec_Y);
Dx = Solution(si).Domain.Dx;
Dy = Solution(si).Domain.Dy;
D2x = Solution(si).Domain.D2x;
D2y = Solution(si).Domain.D2y;
xmat = Solution(si).Domain.mat_X;
ymat = Solution(si).Domain.mat_Y;
beta = Solution.Physics.Beta;

% Parse base flow information
phi = repmat(flip(Base_Flow.phi) , [1 Nx]);
dphi = repmat(flip(Base_Flow.dphi) , [1 Nx]);
ddphi = repmat(flip(Base_Flow.ddphi), [1 Nx]);
Z = zeros(Nx * Ny);
I = eye(Nx * Ny);
iI = 1i * I;
Ux = diag(dphi(:), 0);
Uy = xmat .* ddphi;
Uy = diag(Uy(:), 0);
Vy = diag(-dphi(:), 0);
U = xmat .* dphi;
U = diag(U(:), 0);
V = diag(-phi(:), 0);
Ux = diag(dphi(:), 0);
Uy = xmat .* ddphi;
Uy = diag(Uy(:), 0);
Vx = Z;
Vy = diag(-dphi(:), 0);
L = (D2x + D2y - beta^2*I) - U*Dx - V*Dy;


% Parse solution information
omega = Solution(si).Eigenvalues(evi);
u = reshape(Solution(si).Eigenfunctions.u(:,evi), Ny, Nx, length(evi));
v = reshape(Solution(si).Eigenfunctions.v(:,evi), Ny, Nx, length(evi));
w = reshape(Solution(si).Eigenfunctions.w(:,evi), Ny, Nx, length(evi));
p = reshape(Solution(si).Eigenfunctions.p(:,evi), Ny, Nx, length(evi));

% Calculate differences between left & right sides of the domain
jm = (Nx - 1) / 2 + 1;

min_u = min(real(u), [], 'all');
max_u = max(real(u), [], 'all');

if sign(min_u) ~= sign(max_u) && sign(min_u) ~= 0
    mode_type = 'Anti_Symmetric';
else
    mode_type = 'Symmetric';
end
for j = 1 : (jm - 1)
    jR = j;
    jL = Nx - j + 1;
    switch mode_type
        case 'Anti_Symmetric'
            u_diff(:,j) = u(:,jR) + u(:,jL);
            v_diff(:,j) = v(:,jR) - v(:,jL);
            w_diff(:,j) = w(:,jR) - w(:,jL);
            p_diff(:,j) = p(:,jR) - p(:,jL);
        otherwise
            u_diff(:,j) = u(:,jR) - u(:,jL);
            v_diff(:,j) = v(:,jR) + v(:,jL);
            w_diff(:,j) = w(:,jR) + w(:,jL);
            p_diff(:,j) = p(:,jR) + p(:,jL);
    end
end

u_diff(:,jm) = 0;
v_diff(:,jm) = 0;
w_diff(:,jm) = 0;
p_diff(:,jm) = 0;

u_diff(:,jm + 1 : Nx) = -u_diff(:,flip(1:jm - 1));
v_diff(:,jm + 1 : Nx) = -v_diff(:,flip(1:jm - 1));
w_diff(:,jm + 1 : Nx) = -w_diff(:,flip(1:jm - 1));
p_diff(:,jm + 1 : Nx) = -p_diff(:,flip(1:jm - 1));

% Plot differences
figure('Name', 'velocity differences', 'NumberTitle', 'off')
subplot(1,3,1)
surf(xmat, ymat, abs(u_diff), 'EdgeColor', 'interp', 'FaceColor', 'none')
xlabel('$x$')
ylabel('$y$')
zlabel('$|\Delta u|$')
xlim([-x_limit x_limit])
ylim([0 y_limit])
view(-50, 10)
grid off
subplot(1,3,2)
surf(xmat, ymat, abs(v_diff), 'EdgeColor', 'interp', 'FaceColor', 'none')
xlabel('$x$')
ylabel('$y$')
zlabel('$|\Delta v|$')
xlim([-x_limit x_limit])
ylim([0 y_limit])
view(-50, 10)
grid off
subplot(1,3,3)
surf(xmat, ymat, abs(w_diff), 'EdgeColor', 'interp', 'FaceColor', 'none')
xlabel('$x$')
ylabel('$y$')
zlabel('$|\Delta w|$')
xlim([-x_limit x_limit])
ylim([0 y_limit])
view(-50, 10)
grid off

figure('Name', 'p differences', 'NumberTitle', 'off')
surf(xmat, ymat, abs(p_diff), 'EdgeColor', 'interp', 'FaceColor', 'none')
xlabel('$x$')
ylabel('$y$')
zlabel('$|\Delta p|$')
xlim([-x_limit x_limit])
ylim([0 y_limit])
view(-50, 10)
grid off

% Calculate the differences in equation members
% x-momentum
a11 =  L - Ux; a12 = -Uy    ; a13 = Z        ; a14 = -Dx       ;    b11 = -iI; b12 =  Z ; b13 =  Z ; b14 = Z;

% y-momentum
a21 = -Vx    ; a22 =  L - Vy; a23 = Z        ; a24 = -Dy       ;    b21 =  Z ; b22 = -iI; b23 =  Z ; b24 = Z;

% z-momentum
a31 =  Z     ; a32 =  Z     ; a33 = L        ; a34 = -beta * iI;    b31 =  Z ; b32 =  Z ; b33 = -iI; b34 = Z;

% Continuity
a41 = Dx     ; a42 =  Dy    ; a43 = beta * iI; a44 =  Z        ;    b41 =  Z ; b42 =  Z ; b43 =  Z ; b44 = Z;

for j = 1 : Nx
    a1 = a11 * Solution(si).Eigenfunctions.u(:,evi) + ...
         a12 * Solution(si).Eigenfunctions.v(:,evi) + ...
         a13 * Solution(si).Eigenfunctions.w(:,evi) + ...
         a14 * Solution(si).Eigenfunctions.p(:,evi);
    b1 = b11 * Solution(si).Eigenfunctions.u(:,evi) + ...
         b12 * Solution(si).Eigenfunctions.v(:,evi) + ...
         b13 * Solution(si).Eigenfunctions.w(:,evi) + ...
         b14 * Solution(si).Eigenfunctions.p(:,evi);
    b1 = omega * b1;

    a1_LHS = reshape(a1, Ny, Nx, length(evi));
    b1_RHS = reshape(b1, Ny, Nx, length(evi));
    
    d11(:,j) = a1_LHS(:,j) - b1_RHS(:,j);
end

% d11
figure('Name', 'd11', 'NumberTitle', 'off')
surf(xmat, ymat, abs(d11), 'EdgeColor', 'interp', 'FaceColor', 'none')
xlabel('$x$')
ylabel('$y$')
zlabel('$|\Delta d_{11}|$')
xlim([-x_limit x_limit])
ylim([0 y_limit])
view(-50, 10)
grid off


% Calculate equation members in the domain



end
