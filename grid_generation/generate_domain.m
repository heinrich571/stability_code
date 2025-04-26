function Domain = generate_domain(Problem)

% Expand input structure for convenient code
nx      = Problem.Domain.Nx;
ny      = Problem.Domain.Ny;
xlimit  = Problem.Domain.X_Limit;
ylimit  = Problem.Domain.Y_Limit;
ymedian = Problem.Domain.Y_Median;

% Generate the Chebyshev interval and derivatives
[Dx_cheb, xhat] = cheb(nx);
[Dy_cheb, yhat] = cheb(ny);
xhat(abs(xhat) < 1e-10) = 0;
yhat(abs(yhat) < 1e-10) = 0;
Ix = eye(nx + 1);
Iy = eye(ny + 1);

% x mapping
% f      = 1;
% e      = 1 - f;
% xvec   = xlimit * (e * xhat.^3 + f * xhat);
% dxi_dx = 1 ./ (xlimit * (3 * e * xhat.^2 + f));

s      = 4;
xvec   = xlimit * sinh(s*xhat) / sinh(s);
dxi_dx = sinh(s) / (s * xlimit) * 1 ./ sqrt(1 + sinh(s)^2 * (xvec / xlimit).^2);

Dx_physical_domain = transform_to_physical_domain(dxi_dx, Dx_cheb);
D2x_physical_domain = Dx_physical_domain * Dx_physical_domain;
Dx = kron(Dx_physical_domain, Iy);
D2x = kron(D2x_physical_domain, Iy);

% y mapping

% yvec   = ylimit * yhat;
% dxi_dy = 1 / ylimit;

% ay     = ymedian * ylimit / (ylimit - 2 * ymedian);
% by     = 1 + (2 * ay / ylimit);
% yvec   = ay * (1 + yhat) ./ (by - yhat);
% dxi_dy = ay * (1 + by) ./ ((yvec + ay).^2);

yvec   = ymedian * ylimit * (1 + yhat) ./ (ylimit - yhat * (ylimit - 2 * ymedian));
dxi_dy = 2 * ymedian * ylimit * (ylimit - ymedian) ./ ((yvec * (ylimit - 2 * ymedian) + ymedian * ylimit).^2);

Dy_physical_domain = transform_to_physical_domain(dxi_dy, Dy_cheb);
D2y_physical_domain = Dy_physical_domain * Dy_physical_domain;
Dy = kron(Ix, Dy_physical_domain);
D2y = kron(Ix, D2y_physical_domain);

% Generate meshgrid of x and y grid points
[mat_X, mat_Y] = meshgrid(xvec, yvec);

% Create the output structure
Domain.Dx    = Dx   ; % 1st derivative w.r.t. x matrix
Domain.Dy    = Dy   ; % 1st derivative w.r.t. y matrix
Domain.D2x   = D2x  ; % 2nd derivative w.r.t. x matrix
Domain.D2y   = D2y  ; % 2nd derivative w.r.t. y matrix
Domain.vec_X = xvec; % Vector of the x stations
Domain.vec_Y = yvec; % Vector of the y stations
Domain.mat_X = mat_X; % Matrix of the x stations
Domain.mat_Y = mat_Y; % Matrix of the y stations

end


% Supporting functions
function physderiv = transform_to_physical_domain(strech, chebderiv)

physderiv = diag(strech) * chebderiv;

end
