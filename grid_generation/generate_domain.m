function Domain = generate_domain(X_Limit, Y_Limit, Y_Median, Nx, Ny)

% Generate the Chebyshev interval

[Dx_cheb, xhat] = cheb(Nx);
[Dy_cheb, yhat] = cheb(Ny);


% Define the mapping

a_y = Y_Median*Y_Limit/(Y_Limit-2*Y_Median);
b_y = 1+2*a_y/Y_Limit;
vec_Y   = (a_y*(1+yhat)./(b_y-yhat));
% vec_Y = Y_Limit*(yhat + 1)/2;

f = 0.5;
e = 1-f;
vec_X = X_Limit*(e*xhat.^3 + f*xhat);
% vec_X = X_Limit*xhat;

[mat_X, mat_Y] = meshgrid(vec_X, vec_Y);

Dx_physical_domain = diag(1./(X_Limit*(3*e*xhat.^2+f)))*Dx_cheb;
% Dx_physical_domain = X_Limit*Dx_cheb;
Dy_physical_domain = diag(a_y*(1+b_y)./((vec_Y+a_y).^2))*Dy_cheb;
% Dy_physical_domain = Y_Limit/2*Dy_cheb;
D2x_physical_domain = Dx_physical_domain*Dx_physical_domain;
D2y_physical_domain = Dy_physical_domain*Dy_physical_domain;
Ix = eye(Nx+1);
Iy = eye(Ny+1);

Dx = kron(Dx_physical_domain, Iy);
Dy = kron(Ix, Dy_physical_domain);

D2x = kron(D2x_physical_domain, Iy);
D2y = kron(Ix, D2y_physical_domain);


% Create the output structure

Domain.vec_X = vec_X;                                                       % Vector of the x stations
Domain.vec_Y = vec_Y;                                                       % Vector of the y stations
Domain.mat_X = mat_X;                                                       % Matrix of the x stations
Domain.mat_Y = mat_Y;                                                       % Matrix of the y stations
Domain.Dx    = Dx   ;                                                       % 1st derivative w.r.t. x matrix
Domain.Dy    = Dy   ;                                                       % 1st derivative w.r.t. y matrix
Domain.D2x   = D2x  ;                                                       % 2nd derivative w.r.t. x matrix
Domain.D2y   = D2y  ;                                                       % 2nd derivative w.r.t. y matrix

end
