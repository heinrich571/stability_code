function [mat_A, mat_B] = create_eigenvalue_matrices(Domain, Base_Flow, beta)

% Initialize matrix elements

initializer = zeros(size(Domain.Dx));
Ny  = size(Domain.mat_X, 1);
Nx  = size(Domain.mat_X, 2);
Z   = initializer;
I   = eye(Ny*Nx);


% Compute matrix terms

phi   = flip(Base_Flow.phi);
dphi  = flip(Base_Flow.dphi);
ddphi = flip(Base_Flow.ddphi);

mat_X = Domain.mat_X;
Dx    = Domain.Dx;
Dy    = Domain.Dy;

mat_phi   = repmat(phi  , [1,Nx]);
mat_dphi  = repmat(dphi , [1,Nx]);
mat_ddphi = repmat(ddphi, [1,Nx]);
U  = mat_X.*mat_dphi;
U  = diag(U(:), 0);
V  = diag(-mat_phi(:), 0);
Ux = diag(mat_dphi(:), 0);
Uy = mat_X.*mat_ddphi;
Uy = diag(Uy(:), 0);
Vx = Z;
Vy = diag(-mat_dphi(:), 0);
L  = (U - Dx)*Dx - (-V + Dy)*Dy + beta^2;


% LHS matrix entries
% continuity
a11 = Dx;
a12 = Dy;
a13 = 1i*beta*I;
a14 = Z;

% x-momentum
a21 = L + Ux;
a22 = Uy;
a23 = Z;
a24 = Dx;

% y-momentum
a31 = Vx;
a32 = L + Vy;
a33 = Z;
a34 = Dy;

% z-momentum
a41 = Z;
a42 = Z;
a43 = L;
a44 = 1i*beta*I;


% RHS matrix entries
% continuity
b11 = Z;
b12 = Z;
b13 = Z;
b14 = Z;

% x-momentum
b21 = 1i*I;
b22 = Z;
b23 = Z;
b24 = Z;

% y-momentum

b31 = Z;
b32 = 1i*I;
b33 = Z;
b34 = Z;

% z-momentum
b41 = Z;
b42 = Z;
b43 = 1i*I;
b44 = Z;


% Build eigenvalue problem A and B matrices

mat_A = [a11 a12 a13 a14
         a21 a22 a23 a24
         a31 a32 a33 a34
         a41 a42 a43 a44];

mat_B = [b11 b12 b13 b14
         b21 b22 b23 b24
         b31 b32 b33 b34
         b41 b42 b43 b44];


% Apply boundary conditions - NOTE THAT HENCEFORTH Nx and Ny are the number
% of points and not the number of sections!

dirichlet_factor     = 200;
neumann_factor       = 300;
linear_extrap_factor = 400;

% No-slip on the wall:
% u(y = 0) = 0
row_inds    = get_eqn_bottom_inds('x momentum', Nx, Ny);
column_inds = get_var_bottom_inds('u', Nx, Ny);
linear_inds = get_linear_indices(mat_A, row_inds, column_inds);

mat_A(row_inds(:),:) = 0;
mat_B(row_inds(:),:) = 0;
mat_A(linear_inds)   = dirichlet_factor;
mat_B(linear_inds)   = 1;


% w(y = 0) = 0
row_inds    = get_eqn_bottom_inds('z momentum', Nx, Ny);
column_inds = get_var_bottom_inds('w', Nx, Ny);
linear_inds = get_linear_indices(mat_A, row_inds, column_inds);

mat_A(row_inds(:),:) = 0;
mat_B(row_inds(:),:) = 0;
mat_A(linear_inds)   = dirichlet_factor;
mat_B(linear_inds)   = 1;


% No-penetration on the wall
% v(y = 0) = 0
row_inds    = get_eqn_bottom_inds('y momentum', Nx, Ny);
column_inds = get_var_bottom_inds('v', Nx, Ny);
linear_inds = get_linear_indices(mat_A, row_inds, column_inds);

mat_A(row_inds(:),:) = 0;
mat_B(row_inds(:),:) = 0;
mat_A(linear_inds)   = dirichlet_factor;
mat_B(linear_inds)   = 1;


% % % % Pressure compatibility condition
% % % row_inds_continuity = get_eqn_bottom_inds('continuity', Nx, Ny);
% % % row_inds_x_momentum = get_eqn_bottom_inds('x momentum', Nx, Ny);
% % % row_inds_y_momentum = get_eqn_bottom_inds('y momentum', Nx, Ny);
% % % 
% % % mat_A(row_inds_continuity(:),:) = mat_A(row_inds_x_momentum(:),:);
% % % mat_B(row_inds_continuity(:),:) = 0;

% Disturbances decay as y --> infinity
% u(y --> inf) = 0
row_inds    = get_eqn_top_inds('x momentum', Nx, Ny);
column_inds = get_var_top_inds('u', Nx, Ny);
linear_inds = get_linear_indices(mat_A, row_inds, column_inds);

mat_A(row_inds(:),:) = 0;
mat_B(row_inds(:),:) = 0;
mat_A(linear_inds)   = dirichlet_factor;
mat_B(linear_inds)   = 1;

% v(y --> inf) = 0
row_inds    = get_eqn_top_inds('y momentum', Nx, Ny);
column_inds = get_var_top_inds('v', Nx, Ny);
linear_inds = get_linear_indices(mat_A, row_inds, column_inds);

mat_A(row_inds(:),:) = 0;
mat_B(row_inds(:),:) = 0;
mat_A(linear_inds)   = dirichlet_factor;
mat_B(linear_inds)   = 1;

% w(y --> inf) = 0
row_inds    = get_eqn_top_inds('z momentum', Nx, Ny);
column_inds = get_var_top_inds('w', Nx, Ny);
linear_inds = get_linear_indices(mat_A, row_inds, column_inds);

mat_A(row_inds(:),:) = 0;
mat_B(row_inds(:),:) = 0;
mat_A(linear_inds)   = dirichlet_factor;
mat_B(linear_inds)   = 1;


% No vertical acceleration at the top of the domain
% p(y --> inf) = 0
row_inds    = get_eqn_top_inds('continuity', Nx, Ny);
column_inds = get_var_top_inds('p', Nx, Ny);
linear_inds = get_linear_indices(mat_A, row_inds, column_inds);

mat_A(row_inds(:),:) = 0;
mat_B(row_inds(:),:) = 0;
mat_A(linear_inds)   = dirichlet_factor;
mat_B(linear_inds)   = 1;


% Linear extrapolation of disturbances at the chordwise directions
% RIGHT BOUNDARY
% % d2u/dx2 = 0
% operator_row_inds = get_opr_right_inds(Nx, Ny);
% row_inds          = get_eqn_right_inds('x momentum', Nx, Ny);
% operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
% row_inds          = row_inds(2:end-1);
% 
% j_u_TR = get_var_top_right_ind('u', Nx, Ny);
% j_u_BL = get_var_bottom_left_ind('u', Nx, Ny);
% column_inds = j_u_TR:j_u_BL;
% 
% mat_A(row_inds,:) = 0;
% mat_B(row_inds,:) = 0;
% mat_A(row_inds,column_inds) = linear_extrap_factor*D2x(operator_row_inds,:);
% mat_B(row_inds,column_inds) = 1*D2x(operator_row_inds,:);
% 
% % d2v/dx2 = 0
% operator_row_inds = get_opr_right_inds(Nx, Ny);
% row_inds          = get_eqn_right_inds('y momentum', Nx, Ny);
% operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
% row_inds          = row_inds(2:end-1);
% 
% j_v_TR = get_var_top_right_ind('v', Nx, Ny);
% j_v_BL = get_var_bottom_left_ind('v', Nx, Ny);
% column_inds = j_v_TR:j_v_BL;
% 
% mat_A(row_inds,:) = 0;
% mat_B(row_inds,:) = 0;
% mat_A(row_inds,column_inds) = linear_extrap_factor*D2x(operator_row_inds,:);
% mat_B(row_inds,column_inds) = 1*D2x(operator_row_inds,:);
% 
% % d2w/dx2 = 0
% operator_row_inds = get_opr_right_inds(Nx, Ny);
% row_inds          = get_eqn_right_inds('z momentum', Nx, Ny);
% operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
% row_inds          = row_inds(2:end-1);
% 
% j_w_TR = get_var_top_right_ind('w', Nx, Ny);
% j_w_BL = get_var_bottom_left_ind('w', Nx, Ny);
% column_inds = j_w_TR:j_w_BL;
% 
% mat_A(row_inds,:) = 0;
% mat_B(row_inds,:) = 0;
% mat_A(row_inds,column_inds) = linear_extrap_factor*D2x(operator_row_inds,:);
% mat_B(row_inds,column_inds) = 1*D2x(operator_row_inds,:);
% 
% % d2p/dx2 = 0 - replacing a continuity equation at the right boundary
% operator_row_inds = get_opr_right_inds(Nx, Ny);
% row_inds          = get_eqn_right_inds('continuity', Nx, Ny);
% operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
% row_inds          = row_inds(2:end-1);
% 
% j_p_TR = get_var_top_right_ind('p', Nx, Ny);
% j_p_BL = get_var_bottom_left_ind('p', Nx, Ny);
% column_inds = j_p_TR:j_p_BL;
% 
% mat_A(row_inds,:) = 0;
% mat_B(row_inds,:) = 0;
% mat_A(row_inds,column_inds) = linear_extrap_factor*Dx(operator_row_inds,:);
% mat_B(row_inds,column_inds) = 1*Dx(operator_row_inds,:);


% LEFT BOUNDARY
% d2u/dx2 = 0
% operator_row_inds = get_opr_left_inds(Nx, Ny);
% row_inds          = get_eqn_left_inds('x momentum', Nx, Ny);
% operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
% row_inds          = row_inds(2:end-1);
% 
% j_u_TR = get_var_top_right_ind('u', Nx, Ny);
% j_u_BL = get_var_bottom_left_ind('u', Nx, Ny);
% column_inds = j_u_TR:j_u_BL;
% 
% mat_A(row_inds,:) = 0;
% mat_B(row_inds,:) = 0;
% mat_A(row_inds,column_inds) = linear_extrap_factor*D2x(operator_row_inds,:);
% mat_B(row_inds,column_inds) = 1*D2x(operator_row_inds,:);
% 
% % d2v/dx2 = 0
% operator_row_inds = get_opr_left_inds(Nx, Ny);
% row_inds          = get_eqn_left_inds('y momentum', Nx, Ny);
% operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
% row_inds          = row_inds(2:end-1);
% 
% j_v_TR = get_var_top_right_ind('v', Nx, Ny);
% j_v_BL = get_var_bottom_left_ind('v', Nx, Ny);
% column_inds = j_v_TR:j_v_BL;
% 
% mat_A(row_inds,:) = 0;
% mat_B(row_inds,:) = 0;
% mat_A(row_inds,column_inds) = linear_extrap_factor*D2x(operator_row_inds,:);
% mat_B(row_inds,column_inds) = 1*D2x(operator_row_inds,:);
% 
% % d2w/dx2 = 0
% operator_row_inds = get_opr_left_inds(Nx, Ny);
% row_inds          = get_eqn_left_inds('z momentum', Nx, Ny);
% operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
% row_inds          = row_inds(2:end-1);
% 
% j_w_TR = get_var_top_right_ind('w', Nx, Ny);
% j_w_BL = get_var_bottom_left_ind('w', Nx, Ny);
% column_inds = j_w_TR:j_w_BL;
% 
% mat_A(row_inds,:) = 0;
% mat_B(row_inds,:) = 0;
% mat_A(row_inds,column_inds) = linear_extrap_factor*D2x(operator_row_inds,:);
% mat_B(row_inds,column_inds) = 1*D2x(operator_row_inds,:);
% 
% % d2p/dx2 = 0 - replacing a continuity equation at the left boundary
% operator_row_inds = get_opr_left_inds(Nx, Ny);
% row_inds          = get_eqn_left_inds('continuity', Nx, Ny);
% operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
% row_inds          = row_inds(2:end-1);
% 
% j_p_TR = get_var_top_right_ind('p', Nx, Ny);
% j_p_BL = get_var_bottom_left_ind('p', Nx, Ny);
% column_inds = j_p_TR:j_p_BL;
% 
% mat_A(row_inds,:) = 0;
% mat_B(row_inds,:) = 0;
% mat_A(row_inds,column_inds) = linear_extrap_factor*Dx(operator_row_inds,:);
% mat_B(row_inds,column_inds) = 1*Dx(operator_row_inds,:);


% Linear extrapolation - right boundary
j_u_right_inds = get_var_right_inds('u', Nx, Ny);
j_v_right_inds = get_var_right_inds('v', Nx, Ny);
j_w_right_inds = get_var_right_inds('w', Nx, Ny);
j_p_right_inds = get_var_right_inds('p', Nx, Ny);
i_xmom_right_inds = get_eqn_right_inds('x momentum', Nx, Ny);
i_ymom_right_inds = get_eqn_right_inds('y momentum', Nx, Ny);
i_zmom_right_inds = get_eqn_right_inds('z momentum', Nx, Ny);
i_cont_right_inds = get_eqn_right_inds('continuity', Nx, Ny);
mat_A(i_xmom_right_inds(2:end-1),:) = 0;
mat_B(i_xmom_right_inds(2:end-1),:) = 0;
mat_A(i_ymom_right_inds(2:end-1),:) = 0;
mat_B(i_ymom_right_inds(2:end-1),:) = 0;
mat_A(i_zmom_right_inds(2:end-1),:) = 0;
mat_B(i_zmom_right_inds(2:end-1),:) = 0;
mat_A(i_cont_right_inds(2:end-1),:) = 0;
mat_B(i_cont_right_inds(2:end-1),:) = 0;
for i = 2:Ny-1
    C = (mat_X(i,1)-mat_X(i,2))/(mat_X(i,3)-mat_X(i,2));
    linear_extrap_opr = [1 C-1 -C];
    ind_shift = Ny*[0 1 2];
    % u
    mat_A(i_xmom_right_inds(i),j_u_right_inds(i) + ind_shift) = linear_extrap_factor*linear_extrap_opr;
    mat_B(i_xmom_right_inds(i),j_u_right_inds(i) + ind_shift) = linear_extrap_opr;
    % v
    mat_A(i_ymom_right_inds(i),j_v_right_inds(i) + ind_shift) = linear_extrap_factor*linear_extrap_opr;
    mat_B(i_ymom_right_inds(i),j_v_right_inds(i) + ind_shift) = linear_extrap_opr;
    % w
    mat_A(i_zmom_right_inds(i),j_w_right_inds(i) + ind_shift) = linear_extrap_factor*linear_extrap_opr;
    mat_B(i_zmom_right_inds(i),j_w_right_inds(i) + ind_shift) = linear_extrap_opr;
    % p
    mat_A(i_cont_right_inds(i),j_p_right_inds(i) + ind_shift) = linear_extrap_factor*linear_extrap_opr;
    mat_B(i_cont_right_inds(i),j_p_right_inds(i) + ind_shift) = linear_extrap_opr;
end

% Linear extrapolation - left boundary
j_u_left_inds = get_var_left_inds('u', Nx, Ny);
j_v_left_inds = get_var_left_inds('v', Nx, Ny);
j_w_left_inds = get_var_left_inds('w', Nx, Ny);
j_p_left_inds = get_var_left_inds('p', Nx, Ny);
i_xmom_left_inds = get_eqn_left_inds('x momentum', Nx, Ny);
i_ymom_left_inds = get_eqn_left_inds('y momentum', Nx, Ny);
i_zmom_left_inds = get_eqn_left_inds('z momentum', Nx, Ny);
i_cont_left_inds = get_eqn_left_inds('continuity', Nx, Ny);
mat_A(i_xmom_left_inds(2:end-1),:) = 0;
mat_B(i_xmom_left_inds(2:end-1),:) = 0;
mat_A(i_ymom_left_inds(2:end-1),:) = 0;
mat_B(i_ymom_left_inds(2:end-1),:) = 0;
mat_A(i_zmom_left_inds(2:end-1),:) = 0;
mat_B(i_zmom_left_inds(2:end-1),:) = 0;
mat_A(i_cont_left_inds(2:end-1),:) = 0;
mat_B(i_cont_left_inds(2:end-1),:) = 0;
for i = 2:Ny-1
    C = (mat_X(i,Nx)-mat_X(i,Nx-1))/(mat_X(i,Nx-2)-mat_X(i,Nx-1));
    linear_extrap_opr = [1 C-1 -C];
    ind_shift = -Ny*[0 1 2];
    % u
    mat_A(i_xmom_left_inds(i),j_u_left_inds(i) + ind_shift) = linear_extrap_factor*linear_extrap_opr;
    mat_B(i_xmom_left_inds(i),j_u_left_inds(i) + ind_shift) = linear_extrap_opr;
    % v
    mat_A(i_ymom_left_inds(i),j_v_left_inds(i) + ind_shift) = linear_extrap_factor*linear_extrap_opr;
    mat_B(i_ymom_left_inds(i),j_v_left_inds(i) + ind_shift) = linear_extrap_opr;
    % w
    mat_A(i_zmom_left_inds(i),j_w_left_inds(i) + ind_shift) = linear_extrap_factor*linear_extrap_opr;
    mat_B(i_zmom_left_inds(i),j_w_left_inds(i) + ind_shift) = linear_extrap_opr;
    % p
    mat_A(i_cont_left_inds(i),j_p_left_inds(i) + ind_shift) = linear_extrap_factor*linear_extrap_opr;
    mat_B(i_cont_left_inds(i),j_p_left_inds(i) + ind_shift) = linear_extrap_opr;
end

end


%%% Supporting functions %%%

function linear_inds = get_linear_indices(mat, row_inds, column_inds)

% DESCRIPTION
%   This function returns the linear indices for the specified matrix based on specified row and column indices.
% INPUT
%   mat             matrix          [matrix]
%   row_inds        row indices     [vector]
%   column_inds     column indices  [vector]
% OUTPUT
%   linear_inds     linear indices corresponding to 'row_inds' and 'column_inds' in the input matrix

row_inds    = row_inds(:);
column_inds = column_inds(:);
linear_inds = sub2ind(size(mat), row_inds, column_inds);

end


function equation_number = eqn2num(equation_name)

% DESCRIPTION
%   This function returns the equation number according to their order in the eigenvalue matrix formulation
% INPUT
%   equation_name       name of the equation                                                                    [char]
% OUTPUT
%   equation_number     number of the equation according to its position in the eigenvalue problem formulation  [scalar]

equation_names  = {'continuity' , 'x momentum' , 'y momentum' , 'z momentum'};
equation_number = find(strcmp(equation_names, equation_name));

if isempty(equation_number)
    error('Equation name is invalid')
end

end


function variable_number = var2num(variable_name)

% DESCRIPTION
%   This function returns the variable number according to their order in the eigenvalue vector formulation
% INPUT
%   variable_name       name of the variable                                                                    [char]
% OUTPUT
%   variable_number     number of the variable according to its position in the eigenvalue problem formulation  [scalar]

variable_names  = {'u' , 'v' , 'w' , 'p'};
variable_number = find(strcmp(variable_names, variable_name));

if isempty(variable_number)
    error('Variable name is invalid')
end

end


function i_eqn_TR = get_eqn_top_right_ind(equation_name, Nx, Ny)

% DESCRIPTION
%   This function returns the index of the specified equation at the top-right corner of the domain in the eigenvalue problem matrices.
% INPUT:
%   equation_name       name of the equation                                                                                        [char]
% OUTPUT:
%   i_eqn_TR            row index of specified equation at the top-right corner of the domain in the eigenvalue problem matrices    [scalar]

equation_number = eqn2num(equation_name);
i_eqn_TR = 1 + Nx*Ny*(equation_number - 1);

end


function j_var_TR = get_var_top_right_ind(variable_name, Nx, Ny)

% DESCRIPTION
%   This function returns the index of the specified variable at the top-right corner of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the variable                                                                                        [char]
% OUTPUT:
%   j_var_TR            column index of specified variable at the top-right corner of the domain in the eigenvalue problem matrices [scalar]

variable_number = var2num(variable_name);
j_var_TR = 1 + Nx*Ny*(variable_number - 1);

end


function i_eqn_BL = get_eqn_bottom_left_ind(equation_name, Nx, Ny)

% DESCRIPTION
%   This function returns the index of the specified equation at the bottom-left corner of the domain in the eigenvalue problem matrices.
% INPUT:
%   equation_name       name of the equation                                                                                        [char]
% OUTPUT:
%   i_eqn_BL            row index of specified equation at the bottom-left corner of the domain in the eigenvalue problem matrices  [scalar]

% equation_number = eqn2num(equation_name);
% i_eqn_BL = Nx*Ny*equation_number;

i_eqn_TR = get_eqn_top_right_ind(equation_name, Nx, Ny);
i_eqn_BL = i_eqn_TR + Nx*Ny - 1;

end


function j_var_BL = get_var_bottom_left_ind(variable_name, Nx, Ny)

% DESCRIPTION
%   This function returns the index of the specified variable at the bottom-left corner of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the variable                                                                                            [char]
% OUTPUT:
%   i_eqn_BL            column index of specified variable at the bottom-left corner of the domain in the eigenvalue problem matrices   [scalar]

% variable_number = var2num(variable_name);
% j_var_BL = Nx*Ny*variable_number;

j_var_TR = get_var_top_right_ind(variable_name, Nx, Ny);
j_var_BL = j_var_TR + Nx*Ny - 1;

end


function i_eqn_BR = get_eqn_bottom_right_ind(equation_name, Nx, Ny)

% DESCRIPTION
%   This function returns the index of the specified equation at the bottom-right corner of the domain in the eigenvalue problem matrices.
% INPUT:
%   equation_name       name of the equation                                                                                            [char]
% OUTPUT:
%   i_eqn_BR            row index of specified equation at the bottom-right corner of the domain in the eigenvalue problem matrices     [scalar]

i_eqn_TR = get_eqn_top_right_ind(equation_name, Nx, Ny);
i_eqn_BR = i_eqn_TR + (Ny-1);

end


function j_var_BR = get_var_bottom_right_ind(variable_name, Nx, Ny)

% DESCRIPTION
%   This function returns the index of the specified variable at the bottom-right corner of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the variable                                                                                            [char]
% OUTPUT:
%   j_var_BR            column index of specified variable at the bottom-right corner of the domain in the eigenvalue problem matrices  [scalar]

j_var_TR = get_var_top_right_ind(variable_name, Nx, Ny);
j_var_BR = j_var_TR + (Ny-1);

end


function i_eqn_TL = get_eqn_top_left_ind(equation_name, Nx, Ny)

% DESCRIPTION
%   This function returns the index of the specified equation at the top-left corner of the domain in the eigenvalue problem matrices.
% INPUT:
%   equation_name       name of the equation                                                                                        [char]
% OUTPUT:
%   i_eqn_TL            row index of specified equation at the top-left corner of the domain in the eigenvalue problem matrices     [scalar]

i_eqn_BL = get_eqn_bottom_left_ind(equation_name, Nx, Ny);
i_eqn_TL = i_eqn_BL - (Ny-1);

end


function j_var_TL = get_var_top_left_ind(variable_name, Nx, Ny)

% DESCRIPTION
%   This function returns the index of the specified variable at the top-left corner of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the variable                                                                                            [char]
% OUTPUT:
%   j_var_TL            column index of specified variable at the top-left corner of the domain in the eigenvalue problem matrices      [scalar]

j_var_BL = get_var_bottom_left_ind(variable_name, Nx, Ny);
j_var_TL = j_var_BL - (Ny-1);

end


function i_eqn_T = get_eqn_top_inds(equation_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified equation at the top boundary of the domain in the eigenvalue problem matrices.
% INPUT:
%   equation_name       name of the equation                                                                                    [char]
% OUTPUT:
%   i_eqn_T             row index of specified equation at the top boundary of the domain in the eigenvalue problem matrices    [scalar]

% nx = 1:Nx;
% i_eqn_TR = get_eqn_top_right_ind(equation_name, Nx, Ny);
% i_eqn_T = i_eqn_TR + (nx-1)*Ny;

i_eqn_TR = get_eqn_top_right_ind(equation_name, Nx, Ny);
i_eqn_TL = get_eqn_top_left_ind(equation_name, Nx, Ny);
i_eqn_T  = i_eqn_TR : Ny : i_eqn_TL;

end


function j_var_T = get_var_top_inds(variable_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified variable at the top boundary of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the variable                                                                                        [char]
% OUTPUT:
%   j_var_T             column index of specified variable at the top boundary of the domain in the eigenvalue problem matrices     [scalar]

% nx = 1:Nx;
% j_var_TR = get_var_top_right_ind(variable_name, Nx, Ny);
% j_var_T = i_var_TR + (nx-1)*Ny;

j_var_TR = get_var_top_right_ind(variable_name, Nx, Ny);
j_var_TL = get_var_top_left_ind(variable_name, Nx, Ny);
j_var_T  = j_var_TR : Ny : j_var_TL;

end


function i_eqn_B = get_eqn_bottom_inds(equation_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified equation at the bottom boundary of the domain in the eigenvalue problem matrices.
% INPUT:
%   equation_name       name of the equation                                                                                        [char]
% OUTPUT:
%   i_eqn_B             row index of specified equation at the bottom boundary of the domain in the eigenvalue problem matrices     [scalar]

% nx = 1:Nx;
% i_eqn_BR = get_eqn_bottom_right_ind(equation_name, Nx, Ny);
% i_eqn_B = i_eqn_BR + (nx-1)*Ny;

i_eqn_BR = get_eqn_bottom_right_ind(equation_name, Nx, Ny);
i_eqn_BL = get_eqn_bottom_left_ind(equation_name, Nx, Ny);
i_eqn_B  = i_eqn_BR : Ny : i_eqn_BL;

end


function j_var_B = get_var_bottom_inds(variable_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified variable at the bottom boundary of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the variable                                                                                        [char]
% OUTPUT:
%   j_var_B             column index of specified variable at the bottom boundary of the domain in the eigenvalue problem matrices  [scalar]

% nx = 1:Nx;
% j_var_BR = get_var_bottom_right_ind(variable_name, Nx, Ny);
% j_var_B = j_var_BR + (nx-1)*Ny;

j_var_BR = get_var_bottom_right_ind(variable_name, Nx, Ny);
j_var_BL = get_var_bottom_left_ind(variable_name, Nx, Ny);
j_var_B  = j_var_BR : Ny : j_var_BL;

end


function i_eqn_R = get_eqn_right_inds(equation_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified variable at the right boundary of the domain in the eigenvalue problem matrices.
% INPUT:
%   equation_name       name of the equation                                                                                        [char]
% OUTPUT:
%   i_eqn_R             row index of specified equation at the right boundary of the domain in the eigenvalue problem matrices      [scalar]

% ny = 1:Ny;
% i_eqn_TR = get_eqn_top_right_ind(equation_name, Nx, Ny);
% i_eqn_R = i_eqn_TR + (ny-1);

i_eqn_TR = get_eqn_top_right_ind(equation_name, Nx, Ny);
i_eqn_BR = get_eqn_bottom_right_ind(equation_name, Nx, Ny);
i_eqn_R  = i_eqn_TR : i_eqn_BR;

end


function j_var_R = get_var_right_inds(variable_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified variable at the right boundary of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the variable                                                                                        [char]
% OUTPUT:
%   j_var_R             column index of specified variable at the right boundary of the domain in the eigenvalue problem matrices     [scalar]

% ny = 1:Ny;
% j_var_TR = get_var_top_right_ind(variable_name, Nx, Ny);
% j_var_R = j_var_TR + (ny-1);

j_var_TR = get_var_top_right_ind(variable_name, Nx, Ny);
j_var_BR = get_var_bottom_right_ind(variable_name, Nx, Ny);
j_var_R = j_var_TR : j_var_BR;

end


function i_eqn_L = get_eqn_left_inds(equation_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified variable at the left boundary of the domain in the eigenvalue problem matrices.
% INPUT:
%   equation_name       name of the equation                                                                                    [char]
% OUTPUT:
%   i_eqn_L             row index of specified equation at the left boundary of the domain in the eigenvalue problem matrices   [scalar]

% ny = 1:Ny;
% i_eqn_TL = get_eqn_top_left_ind(equation_name, Nx, Ny);
% i_eqn_L = i_eqn_TL + (ny-1);

i_eqn_TL = get_eqn_top_left_ind(equation_name, Nx, Ny);
i_eqn_BL = get_eqn_bottom_left_ind(equation_name, Nx, Ny);
i_eqn_L  = i_eqn_TL : i_eqn_BL;

end


function j_var_L = get_var_left_inds(variable_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified variable at the left boundary of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the variable                                                                                        [char]
% OUTPUT:
%   j_var_L             column index of specified variable at the left boundary of the domain in the eigenvalue problem matrices     [scalar]

% ny = 1:Ny;
% j_var_TL = get_var_top_left_ind(variable_name, Nx, Ny);
% j_var_L = j_var_TL + (ny-1);

j_var_TL = get_var_top_left_ind(variable_name, Nx, Ny);
j_var_BL = get_var_bottom_left_ind(variable_name, Nx, Ny);
j_var_L = j_var_TL : j_var_BL;

end


function i_opr_TR = get_opr_top_right_ind(~, ~)

% DESCRIPTION
%   This function returns the index of the operator matrix at the top-right corner of the domain.
% INPUT:
%   (~, ~)              two insignificant inputs (ment for compatibility of input method with other similar functions)
% OUTPUT:
%   i_opr_TR            row index of the operator matrix at the top-right corner of the domain      [scalar]

i_opr_TR = 1;

end


function i_opr_BL = get_opr_bottom_left_ind(Nx, Ny)

% DESCRIPTION
%   This function returns the index of the operator matrix at the bottom-left corner of the domain.
% INPUT:
%   Nx                  Number of points along the chordwise direction                              [scalar]
%   Ny                  Number of points along the vertical direction                               [scalar]
% OUTPUT:
%   i_opr_BL            row index of the operator matrix at the bottom-left corner of the domain    [scalar]

i_opr_BL = Nx*Ny;

end


function i_opr_BR = get_opr_bottom_right_ind(Nx, Ny)

% DESCRIPTION
%   This function returns the index of the operator matrix at the bottom-right corner of the domain.
% INPUT:
%   Nx                  Number of points along the chordwise direction                              [scalar]
%   Ny                  Number of points along the vertical direction                               [scalar]
% OUTPUT:
%   i_opr_BR            row index of the operator matrix at the bottom-right corner of the domain   [scalar]

i_opr_TR = get_opr_top_right_ind(Nx, Ny);
i_opr_BR = i_opr_TR + (Ny-1);

end


function i_opr_TL = get_opr_top_left_ind(Nx, Ny)

% DESCRIPTION
%   This function returns the index of the operator matrix at the top-left corner of the domain.
% INPUT:
%   Nx                  Number of points along the chordwise direction                          [scalar]
%   Ny                  Number of points along the vertical direction                           [scalar]
% OUTPUT:
%   i_opr_TL            row index of the operator matrix at the top-left corner of the domain   [scalar]

i_opr_BL = get_opr_bottom_left_ind(Nx, Ny);
i_opr_TL = i_opr_BL - (Ny-1);

end


function i_opr_T = get_opr_top_inds(Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the operator matrix at the top boundary of the domain.
% INPUT:
%   Nx                  Number of points along the chordwise direction                          [scalar]
%   Ny                  Number of points along the vertical direction                           [scalar]
% OUTPUT:
%   i_opr_T             row indices of the operator matrix at the top boundary of the domain    [scalar]

% nx = 1:Nx;
% i_opr_TR = get_opr_top_right_ind(Nx, Ny);
% i_opr_T = i_opr_TR + (nx-1)*Ny;

i_opr_TR = get_opr_top_right_ind(Nx, Ny);
i_opr_BR = get_opr_bottom_right_ind(Nx, Ny);
i_opr_T  = i_opr_TR : Ny : i_opr_BR;

end


function i_opr_B = get_opr_bottom_inds(Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the operator matrix at the bottom boundary of the domain.
% INPUT:
%   Nx                  Number of points along the chordwise direction                              [scalar]
%   Ny                  Number of points along the vertical direction                               [scalar]
% OUTPUT:
%   i_opr_B             row indices of the operator matrix at the bottom boundary of the domain     [scalar]

% nx = 1:Nx;
% i_opr_BR = get_opr_bottom_right_ind(Nx, Ny);
% i_opr_B = i_opr_TR + (nx-1)*Ny;

i_opr_BR = get_opr_bottom_right_ind(Nx, Ny);
i_opr_BL = get_opr_bottom_left_ind(Nx, Ny);
i_opr_B  = i_opr_BR : Ny : i_opr_BL;

end


function i_opr_R = get_opr_right_inds(Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the operator matrix at the right boundary of the domain.
% INPUT:
%   Nx                  Number of points along the chordwise direction                          [scalar]
%   Ny                  Number of points along the vertical direction                           [scalar]
% OUTPUT:
%   i_opr_R             row indices of the operator matrix at the right boundary of the domain  [scalar]

% ny = 1:Ny;
% i_opr_TR = get_opr_top_right_ind(Nx, Ny);
% i_opr_R  = i_opr_TR + (ny-1);

i_opr_TR = get_opr_top_right_ind(Nx, Ny);
i_opr_BR = get_opr_bottom_right_ind(Nx, Ny);
i_opr_R  = i_opr_TR : i_opr_BR;

end


function i_opr_L = get_opr_left_inds(Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the operator matrix at the left boundary of the domain.
% INPUT:
%   Nx                  Number of points along the chordwise direction                          [scalar]
%   Ny                  Number of points along the vertical direction                           [scalar]
% OUTPUT:
%   i_opr_L             row indices of the operator matrix at the left boundary of the domain   [scalar]

% ny = 1:Ny;
% i_opr_TL = get_opr_top_left_ind(Nx, Ny);
% i_opr_L  = i_opr_TL + (ny-1);

i_opr_TL = get_opr_top_left_ind(Nx, Ny);
i_opr_BL = get_opr_bottom_left_ind(Nx, Ny);
i_opr_L  = i_opr_TL : i_opr_BL;

end



