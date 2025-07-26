function [A, B] = evpmatrices(Domain, Base_Flow, Problem)

% Factors for "placing" eigenvalues
dirichlet_factor              = -50;
neumann_factor                = -100;
linear_extrap_factor          = -150;
pressure_compatibility_factor = -200;
lppe_factor                   = -250;
symmetry_factor               = -300;
anti_symmetry_factor          = -350;

% Build eigenvalue problem A and B matrices
[A, B] = bglnsematrices(Domain, Base_Flow, Problem);

% Get spanwise wave number
beta = Problem.Physics.Beta;

% Build zero and identity matrices
Ny  = size(Domain.mat_X, 1);
Nx  = size(Domain.mat_X, 2);
Z   = zeros(Nx * Ny);
I   = eye(Nx * Ny);
iI  = 1i * I;

% Get differentiation matrices
xmat  = Domain.mat_X;
Dx    = Domain.Dx;
Dy    = Domain.Dy;
D2x   = Domain.D2x;
D2y   = Domain.D2y;

% Various matrix elements
dphi  = repmat(flip(Base_Flow.dphi) , [1 Nx]);
ddphi = repmat(flip(Base_Flow.ddphi), [1 Nx]);
Ux    = diag(dphi(:), 0);
Uy    = xmat .* ddphi;
Uy    = diag(Uy(:), 0);
Vy    = diag(-dphi(:), 0);


%% Wall boundary conditions

% Multiple entry options allow for flexible user interaction with the code
Wall_Options.No_Slip        = {'Dirichlet' 'No_Slip' 'No-Slip' 'No Slip'}; % Meant for 'u' and 'w'
Wall_Options.No_Penetration = {'Dirichlet' 'No_Penetration' 'No-Penetration' 'No Penetration'}; % Meant for 'v'
Wall_Options.PC             = {'PC' 'Pressure-Compatibility' 'Pressure Compatibility' 'Pressure_Compatibility'}; % Meant for 'p'
Wall_Options.LPPE           = {'LPPE' 'Linearized-Pressure-Poisson-Equation' 'Linearized Pressure Poisson Equation' 'Linearized_Pressure_Poisson_Equation'}; % Meant for 'p'
Wall_Options.Dirichlet      = {'Dirichlet'}; % JUST FOR TESTS
Wall_Options.Neumann        = {'Neumann'}; % JUST FOR TESTS

% u
switch Problem.Boundary_Conditions.Wall.u
    case Wall_Options.No_Slip % No-slip boundary condition
        row_inds    = get_equation_bottom_indices('x momentum', Nx, Ny);
        column_inds = get_variable_bottom_indices('u', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Wall.u ' for ''u'' at the wall is invalid or not supported'])
end

% v
switch Problem.Boundary_Conditions.Wall.v
    case Wall_Options.No_Penetration % No-penetration boundary condition
        row_inds    = get_equation_bottom_indices('y momentum', Nx, Ny);
        column_inds = get_variable_bottom_indices('v', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Wall.v ' for ''v'' at the wall is invalid or not supported'])
end

% w
switch Problem.Boundary_Conditions.Wall.w
    case Wall_Options.No_Slip % No-slip boundary condition
        row_inds    = get_equation_bottom_indices('z momentum', Nx, Ny);
        column_inds = get_variable_bottom_indices('w', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Wall.w ' for ''w'' at the wall is invalid or not supported'])
end

% p
switch Problem.Boundary_Conditions.Wall.p
    case Wall_Options.PC % Pressure compatibility boundary condition
        row_inds = get_equation_bottom_indices('continuity', Nx, Ny);
        i_opr_B  = get_operator_bottom_indices(Nx, Ny);

        DxDy = Dx * Dy;
        ibIDy = (1i * beta * I) * Dy;
        pressure_compatibility_opr = [-DxDy(i_opr_B,:) , Z(i_opr_B,:) , -ibIDy(i_opr_B,:) , -Dy(i_opr_B,:)];

        % pressure_compatibility_opr = [Z(i_opr_B,:) , D2y(i_opr_B,:) , Z(i_opr_B,:) , -Dy(i_opr_B,:)];

        A(row_inds(:),:) = 0;
        B(row_inds(:),:) = 0;
        A(row_inds(:),:) = pressure_compatibility_factor * pressure_compatibility_opr;
        B(row_inds(:),:) = -1i * pressure_compatibility_opr;
    case Wall_Options.LPPE % LPPE boundary condition
        row_inds = get_equation_bottom_indices('continuity', Nx, Ny);
        i_opr_B  = get_operator_bottom_indices(Nx, Ny);

        lppe_u = Z;
        lppe_v = Z;
        lppe_w = Z;
        lppe_p = D2x + D2y - (beta^2 * I);
        
        lppe_opr = [lppe_u(i_opr_B,:) , lppe_v(i_opr_B,:) , lppe_w(i_opr_B,:) , lppe_p(i_opr_B,:)];

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(row_inds(:),:) =  lppe_factor * lppe_opr;
        B(row_inds(:),:) = -1i * lppe_opr;
    case Wall_Options.Dirichlet % Dirichlet boundary condition - USE FOR TESTING ONLY
        row_inds    = get_equation_bottom_indices('continuity', Nx, Ny);
        column_inds = get_variable_bottom_indices('p', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    case Wall_Options.Neumann % Neumann boundary condition - USE FOR TESTING ONLY
        operator_row_inds = get_operator_bottom_indices(Nx, Ny);
        row_inds          = get_equation_bottom_indices('continuity', Nx, Ny);

        neumann_opr = [Z(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:) Dy(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  neumann_factor * neumann_opr;
        B(row_inds,:) = -1i * neumann_opr;
    otherwise
        warning('No boundary condition applied for the pressure on the wall')
end


%% Top boundary conditions

Top_Options.Decay = {'Dirichlet' 'Decay' 'decay'};
Top_Options.LPPE  = {'LPPE' 'Linearized-Pressure-Poisson-Equation' 'Linearized Pressure Poisson Equation'};
Top_Options.PC    = {'PC' 'Pressure-Compatibility' 'Pressure Compatibility' 'Pressure_Compatibility'};

% u
switch Problem.Boundary_Conditions.Top.u
    case Top_Options.Decay % u(y-->infinity) = 0
        row_inds    = get_equation_top_indices('x momentum', Nx, Ny);
        column_inds = get_variable_top_indices('u', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Wall.u ' for ''u'' at the top is invalid or not supported'])
end

% v
switch Problem.Boundary_Conditions.Top.v
    case Top_Options.Decay % v(y-->infinity) = 0
        row_inds    = get_equation_top_indices('y momentum', Nx, Ny);
        column_inds = get_variable_top_indices('v', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Wall.v ' for ''v'' at the top is invalid or not supported'])
end

% w
switch Problem.Boundary_Conditions.Top.w
    case Top_Options.Decay % w(y-->infinity) = 0
        row_inds    = get_equation_top_indices('z momentum', Nx, Ny);
        column_inds = get_variable_top_indices('w', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Wall.w ' for ''w'' at the top is invalid or not supported'])
end

% p
switch Problem.Boundary_Conditions.Top.p
    case Top_Options.Decay % p(y-->infinity) = 0
        row_inds    = get_equation_top_indices('continuity', Nx, Ny);
        column_inds = get_variable_top_indices('p', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    case Top_Options.PC
        row_inds = get_equation_top_indices('continuity', Nx, Ny);
        i_opr_T  = get_operator_top_indices(Nx, Ny);

        DxDy = Dx * Dy;
        ibIDy = (1i * beta * I) * Dy;
        pressure_compatibility_opr = [-DxDy(i_opr_T,:) , Z(i_opr_T,:) , -ibIDy(i_opr_T,:) , -Dy(i_opr_T,:)];

        % pressure_compatibility_opr = [Z(i_opr_B,:) , D2y(i_opr_B,:) , Z(i_opr_B,:) , -Dy(i_opr_B,:)];

        A(row_inds(:),:) = 0;
        B(row_inds(:),:) = 0;
        A(row_inds(:),:) = pressure_compatibility_factor * pressure_compatibility_opr;
        B(row_inds(:),:) = -1i * pressure_compatibility_opr;
    case Top_Options.LPPE % LPPE boundary condition
        row_inds = get_equation_top_indices('continuity', Nx, Ny);
        i_opr_T  = get_operator_top_indices(Nx, Ny);

        lppe_u = Z;
        lppe_v = Z;
        lppe_w = Z;
        lppe_p = D2x + D2y - (beta^2 * I);

        lppe_opr = [lppe_u(i_opr_T,:) , lppe_v(i_opr_T,:) , lppe_w(i_opr_T,:) , lppe_p(i_opr_T,:)];

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(row_inds(:),:) =  lppe_factor * lppe_opr;
        B(row_inds(:),:) = -1i * lppe_opr;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Wall.p ' for ''p'' at the top is invalid or not supported'])
end


%% Right side

Right_Side_Options.FD_Extrapolation = {'Linear_Extrapolation' 'Linear-Extrapolation' 'Linear Extrapolation' 'FD_Extrapolation' 'FD Extrapolation' 'Finite Difference Extrapolation' 'finite_difference_extrapolation' 'finite difference extrapolation'};
Right_Side_Options.Zero_2nd_Derivative_Extrapolation = {'zero_2nd_derivative' 'zero_2nd_derivative_extrapolation'};
Right_Side_Options.LPPE = {'LPPE' 'Linearized-Pressure-Poisson-Equation' 'Linearized Pressure Poisson Equation'};
Right_Side_Options.Symmetry = {'Symmetry' 'symmetry' 'sym' 's' 'S' 'Sym'};
Right_Side_Options.Anti_Symmetry = {'Anti_Symmetry' 'Anti_Symmetry' 'AntiSymmetry' 'Antisymmetry' 'asym' 'as' 'AS' 'a' 'A' 'ASym'};
Right_Side_Options.Symmetry = {'Symmetry' 'symmetry' 'sym' 's' 'S' 'Sym'};
Right_Side_Options.Dirichlet = {'Dirichlet'};
Right_Side_Options.Neumann  = {'Neumann'};

% u
switch Problem.Boundary_Conditions.Right.u
    case Right_Side_Options.FD_Extrapolation % Linear extrapolation
        j_u_right_inds    = get_variable_right_indices('u', Nx, Ny);
        i_xmom_right_inds = get_equation_right_indices('x momentum', Nx, Ny);

        A(i_xmom_right_inds(2:end-1),:) = 0;
        B(i_xmom_right_inds(2:end-1),:) = 0;
        
        ind_shift = Ny*[0 1 2];
        for i = 2:Ny-1
            C = (xmat(i,1)-xmat(i,2))/(xmat(i,2)-xmat(i,3));
            linear_extrap_opr = [1 -C-1 C];
            
            A(i_xmom_right_inds(i),j_u_right_inds(i) + ind_shift) = linear_extrap_factor*linear_extrap_opr;
            B(i_xmom_right_inds(i),j_u_right_inds(i) + ind_shift) = -1i*linear_extrap_opr;
        end
    case Right_Side_Options.Zero_2nd_Derivative_Extrapolation % Zero 2nd derivative extrapolation
        operator_row_inds = get_operator_right_indices(Nx, Ny);
        row_inds          = get_equation_right_indices('x momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        z_2nd_der_opr = [D2x(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  linear_extrap_factor * z_2nd_der_opr;
        B(row_inds,:) = -1i * z_2nd_der_opr;
    case Right_Side_Options.Symmetry
        % TBD
    case Right_Side_Options.Anti_Symmetry
        % TBD
    case Right_Side_Options.Neumann
        operator_row_inds = get_operator_right_indices(Nx, Ny);
        row_inds          = get_equation_right_indices('x momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        neumann_opr = [Dx(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  neumann_factor * neumann_opr;
        B(row_inds,:) = -1i * neumann_opr;
    case Right_Side_Options.Dirichlet
        row_inds    = get_equation_right_indices('x momentum', Nx, Ny);
        column_inds = get_variable_right_indices('u', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Right.u ' for ''u'' at the right side is invalid or not supported'])
end

% v
switch Problem.Boundary_Conditions.Right.v
    case Right_Side_Options.FD_Extrapolation % Linear extrapolation
        j_v_right_inds    = get_variable_right_indices('v', Nx, Ny);
        i_ymom_right_inds = get_equation_right_indices('y momentum', Nx, Ny);
        
        A(i_ymom_right_inds(2:end-1),:) = 0;
        B(i_ymom_right_inds(2:end-1),:) = 0;
        
        ind_shift = Ny*[0 1 2];
        for i = 2:Ny-1
            C = (xmat(i,1)-xmat(i,2))/(xmat(i,2)-xmat(i,3));
            linear_extrap_opr = [1 -C-1 C];
            
            A(i_ymom_right_inds(i),j_v_right_inds(i) + ind_shift) = linear_extrap_factor*linear_extrap_opr;
            B(i_ymom_right_inds(i),j_v_right_inds(i) + ind_shift) = -1i*linear_extrap_opr;
        end
    case Right_Side_Options.Zero_2nd_Derivative_Extrapolation % Zero 2nd derivative extrapolation
        operator_row_inds = get_operator_right_indices(Nx, Ny);
        row_inds          = get_equation_right_indices('y momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        z_2nd_der_opr = [Z(operator_row_inds,:) D2x(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  linear_extrap_factor * z_2nd_der_opr;
        B(row_inds,:) = -1i * z_2nd_der_opr;
    case Right_Side_Options.Symmetry
        % TBD
    case Right_Side_Options.Anti_Symmetry
        % TBD
    case Right_Side_Options.Neumann
        operator_row_inds = get_operator_right_indices(Nx, Ny);
        row_inds          = get_equation_right_indices('y momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        neumann_opr = [Z(operator_row_inds,:) Dx(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  neumann_factor * neumann_opr;
        B(row_inds,:) = -1i * neumann_opr;
    case Right_Side_Options.Dirichlet
        row_inds    = get_equation_right_indices('y momentum', Nx, Ny);
        column_inds = get_variable_right_indices('v', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Right.v ' for ''v'' at the right side is invalid or not supported'])
end

% w
switch Problem.Boundary_Conditions.Right.w
    case Right_Side_Options.FD_Extrapolation % Linear extrapolation
        j_w_right_inds    = get_variable_right_indices('w', Nx, Ny);
        i_zmom_right_inds = get_equation_right_indices('z momentum', Nx, Ny);
        
        A(i_zmom_right_inds(2:end-1),:) = 0;
        B(i_zmom_right_inds(2:end-1),:) = 0;
        
        ind_shift = Ny*[0 1 2];
        for i = 2:Ny-1
            C = (xmat(i,1)-xmat(i,2))/(xmat(i,2)-xmat(i,3));
            linear_extrap_opr = [1 -C-1 C];
            
            A(i_zmom_right_inds(i),j_w_right_inds(i) + ind_shift) = linear_extrap_factor*linear_extrap_opr;
            B(i_zmom_right_inds(i),j_w_right_inds(i) + ind_shift) = -1i*linear_extrap_opr;
        end
    case Right_Side_Options.Zero_2nd_Derivative_Extrapolation % Zero 2nd derivative extrapolation
        operator_row_inds = get_operator_right_indices(Nx, Ny);
        row_inds          = get_equation_right_indices('z momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        z_2nd_der_opr = [Z(operator_row_inds,:) Z(operator_row_inds,:) D2x(operator_row_inds,:) Z(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  linear_extrap_factor * z_2nd_der_opr;
        B(row_inds,:) = -1i * z_2nd_der_opr;
    case Right_Side_Options.Symmetry
        % TBD
    case Right_Side_Options.Anti_Symmetry
        % TBD
    case Right_Side_Options.Neumann
        operator_row_inds = get_operator_right_indices(Nx, Ny);
        row_inds          = get_equation_right_indices('z momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        neumann_opr = [Z(operator_row_inds,:) Z(operator_row_inds,:) Dx(operator_row_inds,:) Z(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  neumann_factor * neumann_opr;
        B(row_inds,:) = -1i * neumann_opr;
    case Right_Side_Options.Dirichlet
        row_inds    = get_equation_right_indices('z momentum', Nx, Ny);
        column_inds = get_variable_right_indices('w', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Right.w ' for ''w'' at the right side is invalid or not supported'])
end

% p
switch Problem.Boundary_Conditions.Right.p
    case Right_Side_Options.FD_Extrapolation % Linear extrapolation
        j_p_right_inds    = get_variable_right_indices('p', Nx, Ny);
        i_cont_right_inds = get_equation_right_indices('continuity', Nx, Ny);
        
        A(i_cont_right_inds(2:end-1),:) = 0;
        B(i_cont_right_inds(2:end-1),:) = 0;

        ind_shift = Ny*[0 1 2];
        for i = 2:Ny-1
            C = (xmat(i,1)-xmat(i,2))/(xmat(i,2)-xmat(i,3));
            linear_extrap_opr = [1 -C-1 C];

            A(i_cont_right_inds(i),j_p_right_inds(i) + ind_shift) = linear_extrap_factor*linear_extrap_opr;
            B(i_cont_right_inds(i),j_p_right_inds(i) + ind_shift) = -1i*linear_extrap_opr;
        end
    case Right_Side_Options.Zero_2nd_Derivative_Extrapolation % Zero 2nd derivative extrapolation
        operator_row_inds = get_operator_right_indices(Nx, Ny);
        row_inds          = get_equation_right_indices('continuity', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        z_2nd_der_opr = [Z(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:) D2x(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  linear_extrap_factor * z_2nd_der_opr;
        B(row_inds,:) = -1i * z_2nd_der_opr;
    case Right_Side_Options.LPPE % LPPE boundary condition
        row_inds = get_equation_right_indices('continuity', Nx, Ny);
        i_opr_B  = get_operator_right_indices(Nx, Ny);
        row_inds = row_inds(2:end-1);
        i_opr_B  = i_opr_B(2:end-1);

        lppe_u = 2*Ux*Dx;
        lppe_v = 2*(Uy*Dx + Vy*Dy);
        lppe_w = Z;
        lppe_p = D2x + D2y - beta^2*I;

        lppe_opr = [lppe_u(i_opr_B,:) , lppe_v(i_opr_B,:) , lppe_w(i_opr_B,:) , lppe_p(i_opr_B,:)];

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(row_inds(:),:) =  lppe_factor * lppe_opr;
        B(row_inds(:),:) = -1i * lppe_opr;
    case Right_Side_Options.Symmetry
        % TBD
    case Right_Side_Options.Anti_Symmetry
        % TBD
    case Right_Side_Options.Neumann
        operator_row_inds = get_operator_right_indices(Nx, Ny);
        row_inds          = get_equation_right_indices('continuity', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        neumann_opr = [Z(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:) Dx(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  neumann_factor * neumann_opr;
        B(row_inds,:) = -1i * neumann_opr;
    case Right_Side_Options.Dirichlet
        row_inds    = get_equation_right_indices('continuity', Nx, Ny);
        column_inds = get_variable_right_indices('p', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Right.p ' for ''p'' at the right side is invalid or not supported'])
end


%% Left side

Left_Side_Options.FD_Extrapolation = {'Linear_Extrapolation' 'Linear-Extrapolation' 'Linear Extrapolation' 'FD_Extrapolation' 'FD Extrapolation' 'Finite Difference Extrapolation' 'finite_difference_extrapolation' 'finite difference extrapolation'};
Left_Side_Options.Zero_2nd_Derivative_Extrapolation = {'zero_2nd_derivative' 'zero_2nd_derivative_extrapolation'};
Left_Side_Options.LPPE = {'LPPE' 'Linearized-Pressure-Poisson-Equation' 'Linearized Pressure Poisson Equation'};
Left_Side_Options.Symmetry = {'Symmetry' 'symmetry' 'sym' 's' 'S' 'Sym'};
Left_Side_Options.Anti_Symmetry = {'Anti_Symmetry' 'Anti_Symmetry' 'AntiSymmetry' 'Antisymmetry' 'asym' 'as' 'AS' 'a' 'A' 'ASym'};
Left_Side_Options.Dirichlet = {'Dirichlet'};
Left_Side_Options.Neumann = {'Neumann'};

% u
switch Problem.Boundary_Conditions.Left.u
    case Left_Side_Options.FD_Extrapolation % Linear extrapolation
        j_u_left_inds    = get_variable_left_indices('u', Nx, Ny);
        i_xmom_left_inds = get_equation_left_indices('x momentum', Nx, Ny);

        A(i_xmom_left_inds(2:end-1),:) = 0;
        B(i_xmom_left_inds(2:end-1),:) = 0;

        ind_shift = -Ny*[2 1 0];
        for i = 2:Ny-1
            C = (xmat(i,Nx)-xmat(i,Nx-1))/(xmat(i,Nx-1)-xmat(i,Nx-2));
            linear_extrap_opr = [C -C-1 1];
            
            A(i_xmom_left_inds(i),j_u_left_inds(i) + ind_shift) = linear_extrap_factor*linear_extrap_opr;
            B(i_xmom_left_inds(i),j_u_left_inds(i) + ind_shift) = -1i*linear_extrap_opr;
        end
    case Left_Side_Options.Zero_2nd_Derivative_Extrapolation % Zero 2nd derivative extrapolation
        operator_row_inds = get_operator_left_indices(Nx, Ny);
        row_inds          = get_equation_left_indices('x momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        z_2nd_der_opr = [D2x(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  linear_extrap_factor * z_2nd_der_opr;
        B(row_inds,:) = -1i * z_2nd_der_opr;
    case Left_Side_Options.Symmetry
        opr_left_inds = get_operator_whole_left_indices(Nx, Ny);
        eqn_left_inds = get_equation_whole_left_indices('x momentum', Nx, Ny);
        jm = (Nx - 1) / 2 + 1;
        opr = Z;
        for k = 1 : jm - 1
            left_ind = (jm + k - 1) * Ny + (1 : Ny);
            right_ind = (jm - k - 1) * Ny + (1 : Ny);
            for n = 1 : length(left_ind)
                opr(left_ind(n), left_ind(n)) = 1;
                opr(left_ind(n), right_ind(n)) = -1;
            end
        end

        symmetry_opr = [opr(opr_left_inds,:) Z(opr_left_inds,:) Z(opr_left_inds,:) Z(opr_left_inds,:)];
        
        A(eqn_left_inds,:) =  0;
        B(eqn_left_inds,:) =  0;
        A(eqn_left_inds,:) =  symmetry_factor * symmetry_opr;
        B(eqn_left_inds,:) = -1i * symmetry_opr;

        % opr_mid_inds = (jm - 1) * Ny + (1 : Ny);
        % opr_mid_inds = opr_mid_inds(2:end-1);
        % 
        % u_symmetry_opr = [Dx(opr_mid_inds,:) Z(opr_mid_inds,:) Z(opr_mid_inds,:) Z(opr_mid_inds,:)];
        % 
        % eqn_mid_inds = opr_mid_inds + 0 * Nx * Ny;
        % 
        % A(eqn_mid_inds,:) =  0;
        % B(eqn_mid_inds,:) =  0;
        % A(eqn_mid_inds,:) =  symmetry_factor * u_symmetry_opr;
        % B(eqn_mid_inds,:) = -1i * u_symmetry_opr;
    case Left_Side_Options.Anti_Symmetry
        opr_left_inds = get_operator_whole_left_indices(Nx, Ny);
        eqn_left_inds = get_equation_whole_left_indices('x momentum', Nx, Ny);
        jm = (Nx - 1) / 2 + 1;
        opr = Z;
        for k = 1 : jm - 1
            left_ind = (jm + k - 1) * Ny + (1 : Ny);
            right_ind = (jm - k - 1) * Ny + (1 : Ny);
            for n = 1 : length(left_ind)
                opr(left_ind(n), left_ind(n)) = 1;
                opr(left_ind(n), right_ind(n)) = 1;
            end
        end

        symmetry_opr = [opr(opr_left_inds,:) Z(opr_left_inds,:) Z(opr_left_inds,:) Z(opr_left_inds,:)];
        
        A(eqn_left_inds,:) =  0;
        B(eqn_left_inds,:) =  0;
        A(eqn_left_inds,:) =  symmetry_factor * symmetry_opr;
        B(eqn_left_inds,:) = -1i * symmetry_opr;

        % opr_mid_inds = (jm - 1) * Ny + (1 : Ny);
        % opr_mid_inds = opr_mid_inds(2:end-1);
        % 
        % row_inds    = opr_mid_inds;
        % column_inds = opr_mid_inds;
        % linear_inds = get_linear_indices(A, row_inds, column_inds);
        % 
        % A(row_inds(:),:) =  0;
        % B(row_inds(:),:) =  0;
        % A(linear_inds)   =  dirichlet_factor;
        % B(linear_inds)   = -1i;
    case Left_Side_Options.Neumann
        operator_row_inds = get_operator_left_indices(Nx, Ny);
        row_inds          = get_equation_left_indices('x momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        neumann_opr = [Dx(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  neumann_factor * neumann_opr;
        B(row_inds,:) = -1i * neumann_opr;
    case Left_Side_Options.Dirichlet
        row_inds    = get_equation_left_indices('x momentum', Nx, Ny);
        column_inds = get_variable_left_indices('u', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Left.u ' for ''u'' at the left side is invalid or not supported'])
end

% v
switch Problem.Boundary_Conditions.Left.v
    case Left_Side_Options.FD_Extrapolation % Linear extrapolation
        j_v_left_inds    = get_variable_left_indices('v', Nx, Ny);
        i_ymom_left_inds = get_equation_left_indices('y momentum', Nx, Ny);

        A(i_ymom_left_inds(2:end-1),:) = 0;
        B(i_ymom_left_inds(2:end-1),:) = 0;

        ind_shift = -Ny*[2 1 0];
        for i = 2:Ny-1
            C = (xmat(i,Nx)-xmat(i,Nx-1))/(xmat(i,Nx-1)-xmat(i,Nx-2));
            linear_extrap_opr = [C -C-1 1];
            
            A(i_ymom_left_inds(i),j_v_left_inds(i) + ind_shift) = linear_extrap_factor*linear_extrap_opr;
            B(i_ymom_left_inds(i),j_v_left_inds(i) + ind_shift) = -1i*linear_extrap_opr;
        end
    case Left_Side_Options.Zero_2nd_Derivative_Extrapolation % Zero 2nd derivative extrapolation
        operator_row_inds = get_operator_left_indices(Nx, Ny);
        row_inds          = get_equation_left_indices('y momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        z_2nd_der_opr = [Z(operator_row_inds,:) D2x(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  linear_extrap_factor * z_2nd_der_opr;
        B(row_inds,:) = -1i * z_2nd_der_opr;
    case Left_Side_Options.Symmetry
        opr_left_inds = get_operator_whole_left_indices(Nx, Ny);
        eqn_left_inds = get_equation_whole_left_indices('y momentum', Nx, Ny);
        jm = (Nx - 1) / 2 + 1;
        opr = Z;
        for k = 1 : jm - 1
            left_ind = (jm + k - 1) * Ny + (1 : Ny);
            right_ind = (jm - k - 1) * Ny + (1 : Ny);
            for n = 1 : length(left_ind)
                opr(left_ind(n), left_ind(n)) = 1;
                opr(left_ind(n), right_ind(n)) = 1;
            end
        end

        symmetry_opr = [Z(opr_left_inds,:) opr(opr_left_inds,:) Z(opr_left_inds,:) Z(opr_left_inds,:)];
        
        A(eqn_left_inds,:) =  0;
        B(eqn_left_inds,:) =  0;
        A(eqn_left_inds,:) =  symmetry_factor * symmetry_opr;
        B(eqn_left_inds,:) = -1i * symmetry_opr;

        opr_mid_inds = (jm - 1) * Ny + (1 : Ny) + (1 * Nx * Ny);
        opr_mid_inds = opr_mid_inds(2:end-1);
        
        row_inds    = opr_mid_inds;
        column_inds = opr_mid_inds;
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    case Left_Side_Options.Anti_Symmetry
        opr_left_inds = get_operator_whole_left_indices(Nx, Ny);
        eqn_left_inds = get_equation_whole_left_indices('y momentum', Nx, Ny);
        jm = (Nx - 1) / 2 + 1;
        opr = Z;
        for k = 1 : jm - 1
            left_ind = (jm + k - 1) * Ny + (1 : Ny);
            right_ind = (jm - k - 1) * Ny + (1 : Ny);
            for n = 1 : length(left_ind)
                opr(left_ind(n), left_ind(n)) = 1;
                opr(left_ind(n), right_ind(n)) = -1;
            end
        end

        symmetry_opr = [Z(opr_left_inds,:) opr(opr_left_inds,:) Z(opr_left_inds,:) Z(opr_left_inds,:)];
        
        A(eqn_left_inds,:) = 0;
        B(eqn_left_inds,:) = 0;
        A(eqn_left_inds,:) = symmetry_factor * symmetry_opr;
        B(eqn_left_inds,:) = -1i * symmetry_opr;

        % opr_mid_inds = (jm - 1) * Ny + (1 : Ny);
        % opr_mid_inds = opr_mid_inds(2:end-1);
        % 
        % v_symmetry_opr = [Z(opr_mid_inds,:) Dx(opr_mid_inds,:) Z(opr_mid_inds,:) Z(opr_mid_inds,:)];
        % 
        % eqn_mid_inds = opr_mid_inds + (1 * Nx * Ny);
        % 
        % A(eqn_mid_inds,:) =  0;
        % B(eqn_mid_inds,:) =  0;
        % A(eqn_mid_inds,:) =  symmetry_factor * v_symmetry_opr;
        % B(eqn_mid_inds,:) = -1i * v_symmetry_opr;
    case Left_Side_Options.Neumann
        operator_row_inds = get_operator_left_indices(Nx, Ny);
        row_inds          = get_equation_left_indices('y momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        neumann_opr = [Z(operator_row_inds,:) Dx(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  neumann_factor * neumann_opr;
        B(row_inds,:) = -1i * neumann_opr;
    case Left_Side_Options.Dirichlet
        row_inds    = get_equation_left_indices('y momentum', Nx, Ny);
        column_inds = get_variable_left_indices('v', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Left.v ' for ''v'' at the left side is invalid or not supported'])
end

% w
switch Problem.Boundary_Conditions.Left.w
    case Left_Side_Options.FD_Extrapolation % Linear extrapolation
        j_w_left_inds    = get_variable_left_indices('w', Nx, Ny);
        i_zmom_left_inds = get_equation_left_indices('z momentum', Nx, Ny);

        A(i_zmom_left_inds(2:end-1),:) = 0;
        B(i_zmom_left_inds(2:end-1),:) = 0;

        ind_shift = -Ny*[2 1 0];
        for i = 2:Ny-1
            C = (xmat(i,Nx)-xmat(i,Nx-1))/(xmat(i,Nx-1)-xmat(i,Nx-2));
            linear_extrap_opr = [C -C-1 1];
            
            A(i_zmom_left_inds(i),j_w_left_inds(i) + ind_shift) = linear_extrap_factor*linear_extrap_opr;
            B(i_zmom_left_inds(i),j_w_left_inds(i) + ind_shift) = -1i*linear_extrap_opr;
        end
    case Left_Side_Options.Zero_2nd_Derivative_Extrapolation % Zero 2nd derivative extrapolation
        operator_row_inds = get_operator_left_indices(Nx, Ny);
        row_inds          = get_equation_left_indices('z momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        z_2nd_der_opr = [Z(operator_row_inds,:) Z(operator_row_inds,:) D2x(operator_row_inds,:) Z(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  linear_extrap_factor * z_2nd_der_opr;
        B(row_inds,:) = -1i * z_2nd_der_opr;
    case Left_Side_Options.Symmetry
        opr_left_inds = get_operator_left_indices(Nx, Ny);
        eqn_left_inds = get_equation_left_indices('z momentum', Nx, Ny);
        jm = (Nx - 1) / 2 + 1;
        opr = Z;
        for k = 1 : jm - 1
            left_ind = (jm + k - 1) * Ny + (1 : Ny);
            right_ind = (jm - k - 1) * Ny + (1 : Ny);
            for n = 1 : length(left_ind)
                opr(left_ind(n), left_ind(n)) = 1;
                opr(left_ind(n), right_ind(n)) = 1;
            end
        end

        symmetry_opr = [Z(opr_left_inds,:) Z(opr_left_inds,:) opr(opr_left_inds,:) Z(opr_left_inds,:)];
        
        A(eqn_left_inds,:) =  0;
        B(eqn_left_inds,:) =  0;
        A(eqn_left_inds,:) =  symmetry_factor * symmetry_opr;
        B(eqn_left_inds,:) = -1i * symmetry_opr;

        % opr_mid_inds = (jm - 1) * Ny + (1 : Ny) + (2 * Nx * Ny);
        % opr_mid_inds = opr_mid_inds(2:end-1);
        % 
        % row_inds    = opr_mid_inds;
        % column_inds = opr_mid_inds;
        % linear_inds = get_linear_indices(A, row_inds, column_inds);
        % 
        % A(row_inds(:),:) =  0;
        % B(row_inds(:),:) =  0;
        % A(linear_inds)   =  dirichlet_factor;
        % B(linear_inds)   = -1i;
    case Left_Side_Options.Anti_Symmetry
        opr_left_inds = get_operator_whole_left_indices(Nx, Ny);
        eqn_left_inds = get_equation_whole_left_indices('z momentum', Nx, Ny);
        jm = (Nx - 1) / 2 + 1;
        opr = Z;
        for k = 1 : jm - 1
            left_ind = (jm + k - 1) * Ny + (1 : Ny);
            right_ind = (jm - k - 1) * Ny + (1 : Ny);
            for n = 1 : length(left_ind)
                opr(left_ind(n), left_ind(n)) = 1;
                opr(left_ind(n), right_ind(n)) = -1;
            end
        end

        symmetry_opr = [Z(opr_left_inds,:) Z(opr_left_inds,:) opr(opr_left_inds,:) Z(opr_left_inds,:)];
        
        A(eqn_left_inds,:) =  0;
        B(eqn_left_inds,:) =  0;
        A(eqn_left_inds,:) =  symmetry_factor * symmetry_opr;
        B(eqn_left_inds,:) = -1i * symmetry_opr;

        % opr_mid_inds = (jm - 1) * Ny + (1 : Ny);
        % opr_mid_inds = opr_mid_inds(2:end-1);
        % 
        % w_symmetry_opr = [Z(opr_mid_inds,:) Z(opr_mid_inds,:) Dx(opr_mid_inds,:) Z(opr_mid_inds,:)];
        % 
        % eqn_mid_inds = opr_mid_inds + (2 * Nx * Ny);
        % 
        % A(eqn_mid_inds,:) =  0;
        % B(eqn_mid_inds,:) =  0;
        % A(eqn_mid_inds,:) =  symmetry_factor * w_symmetry_opr;
        % B(eqn_mid_inds,:) = -1i * w_symmetry_opr;
    case Left_Side_Options.Neumann
        operator_row_inds = get_operator_left_indices(Nx, Ny);
        row_inds          = get_equation_left_indices('z momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        neumann_opr = [Z(operator_row_inds,:) Z(operator_row_inds,:) Dx(operator_row_inds,:) Z(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  neumann_factor * neumann_opr;
        B(row_inds,:) = -1i * neumann_opr;
    case Left_Side_Options.Dirichlet
        row_inds    = get_equation_left_indices('z momentum', Nx, Ny);
        column_inds = get_variable_left_indices('w', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Left.w ' for ''w'' at the left side is invalid or not supported'])
end

% p
switch Problem.Boundary_Conditions.Left.p
    case Left_Side_Options.FD_Extrapolation % Linear extrapolation
        j_p_left_inds    = get_variable_left_indices('p', Nx, Ny);
        i_cont_left_inds = get_equation_left_indices('continuity', Nx, Ny);

        A(i_cont_left_inds(2:end-1),:) = 0;
        B(i_cont_left_inds(2:end-1),:) = 0;

        ind_shift = -Ny*[2 1 0];
        for i = 2:Ny-1
            C = (xmat(i,Nx)-xmat(i,Nx-1))/(xmat(i,Nx-1)-xmat(i,Nx-2));
            linear_extrap_opr = [C -C-1 1];
            
            A(i_cont_left_inds(i),j_p_left_inds(i) + ind_shift) = linear_extrap_factor*linear_extrap_opr;
            B(i_cont_left_inds(i),j_p_left_inds(i) + ind_shift) = -1i*linear_extrap_opr;
        end
    case Left_Side_Options.Zero_2nd_Derivative_Extrapolation % Zero 2nd derivative extrapolation
        operator_row_inds = get_operator_left_indices(Nx, Ny);
        row_inds          = get_equation_left_indices('continuity', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        z_2nd_der_opr = [Z(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:) D2x(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  linear_extrap_factor * z_2nd_der_opr;
        B(row_inds,:) = -1i * z_2nd_der_opr;
    case Left_Side_Options.LPPE % LPPE boundary condition
        i_opr_B  = get_operator_left_indices(Nx, Ny);
        row_inds = get_equation_left_indices('continuity', Nx, Ny);
        row_inds = row_inds(2:end-1);
        i_opr_B  = i_opr_B(2:end-1);

        lppe_u = 2*Ux*Dx;
        lppe_v = 2*(Uy*Dx + Vy*Dy);
        lppe_w = Z;
        lppe_p = D2x + D2y - beta^2*I;

        lppe_opr = [lppe_u(i_opr_B,:) , lppe_v(i_opr_B,:) , lppe_w(i_opr_B,:) , lppe_p(i_opr_B,:)];

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(row_inds(:),:) =  lppe_factor * lppe_opr;
        B(row_inds(:),:) = -1i * lppe_opr;
    case Left_Side_Options.Symmetry
        opr_left_inds = get_operator_whole_left_indices(Nx, Ny);
        eqn_left_inds = get_equation_whole_left_indices('continuity', Nx, Ny);
        jm = (Nx - 1) / 2 + 1;
        opr = Z;
        for k = 1 : jm - 1
            left_ind = (jm + k - 1) * Ny + (1 : Ny);
            right_ind = (jm - k - 1) * Ny + (1 : Ny);
            for n = 1 : length(left_ind)
                opr(left_ind(n), left_ind(n)) = 1;
                opr(left_ind(n), right_ind(n)) = 1;
            end
        end

        symmetry_opr = [Z(opr_left_inds,:) Z(opr_left_inds,:) Z(opr_left_inds,:) opr(opr_left_inds,:)];
        
        A(eqn_left_inds,:) =  0;
        B(eqn_left_inds,:) =  0;
        A(eqn_left_inds,:) =  symmetry_factor * symmetry_opr;
        B(eqn_left_inds,:) = -1i * symmetry_opr;

        % opr_mid_inds = (jm - 1) * Ny + (1 : Ny) + (3 * Nx * Ny);
        % opr_mid_inds = opr_mid_inds(2:end);
        % 
        % row_inds    = opr_mid_inds;
        % column_inds = opr_mid_inds;
        % linear_inds = get_linear_indices(A, row_inds, column_inds);
        % 
        % A(row_inds(:),:) =  0;
        % B(row_inds(:),:) =  0;
        % A(linear_inds)   =  dirichlet_factor;
        % B(linear_inds)   = -1i;
    case Left_Side_Options.Anti_Symmetry
        opr_left_inds = get_operator_whole_left_indices(Nx, Ny);
        eqn_left_inds = get_equation_whole_left_indices('continuity', Nx, Ny);
        jm = (Nx - 1) / 2 + 1;
        opr = Z;
        for k = 1 : jm - 1
            left_ind = (jm + k - 1) * Ny + (1 : Ny);
            right_ind = (jm - k - 1) * Ny + (1 : Ny);
            for n = 1 : length(left_ind)
                opr(left_ind(n), left_ind(n)) = 1;
                opr(left_ind(n), right_ind(n)) = -1;
            end
        end

        symmetry_opr = [Z(opr_left_inds,:) Z(opr_left_inds,:) Z(opr_left_inds,:) opr(opr_left_inds,:)];
        
        A(eqn_left_inds,:) =  0;
        B(eqn_left_inds,:) =  0;
        A(eqn_left_inds,:) =  symmetry_factor * symmetry_opr;
        B(eqn_left_inds,:) = -1i * symmetry_opr;
        
        % opr_mid_inds = (jm - 1) * Ny + (1 : Ny);
        % opr_mid_inds = opr_mid_inds(2:end-1);
        % 
        % p_symmetry_opr = [Z(opr_mid_inds,:) Z(opr_mid_inds,:) Z(opr_mid_inds,:) Dx(opr_mid_inds,:)];
        % 
        % eqn_mid_inds = opr_mid_inds + (3 * Nx * Ny);
        % 
        % A(eqn_mid_inds,:) =  0;
        % B(eqn_mid_inds,:) =  0;
        % A(eqn_mid_inds,:) =  symmetry_factor * p_symmetry_opr;
        % B(eqn_mid_inds,:) = -1i * p_symmetry_opr;
    case Left_Side_Options.Neumann
        operator_row_inds = get_operator_left_indices(Nx, Ny);
        row_inds          = get_equation_left_indices('continuity', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        neumann_opr = [Z(operator_row_inds,:) Z(operator_row_inds,:) Z(operator_row_inds,:) Dx(operator_row_inds,:)];

        A(row_inds,:) =  0;
        B(row_inds,:) =  0;
        A(row_inds,:) =  neumann_factor * neumann_opr;
        B(row_inds,:) = -1i * neumann_opr;
    case Left_Side_Options.Dirichlet
        row_inds    = get_equation_left_indices('continuity', Nx, Ny);
        column_inds = get_variable_left_indices('p', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) =  0;
        B(row_inds(:),:) =  0;
        A(linear_inds)   =  dirichlet_factor;
        B(linear_inds)   = -1i;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Left.p ' for ''p'' at the left side is invalid or not supported'])
end

% Remove redundant or unnecessary entries
% Rows (equations at which flow variables are determined from boundary conditions)
xmom_top_row_inds    = get_equation_top_indices   ('x momentum', Nx, Ny);
ymom_top_row_inds    = get_equation_top_indices   ('y momentum', Nx, Ny);
zmom_top_row_inds    = get_equation_top_indices   ('z momentum', Nx, Ny);
cont_top_row_inds    = get_equation_top_indices   ('continuity', Nx, Ny);
xmom_bottom_row_inds = get_equation_bottom_indices('x momentum', Nx, Ny);
ymom_bottom_row_inds = get_equation_bottom_indices('y momentum', Nx, Ny);
zmom_bottom_row_inds = get_equation_bottom_indices('z momentum', Nx, Ny);
row_inds_to_remove = sort([xmom_top_row_inds ymom_top_row_inds zmom_top_row_inds cont_top_row_inds xmom_bottom_row_inds ymom_bottom_row_inds zmom_bottom_row_inds]);
A(row_inds_to_remove,:) = []; B(row_inds_to_remove,:) = [];
% Columns (variables which are determined from boundary conditions)
u_top_column_inds    = get_variable_top_indices   ('u', Nx, Ny);
v_top_column_inds    = get_variable_top_indices   ('v', Nx, Ny);
w_top_column_inds    = get_variable_top_indices   ('w', Nx, Ny);
p_top_column_inds    = get_variable_top_indices   ('p', Nx, Ny);
u_bottom_column_inds = get_variable_bottom_indices('u', Nx, Ny);
v_bottom_column_inds = get_variable_bottom_indices('v', Nx, Ny);
w_bottom_column_inds = get_variable_bottom_indices('w', Nx, Ny);
column_inds_to_remove = sort([u_top_column_inds v_top_column_inds w_top_column_inds p_top_column_inds u_bottom_column_inds v_bottom_column_inds w_bottom_column_inds]);
A(:,column_inds_to_remove) = []; B(:,column_inds_to_remove) = [];

end
