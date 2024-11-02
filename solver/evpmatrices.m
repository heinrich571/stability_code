function [A, B] = evpmatrices(Domain, Base_Flow, Problem)

%% Build eigenvalue problem A and B matrices

[A, B] = bglnsematrices(Domain, Base_Flow, Problem);


%% Apply boundary conditions

Ny  = size(Domain.mat_X, 1);
Nx  = size(Domain.mat_X, 2);
Z   = zeros(size(Domain.Dx));
I   = eye(Ny*Nx);

mat_X = Domain.mat_X;
Dx    = Domain.Dx;
Dy    = Domain.Dy;
D2x   = Domain.D2x;
D2y   = Domain.D2y;

beta = Problem.Physics.Beta;

dirichlet_factor     = 200;
neumann_factor       = 300;
linear_extrap_factor = 400;
pressure_compatibility_factor = 500;
lppe_factor = 500;


%% Wall boundary conditions

% Multiple entry options allow for flexible user interaction with the code
Wall_Options.No_Slip        = {'Dirichlet' 'No_Slip' 'No-Slip' 'No Slip'}; % Meant for 'u' and 'w'
Wall_Options.No_Penetration = {'Dirichlet' 'No_Penetration' 'No-Penetration' 'No Penetration'}; % Meant for 'v'
Wall_Options.PC             = {'PC' 'Pressure-Compatibility' 'Pressure Compatibility' 'Pressure_Compatibility'}; % Meant for 'p'
Wall_Options.LPPE           = {'LPPE' 'Linearized-Pressure-Poisson-Equation' 'Linearized Pressure Poisson Equation' 'Linearized_Pressure_Poisson_Equation'}; % Meant for 'p'

% u
switch Problem.Boundary_Conditions.Wall.u
    case Wall_Options.No_Slip % No-slip boundary condition
        row_inds    = get_eqn_bottom_inds('x momentum', Nx, Ny);
        column_inds = get_var_bottom_inds('u', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) = 0;
        B(row_inds(:),:) = 0;
        A(linear_inds)   = dirichlet_factor;
        B(linear_inds)   = 0;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Wall.u ' for ''u'' at the wall is invalid or not supported'])
end

% v
switch Problem.Boundary_Conditions.Wall.v
    case Wall_Options.No_Penetration % No-penetration boundary condition
        row_inds    = get_eqn_bottom_inds('y momentum', Nx, Ny);
        column_inds = get_var_bottom_inds('v', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) = 0;
        B(row_inds(:),:) = 0;

        A(linear_inds)   = dirichlet_factor;
        B(linear_inds)   = 0;

        % % Enforcing dv/dy = 0 for the continuity equation
        % row_inds = get_eqn_bottom_inds('continuity', Nx, Ny);
        % i_opr_B  = get_opr_bottom_inds(Nx, Ny);
        % 
        % continuity_zero_dvdy_opr = [Z(i_opr_B,:) , Dy(i_opr_B,:) , Z(i_opr_B,:) , Z(i_opr_B,:)];
        % 
        % A(row_inds(:),:) = 0;
        % B(row_inds(:),:) = 0;
        % A(row_inds(:),:) = continuity_zero_dvdy_opr;
        % B(row_inds(:),:) = 0;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Wall.v ' for ''v'' at the wall is invalid or not supported'])
end

% w
switch Problem.Boundary_Conditions.Wall.w
    case Wall_Options.No_Slip % No-slip boundary condition
        row_inds    = get_eqn_bottom_inds('z momentum', Nx, Ny);
        column_inds = get_var_bottom_inds('w', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) = 0;
        B(row_inds(:),:) = 0;
        A(linear_inds)   = 1;
        B(linear_inds)   = 0;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Wall.w ' for ''w'' at the wall is invalid or not supported'])
end

% p
switch Problem.Boundary_Conditions.Wall.p
    case Wall_Options.PC % Pressure compatibility boundary condition
        row_inds = get_eqn_bottom_inds('continuity', Nx, Ny);
        i_opr_B  = get_opr_bottom_inds(Nx, Ny);

        pressure_compatibility_opr = [Z(i_opr_B,:) , -D2y(i_opr_B,:) , Z(i_opr_B,:) , Dy(i_opr_B,:)];

        A(row_inds(:),:) = 0;
        B(row_inds(:),:) = 0;
        A(row_inds(:),:) = pressure_compatibility_factor*pressure_compatibility_opr;
        B(row_inds(:),:) = 0;
    case Wall_Options.LPPE % LPPE boundary condition
        row_inds = get_eqn_bottom_inds('continuity', Nx, Ny);
        i_opr_B  = get_opr_bottom_inds(Nx, Ny);

        lppe_u = Z;
        lppe_v = Z;
        lppe_w = Z;
        lppe_p = D2x + D2y - beta^2*I;
        
        lppe_opr = [lppe_u(i_opr_B,:) , lppe_v(i_opr_B,:) , lppe_w(i_opr_B,:) , lppe_p(i_opr_B,:)];

        A(row_inds(:),:) = 0;
        B(row_inds(:),:) = 0;
        A(row_inds(:),:) = lppe_factor*lppe_opr;
        B(row_inds(:),:) = 0;
    otherwise
        warning('No boundary condition applied for the pressure on the wall')
end


%% Top boundary conditions

Top_Options.Decay = {'Dirichlet' 'Decay' 'decay'};
Top_Options.LPPE  = {'LPPE' 'Linearized-Pressure-Poisson-Equation' 'Linearized Pressure Poisson Equation'};

% u
switch Problem.Boundary_Conditions.Top.u
    case Top_Options.Decay % u(y-->infinity) = 0
        row_inds    = get_eqn_top_inds('x momentum', Nx, Ny);
        column_inds = get_var_top_inds('u', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) = 0;
        B(row_inds(:),:) = 0;
        A(linear_inds)   = dirichlet_factor;
        B(linear_inds)   = 0;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Wall.u ' for ''u'' at the top is invalid or not supported'])
end

% v
switch Problem.Boundary_Conditions.Top.v
    case Top_Options.Decay % v(y-->infinity) = 0
        row_inds    = get_eqn_top_inds('y momentum', Nx, Ny);
        column_inds = get_var_top_inds('v', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) = 0;
        B(row_inds(:),:) = 0;
        A(linear_inds)   = dirichlet_factor;
        B(linear_inds)   = 0;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Wall.v ' for ''v'' at the top is invalid or not supported'])
end

% w
switch Problem.Boundary_Conditions.Top.w
    case Top_Options.Decay % w(y-->infinity) = 0
        row_inds    = get_eqn_top_inds('z momentum', Nx, Ny);
        column_inds = get_var_top_inds('w', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) = 0;
        B(row_inds(:),:) = 0;
        A(linear_inds)   = dirichlet_factor;
        B(linear_inds)   = 0;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Wall.w ' for ''w'' at the top is invalid or not supported'])
end

% p
switch Problem.Boundary_Conditions.Top.p
    case Top_Options.Decay % p(y-->infinity) = 0
        row_inds    = get_eqn_top_inds('continuity', Nx, Ny);
        column_inds = get_var_top_inds('p', Nx, Ny);
        linear_inds = get_linear_indices(A, row_inds, column_inds);

        A(row_inds(:),:) = 0;
        B(row_inds(:),:) = 0;
        A(linear_inds)   = dirichlet_factor;
        B(linear_inds)   = 0;
    case Top_Options.LPPE % LPPE boundary condition
        row_inds = get_eqn_top_inds('continuity', Nx, Ny);
        i_opr_B  = get_opr_top_inds(Nx, Ny);

        lppe_u = Z;
        lppe_v = Z;
        lppe_w = Z;
        lppe_p = D2x + D2y - beta^2*I;

        lppe_opr = [lppe_u(i_opr_B,:) , lppe_v(i_opr_B,:) , lppe_w(i_opr_B,:) , lppe_p(i_opr_B,:)];

        A(row_inds(:),:) = 0;
        B(row_inds(:),:) = 0;
        A(row_inds(:),:) = lppe_factor*lppe_opr;
        B(row_inds(:),:) = lppe_opr;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Wall.p ' for ''p'' at the top is invalid or not supported'])
end


%% Right side

Right_Side_Options.FD_Extrapolation = {'Linear_Extrapolation' 'Linear-Extrapolation' 'Linear Extrapolation' 'FD_Extrapolation' 'FD Extrapolation' 'Finite Difference Extrapolation' 'finite_difference_extrapolation' 'finite difference extrapolation'};
Right_Side_Options.Zero_2nd_Derivative_Extrapolation = {'zero_2nd_derivative' 'zero_2nd_derivative_extrapolation'};
Right_Side_Options.LPPE = {'LPPE' 'Linearized-Pressure-Poisson-Equation' 'Linearized Pressure Poisson Equation'};

% u
switch Problem.Boundary_Conditions.Right.u
    case Right_Side_Options.FD_Extrapolation % Linear extrapolation
        j_u_right_inds    = get_var_right_inds('u', Nx, Ny);
        i_xmom_right_inds = get_eqn_right_inds('x momentum', Nx, Ny);

        A(i_xmom_right_inds(2:end-1),:) = 0;
        B(i_xmom_right_inds(2:end-1),:) = 0;
        
        for i = 2:Ny-1
            C = (mat_X(i,1)-mat_X(i,2))/(mat_X(i,3)-mat_X(i,2));
            linear_extrap_opr = [1 C-1 -C];
            ind_shift = Ny*[0 1 2];
            A(i_xmom_right_inds(i),j_u_right_inds(i) + ind_shift) = linear_extrap_opr;
            B(i_xmom_right_inds(i),j_u_right_inds(i) + ind_shift) = 0;
        end
    case Right_Side_Options.Zero_2nd_Derivative_Extrapolation % Zero 2nd derivative extrapolation
        operator_row_inds = get_opr_right_inds(Nx, Ny);
        row_inds          = get_eqn_right_inds('x momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        j_u_TR = get_var_top_right_ind('u', Nx, Ny);
        j_u_BL = get_var_bottom_left_ind('u', Nx, Ny);
        column_inds = j_u_TR:j_u_BL;

        A(row_inds,:) = 0;
        B(row_inds,:) = 0;
        A(row_inds,column_inds) = linear_extrap_factor*D2x(operator_row_inds,:);
        B(row_inds,column_inds) = 0;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Right.u ' for ''u'' at the right side is invalid or not supported'])
end

% v
switch Problem.Boundary_Conditions.Right.v
    case Right_Side_Options.FD_Extrapolation % Linear extrapolation
        j_v_right_inds    = get_var_right_inds('v', Nx, Ny);
        i_ymom_right_inds = get_eqn_right_inds('y momentum', Nx, Ny);
        
        A(i_ymom_right_inds(2:end-1),:) = 0;
        B(i_ymom_right_inds(2:end-1),:) = 0;
        
        for i = 2:Ny-1
            C = (mat_X(i,1)-mat_X(i,2))/(mat_X(i,3)-mat_X(i,2));
            linear_extrap_opr = [1 C-1 -C];
            ind_shift = Ny*[0 1 2];
            A(i_ymom_right_inds(i),j_v_right_inds(i) + ind_shift) = linear_extrap_opr;
            B(i_ymom_right_inds(i),j_v_right_inds(i) + ind_shift) = 0;
        end
    case Right_Side_Options.Zero_2nd_Derivative_Extrapolation % Zero 2nd derivative extrapolation
        operator_row_inds = get_opr_right_inds(Nx, Ny);
        row_inds          = get_eqn_right_inds('y momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        j_v_TR = get_var_top_right_ind('v', Nx, Ny);
        j_v_BL = get_var_bottom_left_ind('v', Nx, Ny);
        column_inds = j_v_TR:j_v_BL;

        A(row_inds,:) = 0;
        B(row_inds,:) = 0;
        A(row_inds,column_inds) = linear_extrap_factor*D2x(operator_row_inds,:);
        B(row_inds,column_inds) = 0;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Right.v ' for ''v'' at the right side is invalid or not supported'])
end

% w
switch Problem.Boundary_Conditions.Right.w
    case Right_Side_Options.FD_Extrapolation % Linear extrapolation
        j_w_right_inds    = get_var_right_inds('w', Nx, Ny);
        i_zmom_right_inds = get_eqn_right_inds('z momentum', Nx, Ny);
        
        A(i_zmom_right_inds(2:end-1),:) = 0;
        B(i_zmom_right_inds(2:end-1),:) = 0;
        
        for i = 2:Ny-1
            C = (mat_X(i,1)-mat_X(i,2))/(mat_X(i,3)-mat_X(i,2));
            linear_extrap_opr = [1 C-1 -C];
            ind_shift = Ny*[0 1 2];
            
            A(i_zmom_right_inds(i),j_w_right_inds(i) + ind_shift) = linear_extrap_opr;
            B(i_zmom_right_inds(i),j_w_right_inds(i) + ind_shift) = 0;
        end
    case Right_Side_Options.Zero_2nd_Derivative_Extrapolation % Zero 2nd derivative extrapolation
        operator_row_inds = get_opr_right_inds(Nx, Ny);
        row_inds          = get_eqn_right_inds('z momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        j_w_TR = get_var_top_right_ind('w', Nx, Ny);
        j_w_BL = get_var_bottom_left_ind('w', Nx, Ny);
        column_inds = j_w_TR:j_w_BL;

        A(row_inds,:) = 0;
        B(row_inds,:) = 0;
        A(row_inds,column_inds) = linear_extrap_factor*D2x(operator_row_inds,:);
        B(row_inds,column_inds) = 0;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Right.w ' for ''w'' at the right side is invalid or not supported'])
end

% p
switch Problem.Boundary_Conditions.Right.p
    case Right_Side_Options.FD_Extrapolation % Linear extrapolation
        j_p_right_inds    = get_var_right_inds('p', Nx, Ny);
        i_cont_right_inds = get_eqn_right_inds('continuity', Nx, Ny);
        
        A(i_cont_right_inds(2:end-1),:) = 0;
        B(i_cont_right_inds(2:end-1),:) = 0;

        for i = 2:Ny-1
            C = (mat_X(i,1)-mat_X(i,2))/(mat_X(i,3)-mat_X(i,2));
            linear_extrap_opr = [1 C-1 -C];
            ind_shift = Ny*[0 1 2];

            A(i_cont_right_inds(i),j_p_right_inds(i) + ind_shift) = linear_extrap_opr;
            B(i_cont_right_inds(i),j_p_right_inds(i) + ind_shift) = 0;
        end
    case Right_Side_Options.Zero_2nd_Derivative_Extrapolation % Zero 2nd derivative extrapolation
        operator_row_inds = get_opr_right_inds(Nx, Ny);
        row_inds          = get_eqn_right_inds('continuity', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        j_p_TR = get_var_top_right_ind('p', Nx, Ny);
        j_p_BL = get_var_bottom_left_ind('p', Nx, Ny);
        column_inds = j_p_TR:j_p_BL;

        A(row_inds,:) = 0;
        B(row_inds,:) = 0;
        A(row_inds,column_inds) = linear_extrap_factor*D2x(operator_row_inds,:);
        B(row_inds,column_inds) = 0;
    case Right_Side_Options.LPPE % LPPE boundary condition
        row_inds = get_eqn_right_inds('continuity', Nx, Ny);
        i_opr_B  = get_opr_right_inds(Nx, Ny);
        row_inds = row_inds(2:end-1);
        i_opr_B  = i_opr_B(2:end-1);

        lppe_u = 2*Ux*Dx;
        lppe_v = 2*(Uy*Dx + Vy*Dy);
        lppe_w = Z;
        lppe_p = D2x + D2y - beta^2*I;

        lppe_opr = [lppe_u(i_opr_B,:) , lppe_v(i_opr_B,:) , lppe_w(i_opr_B,:) , lppe_p(i_opr_B,:)];

        A(row_inds(:),:) = 0;
        B(row_inds(:),:) = 0;
        A(row_inds(:),:) = lppe_factor*lppe_opr;
        B(row_inds(:),:) = lppe_opr;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Right.p ' for ''p'' at the right side is invalid or not supported'])
end


%% Left side

Left_Side_Options.FD_Extrapolation = {'Linear_Extrapolation' 'Linear-Extrapolation' 'Linear Extrapolation' 'FD_Extrapolation' 'FD Extrapolation' 'Finite Difference Extrapolation' 'finite_difference_extrapolation' 'finite difference extrapolation'};
Left_Side_Options.Zero_2nd_Derivative_Extrapolation = {'zero_2nd_derivative' 'zero_2nd_derivative_extrapolation'};
Left_Side_Options.LPPE = {'LPPE' 'Linearized-Pressure-Poisson-Equation' 'Linearized Pressure Poisson Equation'};

% u
switch Problem.Boundary_Conditions.Left.u
    case Left_Side_Options.FD_Extrapolation % Linear extrapolation
        j_u_left_inds    = get_var_left_inds('u', Nx, Ny);
        i_xmom_left_inds = get_eqn_left_inds('x momentum', Nx, Ny);

        A(i_xmom_left_inds(2:end-1),:) = 0;
        B(i_xmom_left_inds(2:end-1),:) = 0;

        for i = 2:Ny-1
            C = (mat_X(i,Nx)-mat_X(i,Nx-1))/(mat_X(i,Nx-2)-mat_X(i,Nx-1));
            linear_extrap_opr = [1 C-1 -C];
            ind_shift = -Ny*[0 1 2];
            
            A(i_xmom_left_inds(i),j_u_left_inds(i) + ind_shift) = linear_extrap_opr;
            B(i_xmom_left_inds(i),j_u_left_inds(i) + ind_shift) = 0;
        end
    case Left_Side_Options.Zero_2nd_Derivative_Extrapolation % Zero 2nd derivative extrapolation
        operator_row_inds = get_opr_left_inds(Nx, Ny);
        row_inds          = get_eqn_left_inds('x momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        j_u_TR = get_var_top_right_ind('u', Nx, Ny);
        j_u_BL = get_var_bottom_left_ind('u', Nx, Ny);
        column_inds = j_u_TR:j_u_BL;

        A(row_inds,:) = 0;
        B(row_inds,:) = 0;
        A(row_inds,column_inds) = linear_extrap_factor*D2x(operator_row_inds,:);
        B(row_inds,column_inds) = 0;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Left.u ' for ''u'' at the left side is invalid or not supported'])
end

% v
switch Problem.Boundary_Conditions.Left.v
    case Left_Side_Options.FD_Extrapolation % Linear extrapolation
        j_v_left_inds    = get_var_left_inds('v', Nx, Ny);
        i_ymom_left_inds = get_eqn_left_inds('y momentum', Nx, Ny);

        A(i_ymom_left_inds(2:end-1),:) = 0;
        B(i_ymom_left_inds(2:end-1),:) = 0;
        for i = 2:Ny-1
            C = (mat_X(i,Nx)-mat_X(i,Nx-1))/(mat_X(i,Nx-2)-mat_X(i,Nx-1));
            linear_extrap_opr = [1 C-1 -C];
            ind_shift = -Ny*[0 1 2];
            
            A(i_ymom_left_inds(i),j_v_left_inds(i) + ind_shift) = linear_extrap_opr;
            B(i_ymom_left_inds(i),j_v_left_inds(i) + ind_shift) = 0;
        end
    case Left_Side_Options.Zero_2nd_Derivative_Extrapolation % Zero 2nd derivative extrapolation
        operator_row_inds = get_opr_left_inds(Nx, Ny);
        row_inds          = get_eqn_left_inds('y momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        j_v_TR = get_var_top_right_ind('v', Nx, Ny);
        j_v_BL = get_var_bottom_left_ind('v', Nx, Ny);
        column_inds = j_v_TR:j_v_BL;

        A(row_inds,:) = 0;
        B(row_inds,:) = 0;
        A(row_inds,column_inds) = linear_extrap_factor*D2x(operator_row_inds,:);
        B(row_inds,column_inds) = 0;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Left.v ' for ''v'' at the left side is invalid or not supported'])
end

% w
switch Problem.Boundary_Conditions.Left.w
    case Left_Side_Options.FD_Extrapolation % Linear extrapolation
        j_w_left_inds    = get_var_left_inds('w', Nx, Ny);
        i_zmom_left_inds = get_eqn_left_inds('z momentum', Nx, Ny);

        A(i_zmom_left_inds(2:end-1),:) = 0;
        B(i_zmom_left_inds(2:end-1),:) = 0;

        for i = 2:Ny-1
            C = (mat_X(i,Nx)-mat_X(i,Nx-1))/(mat_X(i,Nx-2)-mat_X(i,Nx-1));
            linear_extrap_opr = [1 C-1 -C];
            ind_shift = -Ny*[0 1 2];
            
            A(i_zmom_left_inds(i),j_w_left_inds(i) + ind_shift) = linear_extrap_opr;
            B(i_zmom_left_inds(i),j_w_left_inds(i) + ind_shift) = 0;
        end
    case Left_Side_Options.Zero_2nd_Derivative_Extrapolation % Zero 2nd derivative extrapolation
        operator_row_inds = get_opr_left_inds(Nx, Ny);
        row_inds          = get_eqn_left_inds('z momentum', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        j_w_TR = get_var_top_right_ind('w', Nx, Ny);
        j_w_BL = get_var_bottom_left_ind('w', Nx, Ny);
        column_inds = j_w_TR:j_w_BL;

        A(row_inds,:) = 0;
        B(row_inds,:) = 0;
        A(row_inds,column_inds) = linear_extrap_factor*D2x(operator_row_inds,:);
        B(row_inds,column_inds) = 0;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Left.w ' for ''w'' at the left side is invalid or not supported'])
end

% p
switch Problem.Boundary_Conditions.Left.p
    case Left_Side_Options.FD_Extrapolation % Linear extrapolation
        j_p_left_inds    = get_var_left_inds('p', Nx, Ny);
        i_cont_left_inds = get_eqn_left_inds('continuity', Nx, Ny);

        A(i_cont_left_inds(2:end-1),:) = 0;
        B(i_cont_left_inds(2:end-1),:) = 0;

        for i = 2:Ny-1
            C = (mat_X(i,Nx)-mat_X(i,Nx-1))/(mat_X(i,Nx-2)-mat_X(i,Nx-1));
            linear_extrap_opr = [1 C-1 -C];
            ind_shift = -Ny*[0 1 2];
            
            A(i_cont_left_inds(i),j_p_left_inds(i) + ind_shift) = linear_extrap_opr;
            B(i_cont_left_inds(i),j_p_left_inds(i) + ind_shift) = 0;
        end
    case Left_Side_Options.Zero_2nd_Derivative_Extrapolation % Zero 2nd derivative extrapolation
        operator_row_inds = get_opr_left_inds(Nx, Ny);
        row_inds          = get_eqn_left_inds('continuity', Nx, Ny);
        operator_row_inds = operator_row_inds(2:end-1); % exclude top and bottom parts of the domain, as boundary conditions there were already applied
        row_inds          = row_inds(2:end-1);

        j_p_TR = get_var_top_right_ind('p', Nx, Ny);
        j_p_BL = get_var_bottom_left_ind('p', Nx, Ny);
        column_inds = j_p_TR:j_p_BL;

        A(row_inds,:) = 0;
        B(row_inds,:) = 0;
        A(row_inds,column_inds) = linear_extrap_factor*D2x(operator_row_inds,:);
        B(row_inds,column_inds) = 0;
    case Left_Side_Options.LPPE % LPPE boundary condition
        row_inds = get_eqn_left_inds('continuity', Nx, Ny);
        i_opr_B  = get_opr_left_inds(Nx, Ny);
        row_inds = row_inds(2:end-1);
        i_opr_B  = i_opr_B(2:end-1);

        lppe_u = 2*Ux*Dx;
        lppe_v = 2*(Uy*Dx + Vy*Dy);
        lppe_w = Z;
        lppe_p = D2x + D2y - beta^2*I;

        lppe_opr = [lppe_u(i_opr_B,:) , lppe_v(i_opr_B,:) , lppe_w(i_opr_B,:) , lppe_p(i_opr_B,:)];

        A(row_inds(:),:) = 0;
        B(row_inds(:),:) = 0;
        A(row_inds(:),:) = lppe_factor*lppe_opr;
        B(row_inds(:),:) = lppe_opr;
    otherwise
        error(['Boundary condition ' Problem.Boundary_Conditions.Left.p ' for ''p'' at the left side is invalid or not supported'])
end


end


%% Supporting functions %%

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

equation_names  = {'x momentum' , 'y momentum' , 'z momentum' , 'continuity'};
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
%   i_opr_TR            row index of the operator matrix at the top-right corner of the domain      [vector]

i_opr_TR = 1;

end


function i_opr_BL = get_opr_bottom_left_ind(Nx, Ny)

% DESCRIPTION
%   This function returns the index of the operator matrix at the bottom-left corner of the domain.
% INPUT:
%   Nx                  Number of points along the chordwise direction                              [scalar]
%   Ny                  Number of points along the vertical direction                               [scalar]
% OUTPUT:
%   i_opr_BL            row index of the operator matrix at the bottom-left corner of the domain    [vector]

i_opr_BL = Nx*Ny;

end


function i_opr_BR = get_opr_bottom_right_ind(Nx, Ny)

% DESCRIPTION
%   This function returns the index of the operator matrix at the bottom-right corner of the domain.
% INPUT:
%   Nx                  Number of points along the chordwise direction                              [scalar]
%   Ny                  Number of points along the vertical direction                               [scalar]
% OUTPUT:
%   i_opr_BR            row index of the operator matrix at the bottom-right corner of the domain   [vector]

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
%   i_opr_TL            row index of the operator matrix at the top-left corner of the domain   [vector]

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
%   i_opr_T             row indices of the operator matrix at the top boundary of the domain    [vector]

i_opr_TR = get_opr_top_right_ind(Nx, Ny);
i_opr_BR = get_opr_top_left_ind(Nx, Ny);
i_opr_T  = i_opr_TR : Ny : i_opr_BR;

end


function i_opr_B = get_opr_bottom_inds(Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the operator matrix at the bottom boundary of the domain.
% INPUT:
%   Nx                  Number of points along the chordwise direction                              [scalar]
%   Ny                  Number of points along the vertical direction                               [scalar]
% OUTPUT:
%   i_opr_B             row indices of the operator matrix at the bottom boundary of the domain     [vector]

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
%   i_opr_R             row indices of the operator matrix at the right boundary of the domain  [vector]

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
%   i_opr_L             row indices of the operator matrix at the left boundary of the domain   [vector]

i_opr_TL = get_opr_top_left_ind(Nx, Ny);
i_opr_BL = get_opr_bottom_left_ind(Nx, Ny);
i_opr_L  = i_opr_TL : i_opr_BL;

end
