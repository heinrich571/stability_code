function plot_operators(A, B, Domain)

equation_list = {'x momentum' 'y momentum' 'z momentum' 'continuity'};
variable_list = {'u' 'v' 'w' 'p'};

xmat = Domain.mat_X;
ymat = Domain.mat_Y;
Ny = size(xmat, 1);
Nx = size(xmat, 2);

for i = 1 : length(equation_list)
    eqn = equation_list{i};
    
    for j = 1 : length(variable_list)
        var = variable_list{j};
        
        i_TR = get_eqn_top_right_ind(eqn, Nx, Ny);
        i_BL = get_eqn_bottom_left_ind(eqn, Nx, Ny);
        j_TR = get_var_top_right_ind(var, Nx, Ny);
        j_BL = get_var_bottom_left_ind(var, Nx, Ny);
        eqn_inds = i_TR : i_BL;
        var_inds = j_TR : j_BL;

        ALmat = reshape(A(i, var_inds), [Ny Nx]);
        BLmat = reshape(B(i, var_inds), [Ny Nx]);
        
        figure('Name', [eqn ' - L{' var '} '])
        subplot(1,2,1)
        surf(xmat, ymat, abs(ALmat))
        xlabel('x')
        ylabel('y')
        zlabel('L_A')
        subplot(1,2,2)
        surf(xmat, ymat, abs(BLmat))
        xlabel('x')
        ylabel('y')
        zlabel('L_B')
    end
end

end
