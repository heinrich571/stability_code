function Figures = profileAtX(Solution, eigenvalue_inds, X_Stations, Legend_Entry, User_Figures)

if nargin < 5
    User_Figures = figure();
end

if nargin < 4
    Legend_Entry = '';
end

mat_X = Solution.Domain.mat_X;
mat_Y = Solution.Domain.mat_Y;

vec_X = Solution.Domain.vec_X;
vec_Y = Solution.Domain.vec_Y;

Nx = length(vec_X);
Ny = length(vec_Y);

N_Stations = length(X_Stations);

Figures = User_Figures;

count = 1;

for i = eigenvalue_inds
    u_mat = abs(reshape(Solution.Eigenfunctions.u(:,i), Ny, Nx));
    v_mat = abs(reshape(Solution.Eigenfunctions.v(:,i), Ny, Nx));
    w_mat = abs(reshape(Solution.Eigenfunctions.w(:,i), Ny, Nx));
    p_mat = abs(reshape(Solution.Eigenfunctions.p(:,i), Ny, Nx));
    min_u = min(u_mat, [], 'all');
    max_u = max(u_mat, [], 'all');
    min_v = min(v_mat, [], 'all');
    max_v = max(v_mat, [], 'all');
    min_w = min(w_mat, [], 'all');
    max_w = max(w_mat, [], 'all');
    min_p = min(p_mat, [], 'all');
    max_p = max(p_mat, [], 'all');
    for j = 1:4
        if isempty(User_Figures(1).Children)
            switch j
                case 1
                    figureName = ['omega_' num2str(i) ', u - X Stations Profile'];
                case 2
                    figureName = ['omega_' num2str(i) ', v - X Stations Profile'];
                case 3
                    figureName = ['omega_' num2str(i) ', w - X Stations Profile'];
                case 4
                    figureName = ['omega_' num2str(i) ', p - X Stations Profile'];
                otherwise
                    error('Something went wrong with indexing...')
            end
            Figures(count,j) = figure('Name', figureName, 'NumberTitle', 'off');
        else
            figure(Figures(count,j));
        end
        sgtitle(['$\omega_' num2str(i) ' = ' num2str(Solution.Eigenvalues(i)) '$'], 'fontsize', 24)
        for k = 1:N_Stations
            subplot(1, N_Stations, k)
            title(['$x = ' num2str(X_Stations(k)) '$'])
            ylabel('$y$')
            switch j
                case 1
                    min_var = min_u;
                    max_var = max_u;
                    xlabel('$\| u \|$')
                    var = interp2(mat_X, mat_Y, u_mat, X_Stations(k)*ones(size(vec_Y)), vec_Y);
                case 2
                    min_var = min_v;
                    max_var = max_v;
                    xlabel('$\| v \|$')
                    var = interp2(mat_X, mat_Y, v_mat, X_Stations(k)*ones(size(vec_Y)), vec_Y);
                case 3
                    min_var = min_w;
                    max_var = max_w;
                    xlabel('$\| w \|$')
                    var = interp2(mat_X, mat_Y, w_mat, X_Stations(k)*ones(size(vec_Y)), vec_Y);
                case 4
                    min_var = min_p;
                    max_var = max_p;
                    xlabel('$\| p \|$')
                    var = interp2(mat_X, mat_Y, p_mat, X_Stations(k)*ones(size(vec_Y)), vec_Y);
                otherwise
                    error('Something went wrong with indexing...')
            end
            plot(var, vec_Y, 'DisplayName', Legend_Entry)
            xlim([min_var max_var])
        end
    end
    count = count + 1;
end

end
