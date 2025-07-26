function [omega, u, v, w, p, mat_X, mat_Y, Nx_full, Ny] = parse_solution(Solution, si, evi)

Nx = length(Solution(si).Domain.vec_X);
Ny = length(Solution(si).Domain.vec_Y);

mat_X = Solution(si).Domain.mat_X;
mat_Y = Solution(si).Domain.mat_Y;

omega = Solution(si).Eigenvalues(evi);
u = reshape(Solution(si).Eigenfunctions.u(:,evi), Ny, Nx, length(evi));
v = reshape(Solution(si).Eigenfunctions.v(:,evi), Ny, Nx, length(evi));
w = reshape(Solution(si).Eigenfunctions.w(:,evi), Ny, Nx, length(evi));
p = reshape(Solution(si).Eigenfunctions.p(:,evi), Ny, Nx, length(evi));

if mat_X(1,end) == 0 % If solution is for half a domain, unfold it to the full domain
    Nx_full = 2 * Nx - 1;
    
    for iR = Nx-1 : -1 : 1
        iL = Nx_full - iR + 1;

        mat_X(:,iL) = -mat_X(:,iR);
        mat_Y(:,iL) =  mat_Y(:,iR);
    end

    for n = 1 : length(evi)
        for iR = Nx-1 : -1 : 1
            iL = Nx_full - iR + 1;

            if abs(u(end,end,n)) < 1e-6 % Anti-symmetric mode
                u(:,iL) = -u(:,iR,n);
                v(:,iL) =  v(:,iR,n);
                w(:,iL) =  w(:,iR,n);
                p(:,iL) =  p(:,iR,n);
            else % Symmetric mode
                u(:,iL) =  u(:,iR,n);
                v(:,iL) = -v(:,iR,n);
                w(:,iL) = -w(:,iR,n);
                p(:,iL) = -p(:,iR,n);
            end
        end
    end
end

end
