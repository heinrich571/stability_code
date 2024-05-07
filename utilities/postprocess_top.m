eigenvalue_number = 274;
w_name = num2str(evalue(eigenvalue_number));

Nx_ = Nx + 1;
Ny_ = Ny + 1;

u_inds = 1 : Nx_*Ny_;
v_inds = Nx_*Ny_+1 : 2*Nx_*Ny_;
w_inds = 2*Nx_*Ny_+1 : 3*Nx_*Ny_;
p_inds = 3*Nx_*Ny_+1 : 4*Nx_*Ny_;

u = reshape(eigenfunctions_matrix(u_inds,eigenvalue_number), Ny_, Nx_);
v = reshape(eigenfunctions_matrix(v_inds,eigenvalue_number), Ny_, Nx_);
w = reshape(eigenfunctions_matrix(w_inds,eigenvalue_number), Ny_, Nx_);
p = reshape(eigenfunctions_matrix(p_inds,eigenvalue_number), Ny_, Nx_);

figure; contourf(Domain.mat_X, Domain.mat_Y, abs(u), 'linestyle', 'none'); title(['eigenvalue: ' w_name])
colorbar
figure; contourf(Domain.mat_X, Domain.mat_Y, abs(v), 'linestyle', 'none'); title(['eigenvalue: ' w_name])
colorbar
figure; contourf(Domain.mat_X, Domain.mat_Y, abs(w), 'linestyle', 'none'); title(['eigenvalue: ' w_name])
colorbar
figure; contourf(Domain.mat_X, Domain.mat_Y, abs(p), 'linestyle', 'none'); title(['eigenvalue: ' w_name])
colorbar