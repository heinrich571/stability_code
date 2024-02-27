function [fig_baseflow_simple, fig_baseflow_with_streamlines] = show_baseflow(Domain, Base_Flow)

fig_baseflow_simple = figure('Name', 'Base flow solution', 'NumberTitle', 'off');
plot(Base_Flow.phi, Base_Flow.eta  , 'DisplayName', '$\phi(\eta)$')
plot(Base_Flow.dphi, Base_Flow.eta , 'DisplayName', '$\phi''(\eta)$')
plot(Base_Flow.ddphi, Base_Flow.eta, 'DisplayName', '$\phi''''(\eta)$')
ylabel('$\eta = \sqrt{\frac{a}{\nu}} \cdot y$', 'Interpreter', 'LaTeX')
ylim([min(Domain.vec_Y) max(Domain.vec_Y)])
legend('Interpreter', 'LaTeX', 'Location', 'southeast')


fig_baseflow_with_streamlines = figure('Name', 'Base flow with streamlines', 'NumberTitle', 'off');

% plot(Domain.mat_X, Domain.mat_Y, 'k', 'LineWidth', 1, 'Color', 0.9*ones(1, 3))
% plot(Domain.mat_X', Domain.mat_Y', 'k', 'LineWidth', 1, 'Color', 0.9*ones(1, 3))

U      =  Domain.mat_X.*repmat(flip(Base_Flow.dphi), [1 size(Domain.mat_X, 2)]);
V      = -repmat(flip(Base_Flow.phi), [1 size(Domain.mat_X, 2)]);
start_X = min(Domain.vec_X) : 0.25 : -0.25;
start_X = [start_X -flip(start_X)];
start_X = start_X(start_X ~= 0);
start_Y = max(Domain.vec_Y)*ones(size(start_X));

lines = streamline(Domain.mat_X, Domain.mat_Y, U, V, start_X, start_Y);
for i = 1:length(lines)
    lines(i).LineWidth = 1;
    lines(i).Color = [000, 090, 255]/255;
end
lines(1).DisplayName = 'Streamlines';
p_dividing_streamline = plot(zeros(size(Domain.vec_Y)), Domain.vec_Y, '--', 'LineWidth', 1, 'Color', lines(1).Color, 'DisplayName', 'Dividing streamline');
p_wall = plot(Domain.vec_X, zeros(size(Domain.vec_X)), '-k', 'DisplayName', 'Wall');
p_stagnation_point = plot(0, 0, 'ro', 'MarkerSize', 10, 'DisplayName', 'Stagnation point');
delta99 = interp1(Base_Flow.dphi, Base_Flow.eta, 0.99);
p_boundary_layer_height = plot(Domain.vec_X, delta99*ones(size(Domain.vec_X)), '--', 'LineWidth', 1, 'Color', p_stagnation_point.Color);
eta_for_bl_plot = [delta99 ; Base_Flow.eta(Base_Flow.eta < delta99) ; delta99];
bl_for_plot = [Base_Flow.dphi(1) ; Base_Flow.dphi(Base_Flow.eta < delta99) ; interp1(Base_Flow.eta, Base_Flow.dphi, delta99)];
p_boundary_layer_profile_right = plot(bl_for_plot + max(Domain.vec_X)/3, eta_for_bl_plot, 'LineWidth', 1, 'Color', p_stagnation_point.Color, 'DisplayName', 'Boundary layer');
p_boundary_layer_profile_left  = plot(-(bl_for_plot + max(Domain.vec_X)/3), eta_for_bl_plot, 'LineWidth', 1, 'Color', p_stagnation_point.Color);
x_arrows = max(Domain.vec_X)/3*ones(7, 1);
y_arrows = (0:delta99/(length(x_arrows)-1):delta99)';
u_arrows = interp1(Base_Flow.eta, Base_Flow.dphi, y_arrows);
v_arrows = zeros(size(x_arrows));
quiver( x_arrows, y_arrows,  u_arrows, v_arrows, 'LineWidth', 1, 'Color', p_stagnation_point.Color, 'AutoScale', 'off', 'MaxHeadSize', 0.2);
quiver(-x_arrows, y_arrows, -u_arrows, v_arrows, 'LineWidth', 1, 'Color', p_stagnation_point.Color, 'AutoScale', 'off', 'MaxHeadSize', 0.2);
xlabel('x')
ylabel('y')
legend([p_wall lines(1) p_stagnation_point p_boundary_layer_profile_right])
axis equal
grid off
box on

end
