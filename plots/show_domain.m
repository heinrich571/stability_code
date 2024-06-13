function fig = show_domain(Domain)

fig = figure('Name', 'Domain', 'NumberTitle', 'off');
plot(Domain.mat_X, Domain.mat_Y, 'k', 'LineWidth', 1)
plot(Domain.mat_X', Domain.mat_Y', 'k', 'LineWidth', 1)
xlabel('$x$')
ylabel('$y$')
axis equal

end
