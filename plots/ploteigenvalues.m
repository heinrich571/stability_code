function ploteigenvalues(Solution)

figure('Name', 'Eigenvalues', 'NumberTitle', 'off')
scatter(real(Solution.omega), imag(Solution.omega), 'o')
xlabel('Re[\omega]', 'Interpreter', 'latex')
ylabel('Im[\omega]', 'Interpreter', 'latex')

end
