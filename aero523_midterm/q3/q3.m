%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all; clc; close all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Part b
N = 32; omega = 1.5;
D = eye(N+1); L = zeros(N+1); U = L;

for i = 1:N-1
   L(i+1, i) = -1; % Fill the left-matrix
end
for i = 2:N
    U(i, i+1) = -1;% Fill the right-matrix
    D(i,i) = 2; % Fill the rest of the diagonals
end

Sgs = -inv(D + L)*U; % Solve for the iterative matrix
lambdas = omega.*eig(Sgs) + (1-omega);

figure()
scatter(real(lambdas), imag(lambdas), 75, 'ko')
xlabel('Real axis, $\lambda(S_{GS})$', 'fontsize', 16)
ylabel('Imaginary axis, $\lambda(S_{GS})$', 'fontsize', 16)
set(gcf, 'Color', 'w', 'Position', [200 200 800 400]);
export_fig('q3_eigens.eps')
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Part c
num = 1000;
omegas = linspace(1, 2, num);
eigval = eig(Sgs);

data = zeros(num, 1);
for i = 1:num
    vals = omegas(i).*eigval + (1-omegas(i));
    norms = sqrt(real(vals).^2 + imag(vals).^2);
    data(i) = max(norms);
end

idx = find(data == min(data));

figure()
plot(omegas, data, 'k', 'linewidth', 2)
xlabel('Over-relaxation factor, $\omega$', 'fontsize', 16)
ylabel('Largest magnitude of each $\omega$', 'fontsize', 16)
legend(['$\omega = $', num2str(omegas(idx)), newline, '$\lambda_{min} = $', num2str(min(data))], 'location', 'best', 'fontsize', 16)
set(gcf, 'Color', 'w', 'Position', [200 200 800 400]);
export_fig('q3_omegas.eps')


