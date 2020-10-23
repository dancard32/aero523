%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all; clc; close all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

num = 50;
xilin = linspace(0, 1, num);
[xilin, etalin] = meshgrid(xilin, xilin);

xlin = zeros(num);
ylin = zeros(num);
for i = 1:num
    for j = 1:num
        xlin(j,i) = xilin(j,i) + 1/2*etalin(j,i)^2;
        ylin(j,i) = 1+ etalin(j,i) - 1/4*xilin(j,i)*etalin(j,i);
    end
end



figure()
hold on
plot(xlin(1:num, 1),ylin(1:num, 1), 'k', 'linewidth', 1.8)
plot(xlin(1:num, num),ylin(1:num, num), 'k', 'linewidth', 1.8)
plot(xlin(1, 1:num),ylin(1, 1:num), 'k', 'linewidth', 1.8)
plot(xlin(num, 1:num),ylin(num, 1:num), 'k', 'linewidth', 1.8)
xlabel('X-Axis', 'fontsize', 16)
ylabel('Y-Axis', 'fontsize', 16)
set(gcf, 'Color', 'w', 'Position', [200 200 800 400]);
%export_fig('domain.eps')


syms x y xi eta
eqn1 = x == xi + 1/2*eta^2;
eqn2 = y == 1 + eta - 1/4*xi*eta;
sol = solve([eqn1, eqn2],[xi, eta], 'maxdegree', 8);
xi = sol.xi(1);
eta = sol.eta(1);

n = zeros(4,2);
n_eta = [diff(eta, x), diff(eta, y)];
n(3,:) = double(subs(n_eta, [x, y],[0.5+1/2*1^2, 1+1-1/4*0.5*1]));
n(1,:) = double(subs(-n_eta, [x, y],[0.5+1/2*0^2, 1+0-1/4*0.5*0]));

n_xi = [diff(xi, x), diff(xi, y)];
n(2,:) = double(subs(n_xi, [x, y],[1+1/2*0.5^2, 1+0.5-1/4*0.5*1]));
n(4,:) = double(subs(-n_xi, [x, y],[0+1/2*0.5^2, 1+0.5-1/4*0.5*0]));

n(abs(n) < 10^(-10)) = 0; % Correct for output to LaTeX
fid = fopen('n_output','w');
for i = 1:4
    n(i,:) = n(i,:)./norm(n(i,:));
    string = append('n_', num2str(i), ' & = [',num2str(n(i,1)), ', ', num2str(n(i,2)),"] \\");
    fprintf(fid,'%s \n', string);
end
fclose(fid);




















