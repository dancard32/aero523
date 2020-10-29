%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all; clc; close all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% a
syms u ut utt uttt utttt uttttt dt

unp1 = u + dt*ut + 1/2*dt^2*utt + 1/6*dt^3*uttt + 1/24*dt^4*utttt + 1/120*dt^5*uttttt;
un = u;
unm1 = u - dt*ut + 1/2*dt^2*utt - 1/6*dt^3*uttt + 1/24*dt^4*utttt - 1/120*dt^5*uttttt;
unm2 = u - 2*dt*ut + 2*dt^2*utt - 4/3*dt^3*uttt + 2/3*dt^4*utttt + (-2)^5/120*dt^5*uttttt;
fnp1 = ut + dt*utt + 1/2*dt^2*uttt + 1/6*dt^3*utttt + 1/24*dt^4*uttttt;

epsi = 5/3*unp1 - 5/2*un + unm1 - 1/6*unm2 - dt*fnp1;
pretty(simplify(epsi))
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% b
num = 100;
bdf2 = zeros(num, 1);
bdf3 = zeros(num, 1);
avg_der = zeros(num, 1);
thetlin = linspace(0, 2*pi, num);
g = exp(thetlin.*1i);
syms ldt
for i = 1:num
    eqn = g(i) == 4/3 - 1/(3*g(i)) + 2/3*ldt*g(i);
    sol = double(solve(eqn, ldt));
    bdf2(i) = sol;
    
    eqn = g(i) == 18/11 -9/(11*g(i)) + 2/(11*g(i)^2) + 6/11*ldt*g(i);
    sol = double(solve(eqn, ldt));
    bdf3(i) = sol;
    
    eqn = 5/3*g(i) - 5/2 + g(i)^(-1) - 1/6*g(i)^(-2) == ldt*g(i);
    sol = double(solve(eqn, ldt));
    avg_der(i) = sol;
end

figure()
hold on
fill(real(bdf3), imag(bdf3),[1,1,1],'facealpha', 1, 'FaceColor',[1,0,0],'EdgeColor','k','linewidth',1.8)
fill(real(avg_der), imag(avg_der),[1,1,1],'facealpha', 1, 'FaceColor',[0,0,1],'EdgeColor','k','linewidth',1.8)
fill(real(bdf2), imag(bdf2),[1,1,1],'facealpha', 1, 'FaceColor',[0.8,0.8,0.8],'EdgeColor','k','linewidth',1.8)
xlim([-2.5, 0])
axis equal
xlabel('$Re \left( \lambda \Delta T \right)$','fontsize', 18)
ylabel('$Im \left( \lambda \Delta T  \right)$','fontsize', 18)
legend({'BDF3','$\frac{1}{2}\frac{du}{dt}|_{BDF2} + \frac{1}{2}\frac{du}{dt}|_{BDF3}$','BDF2'}, 'fontsize', 18, 'location', 'best', 'interpreter', 'latex')
set(gcf, 'Color', 'w', 'Position', [100 100 1000 500]);
export_fig('eigs.eps')

syms ldt theta
g = exp(1i*theta);
eqn = 0 == 5/3 - 5/2*g^-1 + g^-2 - 1/6*g^-3;
sol = double(solve(eqn, theta));