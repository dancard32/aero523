%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all; clc; close all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
num = 200;
ntlin = 10.^(2:4);
l2norms = zeros(max(size(ntlin)), num);
intervals = zeros(max(size(ntlin)), num);
for j = 1:max(size(ntlin))
    intervals(j,:) = floor(linspace(10, 2*ntlin(j), num));
    for i = 1:max(size(intervals))
        [us, xlin] = BW_method(intervals(j,i), ntlin(j));

        l2norms(j,i) = sqrt(1./ntlin(j).*sum((us(1,:) - us(end,:)).^2));
    end
end

val = floor(num/3);
rate1 = abs(log10(l2norms(1,2*val)/l2norms(1,val))/log10(intervals(1,2*val)/intervals(1,val)));
rate2 = abs(log10(l2norms(2,floor(2.5*val))/l2norms(2,floor(1.1*val)))/log10(intervals(2,floor(2.5*val))/intervals(2,floor(1.1*val))));
rate3 = abs(log10(l2norms(3,floor(2.4*val))/l2norms(3,floor(1.1*val)))/log10(intervals(3,floor(2.4*val))/intervals(3,floor(1.1*val))));

figure()
hold on
plot(intervals(1,:), l2norms(1,:), 'k', 'linewidth', 2)
plot(intervals(2,:), l2norms(2,:), 'g', 'linewidth', 2)
plot(intervals(3,:), l2norms(3,:), 'b', 'linewidth', 2)
xlabel('$N_x$ Intervals', 'fontsize', 16)
ylabel('$L_2$ Residual Error Norm', 'fontsize', 16)
legend({['$N_t$ = 100: $\mathcal{O}$(', num2str(rate1),')'],...
        ['$N_t$ = 1000: $\mathcal{O}$(', num2str(rate2),')'],...
        ['$N_t$ = 10000: $\mathcal{O}$(', num2str(rate3),')']}, 'fontsize', 16, 'location', 'northeast')
set(gca, 'yscale', 'log')
set(gca, 'xscale', 'log')
set(gcf, 'Color', 'w', 'Position', [200 200 800 400]);
export_fig('BW_convergence.eps')
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [us, xlin] = BW_method(nx, nt)
    L = 2; a = 0.5; T = L/a;
    xlin = linspace(0, L, nx + 1);
    tlin = linspace(0, T, nt + 1);
    us = zeros(nt+1, nx+1);
    
    us(1,:) = exp(-100.*(xlin./L - 0.5).^2);
    dx = xlin(2); dt = tlin(2);
    sig = a*dt/dx;
    
    for n = 1:nt
        for j = 1:(nx+1)
            if j-1 == 0
                ujm1 = us(n,nx);
                ujm2 = us(n,nx-1);
            elseif j-1 == 1
                ujm1 = us(n,j-1);
                ujm2 = us(n,nx);
            else
                ujm1 = us(n,j-1);
                ujm2 = us(n,j-2);    
            end

            uj = us(n,j);
            us(n+1,j) = us(n,j) - sig/2*(3*uj - 4*ujm1 + ujm2) + sig^2/2*(ujm2 - 2*ujm1 + uj);
        end
    end
end

