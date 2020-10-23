%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all; clc; close all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%set(gcf, 'Color', 'w', 'Position', [200 200 800 400]);
%export_fig('figs/test.eps')
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% num = 50, Nt = 50     % sig = 1
% Nt = 50; num = 2*Nt - 1;% sig = 2
num = 50; Nt = 500;
L = 2; a = 0.5; T = L/a;
xlin = linspace(0, 2, num);
u0 = exp(-100.*(xlin./L - 0.5).^2);
us = zeros(Nt, num);
us(1,:) = u0;

dx = L/(num-1);
dt = T/(Nt -1);
sig = a*dt/dx;


figure()
pause(5)
hold on
for n = 1:Nt
    dt*(n-1)
    plot(xlin, us(n,:))
    pause(0.005)
    if n < Nt-1
        for j = 3:num
            uj = us(n,j);
            ujm1 = us(n,j-1);
            ujm2 = us(n,j-2);
            us(n+1,j) = us(n,j) - sig/2*(3*uj - ujm1 + ujm2) + sig^2/2*(ujm2 - 2*ujm1 + uj);
        end
    end
end

