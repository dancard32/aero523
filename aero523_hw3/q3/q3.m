%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all; clc; close all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Input Variables
L = 2; a = 0.5; T = L/a;
nu = 0.1; sig = 0.5;
N = 64;

% Variables for input
xlin = linspace(0, L, N + 1);
dx = xlin(2);
dt = sig * dx/a;
nt = max(size(0:dt:T));

A = sparse(N + 1);  % Pre-allocate A matrix
for i = 1:(N+1)
   A(i,i)  = -2*nu/dx^2; % Diagonal
   
   if i > 1 % U, D entries
        A(i,i-1) = nu/dx^2 + a/(2*dx);
   end
   if i < N+1
        A(i,i+1) = nu/dx^2 - a/(2*dx);
   end
   if i == 1   % Handle Periodic Boundaries
       A(i,N+1) = nu/dx^2 + a/(2*dx);
   elseif i == N+1
       A(i,1) = nu/dx^2 - a/(2*dx);
   end
end
% Create super matrix for trapezoidal method
bigA = (eye(N+1) + dt/2*A)*inv((eye(N+1) - dt/2*A));

us = zeros(N+1, nt+1);  % Pre-allocation and initial condition
us(:,1) = exp(-100.*(xlin./L - 0.5).^2);
for j = 1:nt
    us(:,j+1) = bigA*us(:,j);   % Iterate with trapozoidal iteration
end

figure()
hold on
plot(xlin, us(:,end), 'k', 'linewidth', 2)
xlabel('$x$ axis', 'fontsize', 18)
ylabel('$u$ approximation', 'fontsize', 18)
legend('$u(x,T)$', 'location', 'best','fontsize', 18)
set(gcf, 'Color', 'w', 'Position', [200 200 800 400]);
export_fig('q3a.eps')
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% c
L = 2; a = 0.5; T = L/a;
nu = 0.1; sig = 0.5;
k = 4*pi/L; omega = 5*a/L;
% Variables for input
nums = [10, 64, 100, 500, 1000];
ntvals = zeros(1, max(size(nums))); 
l2err = zeros(1, max(size(nums))); 
l2errnt= zeros(1, max(size(nums))); idx = 1;
for N = nums
    fprintf('Iteration - %.f\n', N)
    xlin = linspace(0, L, N + 1);
    dx = xlin(2);
    dt = sig * dx/a;
    nt = N/a;
    tlin = linspace(0, T, nt + 1);
    ntvals(idx) = nt;

    us = zeros(N+1, nt+1);  % Pre-allocation and initial condition
    us(:,1) = sin(k.*xlin - omega*0);
    A = buildA(N, nu, dx, a);
    for j = 1:nt
        t = tlin(j);
        t2 = tlin(j+1);
        source1 = (a*k - omega).*cos(k.*xlin - omega*t) + nu*k^2*sin(k.*xlin - omega*t);
        source2 = (a*k - omega).*cos(k.*xlin - omega*t2) + nu*k^2*sin(k.*xlin - omega*t2);
        
        source_contrib = (eye(N+1) - dt/2*A)^(-1)*dt/2*(source1' + source2');
        us(:,j+1) = (eye(N+1) - dt/2*A)^(-1)*(eye(N+1) + dt/2*A)*us(:,j) + source_contrib;
    end

    uexact = sin(k.*xlin - omega*T);
    l2err(idx) = sqrt(1/N*sum((uexact' - us(:,end)).^2));
    l2errnt(idx) = sqrt(1/nt*sum((uexact' - us(:,end)).^2));
    idx = idx + 1;
    if N == 64
        figure()
        hold on
        plot(xlin, us(:,end), 'k-', 'linewidth', 2)
        plot(xlin, uexact, 'b', 'linewidth', 2)
        xlabel('$x$ axis', 'fontsize', 18)
        ylabel('$u$ approximation', 'fontsize', 18)
        legend({'$u(x,T)$', '$u^{MS}(x,T)$'}, 'location', 'northeast','fontsize', 18)
        set(gcf, 'Color', 'w', 'Position', [200 200 800 400]);
        export_fig('q3c.eps')
    end
end
ratesnx = zeros(1, (max(size(l2err)) - 1));
fid = fopen('rates_nx','w');
for i = 1:(max(size(l2err)) - 1)
    ratesnx(i) = abs(log10(l2err(i+1)/l2err(i))/log10(nums(i+1)/nums(i)));
    fprintf(fid,'$ \\text{rate}_{N = %.f} $ &  = %.4f\\\\ \n',nums(i), ratesnx(i));
end
fclose(fid);
ratesnt = zeros(1, (max(size(l2err)) - 1));
fid = fopen('rates_nt','w');
for i = 1:(max(size(l2err)) - 1)
    ratesnt(i) = abs(log10(l2errnt(i+1)/l2errnt(i))/log10(ntvals(i+1)/ntvals(i)));
    fprintf(fid,'$ \\text{rate}_{N = %.f} $ &  = %.4f\\\\ \n',ntvals(i), ratesnt(i));
end
fclose(fid);

figure()
hold on
plot(nums, l2err, 'k-', 'linewidth', 2)
plot(ntvals, l2errnt, 'b-', 'linewidth', 2)
xlabel('Number of intervals, $N_x,\ N_t$', 'fontsize', 18)
ylabel('$L_2$ Residual Error', 'fontsize', 18)
legend({['$N_x,\ L_2$ Solution Error', newline, 'Rate = ', num2str(ratesnx(1))], ...
    ['$N_t,\ L_2$ Solution Error', newline, 'Rate = ', num2str(ratesnt(1))]}, 'location', 'northeast','fontsize', 18)
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
set(gcf, 'Color', 'w', 'Position', [200 200 800 400]);
export_fig('q3d.eps')
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function A = buildA(N, nu, dx, a)
    A = sparse(N + 1);  % Pre-allocate A matrix
    for i = 1:(N+1)
       A(i,i)  = -2*nu/dx^2; % Diagonal

       if i > 1 % U, D entries
            A(i,i-1) = nu/dx^2 + a/(2*dx);
       end
       if i < N+1
            A(i,i+1) = nu/dx^2 - a/(2*dx);
       end

       if i == 1   % Handle Periodic Boundaries
           A(i,N+1) = nu/dx^2 + a/(2*dx);
       elseif i == N+1
           A(i,1) = nu/dx^2 - a/(2*dx);
       end
    end
end