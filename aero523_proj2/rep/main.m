%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all; clc; close all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'defaultaxesfontsize', 16);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


mesh = readgri('mesh0.gri');
V = mesh.Node; E = mesh.Elem; BE = mesh.B.nodes;
[IE, BE] = edgehash(E);

uinf = getIC(1, 2.2);
u0 = ones(mesh.nElem, 4).*uinf; u = u0;

R = zeros(mesh.nElem, 4); dta = zeros(mesh.nElem, 1); err = [];

figure()
fprintf('Residual Norm\n ------------------')
for iter = 1:10
   R = R.*0; dta = dta.*0;
    
   % Loop over Interior edges
   for i=1:max(size(IE))
       temp = IE(i,:);
       x1 = V(temp(1),:); x2 = V(temp(2),:);
       e1 = temp(3); e2 = temp(4);
       u1 = u(e1,:); u2 = u(e2,:);
       
       dx = x2 - x1; deltal = norm(dx);
       nhat = [-dx(2), dx(1)]./deltal;
       
       [F,FL, FR, ls] = flux(u1, u2, nhat);
       R(e1,:) = R(e1,:) + F*deltal;
       R(e2,:) = R(e2,:) - F*deltal;
       dta(e1) = dta(e1) + ls*deltal;
       dta(e2) = dta(e2) + ls*deltal;
   end
   
   % Loop over Boundary edges
   for i=1:max(size(BE))
       temp = IE(i,:);
       x1 = V(temp(1),:); x2 = V(temp(2),:);
       e1 = temp(3);
       uedge = u(e1,:);
       
       dx = x2 - x1; deltal = norm(dx);
       nhat = [dx(2), -dx(1)]./deltal;
       
       [F,FL, FR, ls] = flux(u1, u2, nhat);
       R(e1,:) = R(e1,:) - F*deltal;
       dta(e1) = dta(e1) + ls*deltal;
   end
   
   CFL = 1;
   dta = 2*CFL./dta;
   u = u - dta.*R;
   err(end+1) = sum(sum(abs(R))); err(end)
   
   P = (1.4-1)*(u(:,4) - 0.5.*u(:,1).*( (u(:,2)./u(:,1)).^2 + (u(:,3)./u(:,1)).^2 ));
   
   pdeplot(V, E, 'XYData', P);
   axis equal; axis tight; axis off;
   colormap jet
end