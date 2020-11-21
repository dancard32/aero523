function [F,FL, FR, ls] = flux(Ul, Ur, n)
%FLUX Roe-Flux Function
gam = 1.4;    

% Left side arguments
rhol = Ul(1); ul = Ul(2)/rhol; vl = Ul(3)/rhol; rhoEl = Ul(4);
pl = (gam-1)*(rhoEl - 0.5*rhol*(ul^2 + vl^2));
Hl = (rhoEl + pl)/rhol;

% Right side arguments
rhor = Ur(1); ur = Ur(2)/rhor; vr = Ur(3)/rhor; rhoEr = Ur(4);
pr = (gam-1)*(rhoEr - 0.5*rhor*(ur^2 + vr^2));
Hr = (rhoEr + pr)/rhor;

% Left and Right side fluxes
FL = [dot(Ul(2:3), n), dot([Ul(2)*ul+pl, Ul(3)*ul],n), dot([Ul(2)*vl, Ul(3)*vl+pl],n), Hl*dot([Ul(2), Ul(3)],n)];
FR = [dot(Ur(2:3), n), dot([Ur(2)*ur+pr, Ur(3)*ur],n), dot([Ur(2)*vr, Ur(3)*vr+pr],n), Hr*dot([Ur(2), Ur(3)],n)];

% Roe-Averages
vell = [ul, vl]; velr = [ur, vr];
v = (sqrt(rhol)*vell + sqrt(rhor)*velr)/(sqrt(rhol) + sqrt(rhor));
H = (sqrt(rhol)*Hl + sqrt(rhor)*Hr)/(sqrt(rhol) + sqrt(rhor));

% Calculating Eigenvalues
q = norm(v);
c = sqrt((gam-1)*(H - 0.5*q^2));
u = dot(v, n);
ls = abs([u+c, u-c, u]);

% Entropy fix
for i=1:3
   if abs(ls(i))  < 0.1*c
       ls(i) = ((0.1*c)^2 + ls(i)^2)/(2*0.1*c);
   end
end

delrho = rhor - rhol; delmo = [rhor*ur-rhol*ul, rhor*vr-rhol*vl]; dele = rhoEr - rhoEl;
s1 = 0.5*(ls(1) + ls(2)); s2 = 0.5*(ls(1) - ls(2));
G1 = (gam-1)*(0.5*q^2*delrho - dot(v,delmo) + dele); G2 = -u*delrho + dot(delmo, n);
C1 = G1/c^2 *(s1 - ls(3)) + G2/c*s2; C2 = G1/c*s2 + (s1-ls(3))*G2;

RHS = [ls(3)*delrho+C1, ls(3)*delmo(1)+C1*v(1)+C2*n(1), ls(3)*delmo(2)+C1*v(2)+C2*n(2), ls(3)*dele+C1*H+C2*u];

F = 0.5*(FL + FR) - 0.5*RHS;
ls = max(ls);

end

