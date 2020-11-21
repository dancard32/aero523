function uinf = getIC(alpha, mach)
    gam = 1.4;
    uinf = [1, mach*cosd(alpha), mach*sin(alpha), 1/(gam*(gam-1)) + mach^2/2];
    
end