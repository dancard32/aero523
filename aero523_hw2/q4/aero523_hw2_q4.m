%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all; clc; close all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Read in values
E = dlmread('E.txt');
V = dlmread('V.txt');

fxsum = 0; fysum = 0; % Initialize the summations
for i = 1:length(E)
    node1 =  V(E(i,1),:);   % Nodal values
    node2 =  V(E(i,2),:);   % Nodal values
   
    lvec = node1 - node2;   % Direction of L
    deltal = norm(lvec);    % Calculate the total length of L
    midpoint = (node1 + node2)/2;      % Determine mid-point of L  
    
    nhat = lvec * [0, 1; -1, 0];    % Apply 90degree rotation
    nhat = nhat./norm(nhat);        % Normalize
    fxsum = fxsum + dot([midpoint(1), 0], nhat)*deltal; % Apply divergence
    fysum = fysum + dot([0, midpoint(2)], nhat)*deltal; % Apply divergence
end

% Output to latex
fid = fopen('fx_out','w');
fprintf(fid,'& = %.2f', fxsum);
fclose(fid);

fid = fopen('fy_out','w');
fprintf(fid,'& = %.2f', fysum);
fclose(fid);