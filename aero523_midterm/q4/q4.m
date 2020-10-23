%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all; clc; close all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N = [4, 5, 7, nan, nan, nan; 7, 8, 9, nan, nan, nan;
     9, 10, nan, nan, nan, nan; 10,11,13, nan, nan, nan;
     13,16, nan, nan, nan, nan; 15,16, nan, nan, nan, nan;
     14,15, nan, nan, nan, nan; 12,14, nan, nan, nan, nan;
     8:12, nan; 5:8, nan, nan;
     2,6, nan, nan, nan, nan; 1,2, nan, nan, nan, nan;
     1,3, nan, nan, nan, nan; 3,4, nan, nan, nan, nan;
     1:6; 11:16];
E = [12, 13, 15; 11, 12, 15; 13, 14, 15; 1, 15, 14; 1, 10, 15;
    10, 11, 15; 1, 2 , 10; 2, 9 , 10; 2, 3 , 9; 3, 4 , 9;
    4, 16,  9; 8, 9 , 16; 4, 5 , 16; 7, 8 , 16; 6, 7 , 16; 5, 6 , 16];

counter = 0;
for j = 1:max(size(N))
    for i = 1:max(size(N(j, :)))
        if isnan(N(j, i)) ~= 1
            [row, col] = find(N == N(j, i));
            N(row(end), col(end)) = nan;
            counter = counter + 1;
        end
    end
end
counter = counter - 1;
    
