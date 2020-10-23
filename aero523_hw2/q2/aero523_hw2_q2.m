%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all; clc; close all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

syms a_3 a_2 a_1 a_0 h

eqn1 = 0 == a_2 + a_1 + a_0;
eqn2 = 0 == 2*a_2 + a_1;
eqn3 = 1/h^2 == 2*a_2 + 1/2*a_1;
sol = solve([eqn1, eqn2, eqn3],[a_2, a_1, a_0]);

tolatex('p1_a2',sol.a_2)
tolatex('p1_a1',sol.a_1)
tolatex('p1_a0',sol.a_0)

eqn1 = 0 == a_3 + a_2 + a_1 + a_0;
eqn2 = 0 == 3*a_3 + 2*a_2 + a_1;
eqn3 = 1/h^2 == 9/2*a_3 + 2*a_2 + 1/2*a_1;
eqn4 = 0 == 9/2*a_3 + 4/3*a_2 + 1/6*a_1;
sol = solve([eqn1, eqn2, eqn3, eqn4],[a_3, a_2, a_1, a_0]);

tolatex('p2_a3',sol.a_3)
tolatex('p2_a2',sol.a_2)
tolatex('p2_a1',sol.a_1)
tolatex('p2_a0',sol.a_0)