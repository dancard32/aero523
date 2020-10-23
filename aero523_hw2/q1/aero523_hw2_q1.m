%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all; clc; close all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

syms x_0 x_p1 x_m1 y_0 y_p1 y_m1 Delta_x Delta_y x y u_m1p1 u_m1m1 u_00 u_p1p1 u_p1m1

x_nodep1 = ((x - x_m1) * (x - x_0)) / ((2*Delta_x) * (Delta_x));
x_nodem1 = ((x - x_0) * (x - x_p1)) / ( (-Delta_x) * (-2*Delta_x));

y_nodep1 = ((y - y_m1) * (y - y_0)) / ((2*Delta_y) * (Delta_y));
y_nodem1 = ((y - y_0) * (y - y_p1)) / ( (-Delta_y) * (-2*Delta_y));

pretty(simplify(diff(diff(x_nodep1*y_nodep1, x), y))) % xp1yp1
pretty(simplify(diff(diff(x_nodem1*y_nodem1, x), y))) % xm1ym1
pretty(simplify(diff(diff(x_nodep1*y_nodem1, x), y))) % xp1ym1
pretty(simplify(diff(diff(x_nodem1*y_nodep1, x), y))) % xm1yp1