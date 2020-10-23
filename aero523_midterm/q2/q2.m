%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all; clc; close all
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
syms x
f1 = x^4; f2 = 2*x;

u1 = f1;
u2 = f2 - int(u1*f2, x, 0, 1)/int(u1*u1, x, 0, 1) * u1;

int(u1*u2, x, 0, 1)
tolatex('u1norm', u1/sqrt(int(u1*u1, x, 0, 1)))
tolatex('u2norm', u2/sqrt(int(u2*u2, x, 0, 1)))

g = x^2;
proju1g = int(u1*g, x, 0, 1) / int(u1*u1, x, 0, 1)* u1;
proju2g = int(u2*g, x, 0, 1) / int(u2*u2, x, 0, 1)* u2;

