%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File: solve_lp.m
% Author: David Fridovich-Keil ( dfk@eecs.berkeley.edu )
% 
% Sets up and solves the following LP, where P and h are loaded from disk.
%                         p* = arg max p' P' h
%                            s.t. p(i) >= 0
%                                 p' 1 == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% File names.
pzx_file = 'pzx_5x5_1000.csv';
hmz_file = 'hmz_5x5_1000.csv';

% Load matrices from file.
P = csvread(pzx_file);
h = csvread(hmz_file);

size(P)
size(h)

% Determine dimension of optimization variable.
n = size(P, 2);

% Setup and solve LP.
cvx_begin
    variable p(n)
    maximize transpose(h) * P * p
    subject to
        p >= 0
        sum(p) == 1    
cvx_end