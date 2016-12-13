%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Author: David Fridovich-Keil
% File: compare_explorers.m
%
% Usage: Set the file names for the RW and LP recorded data and run.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

rw_data = csvread('size10x10_srcs2_steps1_samples1k_rw.csv');
lp_data = csvread('size10x10_srcs2_steps1_samples1k_lp.csv');
socp_data = csvread('size10x10_srcs2_steps1_samples1k_socp.csv');

num_iterations = size(lp_data, 2);

socp_means = zeros(1, num_iterations); socp_stddevs = zeros(1, num_iterations);
lp_means = zeros(1, num_iterations); lp_stddevs = zeros(1, num_iterations);
rw_means = zeros(1, num_iterations); rw_stddevs = zeros(1, num_iterations);

for jj = 1:num_iterations
    socp_means(jj) = mean(socp_data(:, jj)); 
    socp_stddevs(jj) = std(socp_data(:, jj));
    lp_means(jj) = mean(lp_data(:, jj)); 
    lp_stddevs(jj) = std(lp_data(:, jj));
    rw_means(jj) = mean(rw_data(:, jj));
    rw_stddevs(jj) = std(rw_data(:, jj));
end

figure;
hold on;
grid on;
errorbar((0:(num_iterations-1)), socp_means, 1.0 * socp_stddevs, 'kX');
errorbar((0:(num_iterations-1)) - 0.1, lp_means, 1.0 * lp_stddevs, 'r*');
errorbar((0:(num_iterations-1)) + 0.1, rw_means, 1.0 * rw_stddevs, 'bo');
xlim([-0.5, num_iterations-0.5]);
set(gca, 'fontsize', 24);
legend('SOCP explorer', 'LP explorer', 'RW explorer');
xlabel('Time step');
ylabel('Entropy (nats)');
title(sprintf('Empirical Comparison of (Robust) Information Maximization \n and a Random Walk'));