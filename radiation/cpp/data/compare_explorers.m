%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Author: David Fridovich-Keil
% File: compare_explorers.m
%
% Usage: Set the file names for the RW and LP recorded data and run.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

lp_data = csvread('size5x5_srcs2_steps2_samples20k_lp.csv');
rw_data = csvread('size5x5_srcs2_steps2_samples20k_rw.csv');

num_iterations = size(lp_data, 2);

lp_means = zeros(1, num_iterations); lp_stddevs = zeros(1, num_iterations);
rw_means = zeros(1, num_iterations); rw_stddevs = zeros(1, num_iterations);

for jj = 1:num_iterations
    lp_means(jj) = mean(lp_data(:, jj)); 
    lp_stddevs(jj) = std(lp_data(:, jj));
    rw_means(jj) = mean(rw_data(:, jj));
    rw_stddevs(jj) = std(rw_data(:, jj));
end

figure;
hold on;
grid on;
errorbar((0:(num_iterations-1)) - 0.1, lp_means, 1.0 * lp_stddevs, 'r*');
errorbar((0:(num_iterations-1)) + 0.1, rw_means, 1.0 * rw_stddevs, 'bo');
xlim([-0.5, num_iterations-0.5]);
set(gca, 'fontsize', 24);
legend('LP explorer', 'RW explorer');
xlabel('Time step');
ylabel('Entropy (nats)');
title(sprintf('Empirical Comparison of Information Maximization \n and a Random Walk'));