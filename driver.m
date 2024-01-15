clc; clear all; close all; 

% Add the path to the class file if it's not in the current directory
% addpath('path_to_class_file');

% Define the function and parameters for PCE
func = @(X) log(1 + X(:,1).^2) .* sin(5 * X(:,2));

dim = 2;
max_order = 10;
a = [1.30748, 1.30748]; % Lower bounds for each dimension
b = [2.69253, 2.69253]; % Upper bounds for each dimension
nsamp_pce = 40;
nsamp_mcs = 1e4;

% Create an instance of the PolynomialChaosExpansion class
pce = PolynomialChaosExpansion(func, dim, max_order, a, b);

% Run the PCE process
pce.runPCE(nsamp_pce, nsamp_mcs);

% Optionally, save the results or handle them as needed