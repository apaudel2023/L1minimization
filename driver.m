clc; clear all; close all; 

% Add the path to the class file if it's not in the current directory
% addpath('path_to_class_file');

% ======================================
% Define the function and parameters for PCE
func = @(X) log(1 + X(:,1).^2) .* sin(5 * X(:,2));

dim = 2;

a = [1.30748, 1.30748]; % Lower bounds for each dimension
b = [2.69253, 2.69253]; % Upper bounds for each dimension
max_order = 10;
nsamp_pce = 40;
nsamp_mcs = 1e4;
result_file_name = 'PCE_l2';

% ======================================
% Create an instance of the PolynomialChaosExpansion class
pce = PolynomialChaosExpansion(func, dim, max_order, a, b, result_file_name);

% Run the PCE process
pce.runPCE(nsamp_pce, nsamp_mcs);

% ======================================
plotter = PlotResults(result_file_name, 'FontSize', 16, 'LineWidth', 3);

% Load results
results = plotter.loadResults();

% Plot PDF (Assuming 'results' has the necessary data for plotting)
plotter.pdf_pce_mcs(results);


% Plot KL Divergence Convergence
plotter.plotKLConvergence(results.PCE.kld);

% ======================================
% move result to a seperate directory
result_dir_name = 'PCE-RESULTS';
result_dir = moveResults(result_dir_name);
fprintf('Results are saved in the directory below:\n');
fprintf(result_dir);

