% =========================================================================
% Script Name: sobol_stage1_planner.m
% Description: Sobol Variance-Based Sensitivity Analysis - Stage 1 (Planner)
%              Generates the Saltelli quasi-random cross-sampling design matrix
%              for the reduced subset of highly sensitive parameters.
%
% Reference / Attribution:
% Developed, refactored, and parallelized in collaboration with 
% Google Gemini (Large Language Model, Google LLC) for Master's Thesis 
% Supplementary Materials.
% =========================================================================

clc; clear; close all;

% 1. Configuration & Parameter Selection (Reduced Space)
% Based on Morris screening, only highly sensitive interacting parameters are included
params = { ...
    'Dcoeff', ...          % Protein Diffusivity [cm^2/s]
    'ka_drug', ...         % Drug Adsorption Constant [1/s]
    'kd_drug', ...         % Drug Desorption Constant [1/s]
};

baselines = [ ...
    2.5e-7, ...            % Dcoeff
    0.05, ...              % ka_drug
    0.005, ...             % kd_drug
];

% 2. Sampling Space Bounds & Logarithmic Mapping Setup
p = numel(params);         % Number of active parameters
N = 500;                   % Base sample size (Generates N*(p+2) total runs)
factor = 1.5;              % Variation factor around baseline (+/-)

% Define lower and upper parameter limits
lb = baselines ./ factor;
ub = baselines .* factor;

% Function mapping normalized coordinates [0, 1] to actual log-scale parameters
mapToActual = @(z, low, high) low .* (high ./ low).^z;

% 3. Saltelli / Sobol Design Matrix Generation
total_sims = N * (p + 2);
design_matrix = zeros(total_sims, p);

rng(42, 'twister'); % Set seed for deterministic reproducibility

fprintf('--- Generating Sobol/Saltelli Design Matrix (%d Total Runs) ---\n', total_sims);

% Generate unified random base for A and B to maintain low discrepancy
rand_base = rand(N, 2*p);
A_norm = rand_base(:, 1:p);
B_norm = rand_base(:, p+1:2*p);

% Map normalized coordinates to physical dimensional values
A = mapToActual(A_norm, lb, ub);
B = mapToActual(B_norm, lb, ub);

% 4. Matrix Assembly (A, B, and Cross-Matrices AB_i)
design_matrix(1:N, :) = A;
design_matrix(N+1:2*N, :) = B;

row_idx = 2*N + 1;
for i = 1:p
    AB_i = A;
    AB_i(:, i) = B(:, i); % Swap column i from matrix B into matrix A
    design_matrix(row_idx : row_idx+N-1, :) = AB_i;
    row_idx = row_idx + N;
end

% 5. Save Master Design Configuration
save('sobol_design.mat', 'params', 'baselines', 'lb', 'ub', ...
    'design_matrix', 'N', 'p', 'total_sims');

fprintf('Successfully compiled and saved "sobol_design.mat".\n');
