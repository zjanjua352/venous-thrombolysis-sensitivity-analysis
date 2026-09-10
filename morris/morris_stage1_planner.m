% =========================================================================
% Script Name: morris_stage1_planner.m
% Description: Morris Global Sensitivity Analysis - Stage 1 (Trajectory Planner)
%              Generates randomized One-Factor-At-A-Time (OFAT) sampling paths,
%              applies logarithmic scaling, and exports the static execution matrix.
%
% Reference / Attribution:
% Developed, refactored, and parallelized in collaboration with 
% Google Gemini (Large Language Model, Google LLC) for Master's Thesis 
% Supplementary Materials.
% =========================================================================

clc; clear; close all;

% 1. Configuration & Parameter Selection
% Define parameter names and baseline physiological/physical values
params = { ...
    'Dcoeff', ...          % Protein Diffusivity [cm^2/s]
    'k10', ...             % Drug Elimination Rate Constant [1/s]
    'ka_drug', ...         % Drug Adsorption Constant [1/s]
    'kd_drug', ...         % Drug Desorption Constant [1/s]
    'epsilon_0', ...       % Initial Clot Porosity [-]
    'RG', ...              % Clot Retraction Ratio [-]
    'L_clot', ...          % Clot Length [cm]
    'H', ...               % Clot Hematocrit [-]
    'kPAI', ...            % Drug Inhibition Constant [1/(uM*s)]
    'drugDose', ...        % Bolus/Infusion Dose [mg]
    'R_f0', ...            % Initial Fibrin Fiber Radius [cm]
    'KM_PLG', ...          % Plasminogen Activation Constant [uM]
    't_infusion', ...      % Infusion Duration [s]
    'dpdx_clot', ...       % Clot Pressure Drop [mmHg/cm]
    'D_vein' ...           % Vessel Diameter [cm]
};

baselines = [ ...
    2.5e-7, ...            % Dcoeff
    0.0019, ...            % k10
    0.05, ...              % ka_drug
    0.005, ...             % kd_drug
    0.99, ...              % epsilon_0
    1.0, ...               % RG
    1.0, ...               % L_clot
    0.0, ...               % H
    0.001, ...             % kPAI
    100, ...               % drugDose
    7.0e-6, ...            % R_f0
    0.15, ...              % KM_PLG
    120*60, ...            % t_infusion
    1.0, ...               % dpdx_clot
    0.5 ...                % D_vein
];

% 2. Sampling Space Bounds & Logarithmic Mapping Setup
p = numel(params);         % Number of active parameters
r = 10;                    % Number of independent sampling trajectories
delta = 0.5;               % Normalized grid step size
factor = 1.5;              % Variation factor around baseline (+/-)

% Define lower and upper parameter limits
lb = baselines ./ factor;
ub = baselines .* factor;

% Apply physical boundary clipping
eps_idx = find(strcmp(params, 'epsilon_0'));
if ~isempty(eps_idx)
    ub(eps_idx) = min(ub(eps_idx), 0.999); % Porosity upper bound ceiling
end

% Function mapping normalized coordinates [0, 1] to actual log-scale parameters
mapToActual = @(z, low, high) low .* (high ./ low).^z;

% 3. Morris Trajectory Generation (OFAT Sampling)
runs_per_traj = p + 1;
total_sims = r * runs_per_traj;
design_matrix = zeros(total_sims, p);
trajectory_tracking = zeros(r, runs_per_traj, p);

rng(42, 'twister'); % Set seed for deterministic reproducibility

fprintf('--- Generating %d Morris Trajectories (%d Total Runs) ---\n', r, total_sims);

row_counter = 1;
for iTraj = 1:r
    % Sample randomized base point in restricted hypercube
    x0 = rand(1, p) * (1 - delta);
    
    % Randomized permutation of parameter perturbation sequence
    order = randperm(p);
    
    % Build trajectory steps in normalized unit space
    X_norm = zeros(runs_per_traj, p);
    X_norm(1, :) = x0;
    
    for j = 1:p
        X_norm(j + 1, :) = X_norm(j, :);
        X_norm(j + 1, order(j)) = X_norm(j, order(j)) + delta;
    end
    
    % Map normalized trajectory coordinates to actual physical dimensions
    for k = 1:runs_per_traj
        actual_vals = mapToActual(X_norm(k, :), lb, ub);
        design_matrix(row_counter, :) = actual_vals;
        trajectory_tracking(iTraj, k, :) = actual_vals;
        row_counter = row_counter + 1;
    end
end

% 4. Save Master Design Configuration
save('morris_design.mat', 'params', 'baselines', 'lb', 'ub', ...
    'design_matrix', 'trajectory_tracking', 'r', 'p', 'delta', ...
    'runs_per_traj', 'total_sims');

fprintf('Successfully compiled and saved "morris_design.mat".\n');
