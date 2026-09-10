% =========================================================================
% Script Name: sobol_stage3_analyzer.m
% Description: Sobol Variance-Based Sensitivity Analysis - Stage 3 (Analyzer)
%              Features auto-truncation for interrupted Saltelli blocks,
%              physiological dead-zone filtering, Jansen estimator variance
%              math for S_i and S_Ti, and interaction gap plotting.
%
% Reference / Attribution:
% Developed, refactored, and parallelized in collaboration with 
% Google Gemini (Large Language Model, Google LLC) for Master's Thesis 
% Supplementary Materials.
% =========================================================================

clc; clear all; close all;

% 1. Load the Blueprint & Initialize
if ~exist('sobol_design.mat', 'file')
    error('Design file not found. Please ensure Stage 1 was run.');
end
load('sobol_design.mat'); % Loads params, N (base samples), p (num params)
p = numel(params);
block_size = p + 2; % Saltelli block constraint

% 2. Auto-Truncation / File Counting (The Recovery Logic)
results_dir = 'sobol_results/';
files = dir([results_dir, 'run_*.mat']);
num_files_found = length(files);

% Calculate mathematically complete Saltelli blocks
complete_blocks = floor(num_files_found / block_size);
valid_runs = complete_blocks * block_size;

fprintf('Found %d completed runs in the directory.\n', num_files_found);
fprintf('Auto-truncating to %d complete Saltelli blocks (%d valid runs).\n', complete_blocks, valid_runs);

% 3. Load Data & Apply "Dead Zone" Filter
Y = NaN(valid_runs, 1);
for i = 1:valid_runs
    filename = sprintf('%srun_%04d.mat', results_dir, i);
    data = load(filename);
    
    % Dead-Zone Filter: Stagnant physiology tracking
    if data.metric < 0.005
        Y(i) = NaN;
    else
        Y(i) = data.metric;
    end
end

% 4. Reshape & Clean Matrix Blocks
% Reshape array into [complete_blocks x (p + 2)]
Y_blocks = reshape(Y, [complete_blocks, block_size]);

% Find blocks where ANY run resulted in NaN and drop them
broken_blocks = any(isnan(Y_blocks), 2);
Y_clean = Y_blocks(~broken_blocks, :);
N_clean = size(Y_clean, 1);

fprintf('Dropped %d broken blocks (physiological dead-zones).\n', sum(broken_blocks));
fprintf('Calculating variance using %d fully intact physiological blocks.\n', N_clean);

if N_clean == 0
    error('No fully intact blocks remained after dead-zone filtering.');
end

% 5. Decompose into A, B, and A_B Matrices
YA = Y_clean(:, 1);
YB = Y_clean(:, 2);
Y_AB = Y_clean(:, 3:end); % [N_clean x p] matrix of swapped rows

% Calculate Global Variance V(Y) across all independent baseline samples
Y_all = [YA; YB];
V_Y = var(Y_all);

% 6. Jansen Variance Math (Sobol Indices)
Si = zeros(p, 1);
STi = zeros(p, 1);

for i = 1:p
    % First-Order Index (Si): Main Effect
    V_i = V_Y - (1 / (2 * N_clean)) * sum((YB - Y_AB(:, i)).^2);
    Si(i) = max(0, V_i / V_Y); % max(0) corrects for tiny numerical drift
    
    % Total-Order Index (STi): Total Authority + Interactions
    V_Ti = (1 / (2 * N_clean)) * sum((YA - Y_AB(:, i)).^2);
    STi(i) = max(0, V_Ti / V_Y);
end

% 7. Console Output & Plotting
fprintf('\n--- Sobol Sensitivity Results (Jansen Estimator) ---\n');
for i = 1:p
    fprintf('%-15s: Si = %6.4f  |  STi = %6.4f  |  Interaction Gap = %6.4f\n', ...
        params{i}, Si(i), STi(i), STi(i) - Si(i));
end

% Generate Grouped Bar Chart
figure('Name','Sobol Sensitivity Indices','NumberTitle','off', 'Position', [150, 150, 900, 600]);

plot_data = [Si, STi];
b = bar(plot_data, 'grouped');
b(1).FaceColor = [0 0.4470 0.7410]; % Blue (Si)
b(2).FaceColor = [0.8500 0.3250 0.0980]; % Red (STi)

set(gca, 'XTickLabel', strrep(params, '_', '\_'), 'TickLabelInterpreter', 'tex', 'FontSize', 12);
ylabel('Percentage of Model Variance', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Sobol Variance Decomposition (Calculated from %d valid blocks)', N_clean), 'FontSize', 16, 'FontWeight', 'bold');

lgd = legend({'First-Order Main Effect (S_i)', 'Total-Order Authority (S_{Ti})'}, 'Location', 'northeast', 'FontSize', 12);
grid on; box on;
