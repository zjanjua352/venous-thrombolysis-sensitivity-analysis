% =========================================================================
% Script Name: morris_stage3_analyzer.m
% Description: Morris Global Sensitivity Analysis - Stage 3 (Analyzer)
%              Aggregates simulation checkpoints, reconstructs OFAT paths,
%              applies physiological dead-zone filtering, and calculates
%              Elementary Effects (EE), mu, mu*, and sigma indices.
%
% Reference / Attribution:
% Developed, refactored, and parallelized in collaboration with 
% Google Gemini (Large Language Model, Google LLC) for Master's Thesis 
% Supplementary Materials.
% =========================================================================

clc; clear; close all;

% 1. Load Blueprint & Checkpoints
if ~exist('morris_design.mat', 'file')
    error('Master design file missing.');
end
load('morris_design.mat');

results_dir = 'morris_results/';
Y_all = nan(total_sims, 1);

for i = 1:total_sims
    filename = sprintf('%srun_%04d.mat', results_dir, i);
    if exist(filename, 'file')
        data = load(filename);
        % Physiological dead-zone threshold filter
        if ~isnan(data.metric) && data.metric >= 0.005
            Y_all(i) = data.metric;
        else
            Y_all(i) = NaN; % Stagnant lysis flagged for filtering
        end
    end
end

% 2. Reconstruct Trajectories & Calculate Elementary Effects
EE = nan(r, p);
valid_traj_count = 0;

for iTraj = 1:r
    idx_start = (iTraj - 1) * runs_per_traj + 1;
    idx_end   = iTraj * runs_per_traj;
    
    Y_traj = Y_all(idx_start:idx_end);
    
    % Drop trajectory if any point encountered a dead-zone
    if any(isnan(Y_traj))
        continue;
    end
    
    valid_traj_count = valid_traj_count + 1;
    X_actual = squeeze(trajectory_tracking(iTraj, :, :));
    
    % Calculate Elementary Effect (EE) for each step in the trajectory
    for step = 1:p
        step_idx1 = step;
        step_idx2 = step + 1;
        
        dX = X_actual(step_idx2, :) - X_actual(step_idx1, :);
        moved_param = find(abs(dX) > 1e-12); % Identify perturbed parameter
        
        if ~isempty(moved_param)
            delta_actual = dX(moved_param);
            dY = Y_traj(step_idx2) - Y_traj(step_idx1);
            EE(iTraj, moved_param) = dY / delta_actual; % Dimensional EE ratio
        end
    end
end

% 3. Compute Sensitivity Indices
mu      = mean(EE, 1, 'omitnan');       % Algebraic Directional Mean
mu_star = mean(abs(EE), 1, 'omitnan');  % Absolute Overall Importance Index
sigma   = std(EE, 0, 1, 'omitnan');     % Interaction / Non-linearity Index

% 4. Display Results in Formatted Table
fprintf('\n================== Morris Sensitivity Results ==================\n');
fprintf('Valid Trajectories Analyzed: %d / %d\n', valid_traj_count, r);
fprintf('%-15s | %-12s | %-12s | %-12s\n', 'Parameter', 'mu*', 'mu', 'sigma');
fprintf('----------------------------------------------------------------\n');
for i = 1:p
    fprintf('%-15s | %12.4e | %12.4e | %12.4e\n', ...
        params{i}, mu_star(i), mu(i), sigma(i));
end
fprintf('================================================================\n\n');

% 5. Export Compiled Results Struct
MorrisResults.paramNames = params;
MorrisResults.baselines  = baselines;
MorrisResults.lb         = lb;
MorrisResults.ub         = ub;
MorrisResults.EE         = EE;
MorrisResults.mu         = mu;
MorrisResults.mu_star    = mu_star;
MorrisResults.sigma      = sigma;
MorrisResults.valid_traj = valid_traj_count;

save('morris_results_final.mat', 'MorrisResults');
