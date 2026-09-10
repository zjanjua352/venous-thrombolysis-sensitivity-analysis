% =========================================================================
% Script Name: morris_stage2_runner.m
% Description: Morris Global Sensitivity Analysis - Stage 2 (Parallel Runner)
%              Coordinates batch profile pre-generation, thread-safe parfor
%              execution, file-based checkpointing, and dynamic progress monitoring.
%
% Reference / Attribution:
% Developed, refactored, and parallelized in collaboration with 
% Google Gemini (Large Language Model, Google LLC) for Master's Thesis 
% Supplementary Materials.
% =========================================================================

clc; clear; close all;

% 1. Load Design Configuration
if ~exist('morris_design.mat', 'file')
    error('Blueprint missing: Please execute "morris_stage1_planner.m" first.');
end
load('morris_design.mat');

% 2. Checkpoint Infrastructure & Output Directory
results_dir = 'morris_results/';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% 3. Phase 1: Sequential Boundary Profile Pre-generation
% Pre-computes 0D systemic plasma profiles to prevent parallel I/O race conditions
fprintf('--- Phase 1: Pre-generating Systemic Boundary Profiles ---\n');
for batchRun = 1:total_sims
    profile_file = sprintf('%sCsolve_%04d.mat', results_dir, batchRun);
    if ~exist(profile_file, 'file')
        this_run_params = design_matrix(batchRun, :);
        
        % Inject parameter perturbations into compartmental model
        for ii = 1:numel(params)
            name = params{ii};
            val = this_run_params(ii);
            eval([name, ' = val;']);
        end
        
        % Execute compartmental solver
        run('compartmentalModel_MorrisSA.m');
    end
end
fprintf('All boundary condition files successfully generated and cached.\n\n');

% 4. Phase 2: Parallel Pool Configuration & Progress Bar Setup
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local', 4); % Thread allocation balancing RAM footprint
end

fprintf('Initializing Parallel Progress Monitor (%d total evaluations)...\n', total_sims);
q = parallel.pool.DataQueue;
afterEach(q, @(x) update_progress_bar(total_sims));

% 5. Parallel Execution Loop
fprintf('--- Phase 2: Executing Parallel Morris Simulations ---\n');
tic;
parfor iRun = 1:total_sims
    output_filename = sprintf('%srun_%04d.mat', results_dir, iRun);
    
    % Checkpoint verification: Skip completed simulations
    if exist(output_filename, 'file')
        send(q, iRun);
        continue;
    end
    
    this_run_params = design_matrix(iRun, :);
    
    % Execute isolated simulation wrapper
    metric = execute_worker_simulation(this_run_params, params, iRun);
    
    % Save evaluation output via parfor-safe function
    save_parfor(output_filename, metric, this_run_params);
    send(q, iRun);
end
fprintf('\n');
toc;
fprintf('Morris execution queue complete. Results stored in "%s".\n', results_dir);

% ================== HELPER FUNCTIONS ==================

function out_metric = execute_worker_simulation(this_run_params, params, iRun)
    % Dynamically inject sampled parameter values into worker workspace
    for ii = 1:numel(params)
        name = params{ii};
        val = this_run_params(ii);
        eval([name, ' = val;']);
    end
    
    % Execute 1D transport-reaction solver
    run('integratedModel_new_main5venous_morrisSA.m');
    
    % Extract Extent of Lysis metric (dimensionally consistent fraction in [0, 1])
    if exist('extent_of_progress', 'var')
        out_metric = extent_of_progress;
    else
        out_metric = NaN; % Trap numerical/solver stagnation
    end
end

function save_parfor(filename, metric, params_val)
    % Thread-safe disk writer bypassing parfor workspace restrictions
    save(filename, 'metric', 'params_val');
end

function update_progress_bar(total_sims)
    % Stateful Command Window progress tracking
    persistent current_count;
    if isempty(current_count)
        current_count = 0;
    end
    current_count = current_count + 1;
    pct = (current_count / total_sims) * 100;
    fprintf('\rSimulation Progress: [%-30s] %5.1f%% (%d/%d)', ...
        repmat('=', 1, floor(pct / 3.33)), pct, current_count, total_sims);
    if current_count == total_sims
        current_count = 0; % Reset counter on completion
    end
end
