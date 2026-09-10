% =========================================================================
% Script Name: plot_morris_ranked_comparison.m
% Description: Generates publication-ready horizontal bar charts for Morris
%              sensitivity indices (mu* and sigma), formatted with group-based
%              coloring, automated LaTeX parameter mapping, and custom margins.
%
% Reference / Attribution:
% Developed, refactored, and parallelized in collaboration with 
% Google Gemini (Large Language Model, Google LLC) for Master's Thesis 
% Supplementary Materials.
% =========================================================================

clc; clear; close all;

% 1. Load Morris Results Datasets
data_files = { ...
    'morris_analysis/morris_kinetic_tnk.mat', ...
    'morris_analysis/morris_transport_.mat', ...
    'morris_analysis/morris_dosage_.mat', ...
    'morris_analysis/morris_permb_.mat' ...
};
labels = {'Kinetic Parameters'; 'Transport/Flow Parameters'; 'Dosage Parameters'; 'Clot Permeability Parameters'};

nDatasets = numel(data_files);
datasets = cell(nDatasets, 1);
for i = 1:nDatasets
    datasets{i} = load(data_files{i});
end

% 2. Build Master Parameter List & Map Indices
paramSets = cell(nDatasets, 1);
for i = 1:nDatasets
    MR = datasets{i}.MorrisResults;
    paramSets{i} = MR.paramNames(:);
end
masterParams = unique(cat(1, paramSets{:}));
nparams = numel(masterParams);

mu_star_mat = nan(nDatasets, nparams);
sigma_mat   = nan(nDatasets, nparams);

for i = 1:nDatasets
    MR = datasets{i}.MorrisResults;
    [~, loc] = ismember(MR.paramNames(:), masterParams);
    mu_star_mat(i, loc) = MR.mu_star;
    sigma_mat(i, loc)   = MR.sigma;
end

% 3. Flatten, Group-Color, and Rank
[flat_mu_star, group_idx] = max(mu_star_mat, [], 1, 'omitnan');
[flat_sigma, ~]           = max(sigma_mat, [], 1, 'omitnan');

[sorted_mu_star, sortIdx] = sort(flat_mu_star, 'ascend');
sortedParams              = masterParams(sortIdx);
sorted_sigma              = flat_sigma(sortIdx);
sorted_group              = group_idx(sortIdx);

% 4. Automated LaTeX Parameter Dictionary Mapping
keys = {'Dcoeff', 'k10', 'ka_drug', 'kd_drug', 'epsilon_0', 'RG', ...
        'L_clot', 'H', 'kPAI', 'drugDose', 'R_f0', 'KM_PLG', ...
        't_infusion', 'dpdx_clot', 'D_vein'};

values = { ...
    '\textit{Protein Diffusivity}, $D_{coeff}$', ...
    '\textit{Drug Elimination Rate}, $k_{10}$', ...
    '\textit{Drug Adsorption Constant}, $k_{a,drug}$', ...
    '\textit{Drug Desorption Constant}, $k_{d,drug}$', ...
    '\textit{Initial Clot Porosity}, $\epsilon_0$', ...
    '\textit{Clot Retraction}, $R_G$', ...
    '\textit{Clot Length}, $L_{clot}$', ...
    '\textit{Clot Hematocrit}, $H$', ...
    '\textit{Drug Inhibition Constant}, $k_{PAI}$', ...
    '\textit{Bolus/Infusion Dose}, $Dose$', ...
    '\textit{Fibrin Radius}, $R_{f0}$', ...
    '\textit{Plasminogen Activation}, $K_{M,PLG}$', ...
    '\textit{Infusion Duration}, $t_{inf}$', ...
    '\textit{Clot Pressure Drop}, $\Delta P / \Delta x$', ...
    '\textit{Vessel Diameter}, $D_{vein}$' ...
};

paramMap = containers.Map(keys, values);
displayParams = cell(nparams, 1);
for p = 1:nparams
    rawName = sortedParams{p};
    if isKey(paramMap, rawName)
        displayParams{p} = paramMap(rawName);
    else
        displayParams{p} = strrep(rawName, '_', '\_');
    end
end

% 5. Plotting Figure 1: mu* (Overall Importance)
figure('Name', 'Morris mu* Ranking', 'NumberTitle', 'off', 'Position', [100, 100, 900, 750]);
dsColors = lines(nDatasets);
yPos = 1:nparams;
barWidth = 0.75;
hold on;

for k = 1:nparams
    yline(k, ':', 'Color', [0.85 0.85 0.85], 'HandleVisibility', 'off');
end

for p = 1:nparams
    val = sorted_mu_star(p);
    g   = sorted_group(p);
    if ~isnan(val)
        rectangle('Position', [min(val, eps), yPos(p) - barWidth/2, abs(val - min(val, eps)), barWidth], ...
            'FaceColor', dsColors(g, :), 'EdgeColor', 'k', 'LineWidth', 0.5);
    end
end

set(gca, 'XScale', 'log', 'YTick', yPos, 'YTickLabel', displayParams, ...
    'TickLabelInterpreter', 'latex', 'FontSize', 12);
xlabel('$\mu^*$ (Overall Importance)', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
title('Global Sensitivity: $\mu^*$ Ranking', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
xlim([eps, max(sorted_mu_star) * 2]);
ylim([0.5, nparams + 0.5]);
grid on; box on;
set(gca, 'Position', [0.35, 0.15, 0.60, 0.75]); % Left margin allocation for LaTeX text

hDS = gobjects(nDatasets, 1);
for d = 1:nDatasets
    hDS(d) = patch(nan, nan, dsColors(d, :), 'EdgeColor', 'k');
end
lgd1 = legend(hDS, labels, 'Orientation', 'horizontal', 'FontSize', 12, 'Interpreter', 'tex');
lgd1.Position = [0.35, 0.02, 0.60, 0.05];
hold off;

% 6. Plotting Figure 2: sigma (Non-Linearity / Interactions)
figure('Name', 'Morris Sigma Ranking', 'NumberTitle', 'off', 'Position', [1050, 100, 900, 750]);
hold on;

for k = 1:nparams
    yline(k, ':', 'Color', [0.85 0.85 0.85], 'HandleVisibility', 'off');
end

for p = 1:nparams
    val = sorted_sigma(p);
    g   = sorted_group(p);
    if ~isnan(val)
        rectangle('Position', [min(val, eps), yPos(p) - barWidth/2, abs(val - min(val, eps)), barWidth], ...
            'FaceColor', dsColors(g, :), 'EdgeColor', 'k', 'LineWidth', 0.5);
    end
end

set(gca, 'XScale', 'log', 'YTick', yPos, 'YTickLabel', displayParams, ...
    'TickLabelInterpreter', 'latex', 'FontSize', 12);
xlabel('$\sigma$ (Non-Linearity / Interactions)', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
title('Global Sensitivity: $\sigma$ Ranking', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
xlim([eps, max(sorted_sigma) * 2]);
ylim([0.5, nparams + 0.5]);
grid on; box on;
set(gca, 'Position', [0.35, 0.15, 0.60, 0.75]);

for d = 1:nDatasets
    hDS(d) = patch(nan, nan, dsColors(d, :), 'EdgeColor', 'k');
end
lgd2 = legend(hDS, labels, 'Orientation', 'horizontal', 'FontSize', 12, 'Interpreter', 'tex');
lgd2.Position = [0.35, 0.02, 0.60, 0.05];
hold off;
