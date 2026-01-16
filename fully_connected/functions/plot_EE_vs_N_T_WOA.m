function plot_EE_vs_N_T_WOA(para, H, user_r, user_theta, N_T_range)
%Plot energy efficiency versus number of TTDs for WOA methods
%  plot_EE_vs_N_T_WOA(para, H, user_r, user_theta, N_T_range)
%Inputs:
%   para: structure of the initial parameters
%   H: channel for all users
%   user_r: distance for all users
%   user_theta: angle for all users
%   N_T_range: array of number of TTDs values (e.g., [4, 8, 12, 16])
%Date: 17/01/2026

if nargin < 5
    N_T_range = [4, 8, 12, 16, 20]; % Default N_T range
end

EE_WOA_PNF = zeros(length(N_T_range), 1);
EE_WOA_robust = zeros(length(N_T_range), 1);
EE_orig_PNF = zeros(length(N_T_range), 1);
EE_orig_robust = zeros(length(N_T_range), 1);

for idx = 1:length(N_T_range)
    N_T = N_T_range(idx);
    para_temp = para;
    para_temp.N_T = N_T;
    
    % Run original HTS PNF
    SE_orig_PNF = algorithm_HTS_PNF(para_temp, H, user_r, user_theta);
    EE_orig_PNF(idx) = SE_orig_PNF / para_temp.Pt;
    
    % Run original HTS robust
    SE_orig_robust = algorithm_HTS_robust(para_temp, H, user_r, user_theta);
    EE_orig_robust(idx) = SE_orig_robust / para_temp.Pt;
    
    % Run WOA HTS PNF
    SE_WOA_PNF = algorithm_HTS_PNF_WOA(para_temp, H, user_r, user_theta);
    EE_WOA_PNF(idx) = SE_WOA_PNF / para_temp.Pt;
    
    % Run WOA HTS robust
    SE_WOA_robust = algorithm_HTS_robust_WOA(para_temp, H, user_r, user_theta);
    EE_WOA_robust(idx) = SE_WOA_robust / para_temp.Pt;
end

% Plot
figure;
plot(N_T_range, EE_orig_PNF, '-b', 'LineWidth', 1.5);
hold on;
plot(N_T_range, EE_orig_robust, '-.r', 'LineWidth', 1.5);
plot(N_T_range, EE_WOA_PNF, '--g', 'LineWidth', 1.5);
plot(N_T_range, EE_WOA_robust, ':m', 'LineWidth', 1.5);
xlabel('Number of TTDs', 'Interpreter', 'Latex');
ylabel('Energy Efficiency (bit/s/Hz/W)', 'Interpreter', 'Latex');
legend('HTS PNF (Original)', 'HTS Robust (Original)', 'HTS PNF (WOA)', 'HTS Robust (WOA)', 'Interpreter', 'Latex');
title('Energy Efficiency vs Number of TTDs', 'Interpreter', 'Latex');
grid on;
box on;
hold off;