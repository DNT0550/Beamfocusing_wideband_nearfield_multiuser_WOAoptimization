function plot_SE_vs_t_max_WOA(para, H, user_r, user_theta, t_max_range)
%Plot spectral efficiency versus maximum time delay for WOA methods
%  plot_SE_vs_t_max_WOA(para, H, user_r, user_theta, t_max_range)
%Inputs:
%   para: structure of the initial parameters
%   H: channel for all users
%   user_r: distance for all users
%   user_theta: angle for all users
%   t_max_range: array of maximum time delay values (e.g., [1e-9, 2e-9, 3e-9])
%Date: 17/01/2026

if nargin < 5
    t_max_range = [1e-9, 2e-9, 3e-9, 4e-9, 5e-9]; % Default t_max range in seconds
end

SE_WOA_PNF = zeros(length(t_max_range), 1);
SE_WOA_robust = zeros(length(t_max_range), 1);
SE_orig_PNF = zeros(length(t_max_range), 1);
SE_orig_robust = zeros(length(t_max_range), 1);

for idx = 1:length(t_max_range)
    t_max = t_max_range(idx);
    para_temp = para;
    para_temp.t_max = t_max;
    
    % Run original HTS PNF
    SE_orig_PNF(idx) = algorithm_HTS_PNF(para_temp, H, user_r, user_theta);
    
    % Run original HTS robust
    SE_orig_robust(idx) = algorithm_HTS_robust(para_temp, H, user_r, user_theta);
    
    % Run WOA HTS PNF
    SE_WOA_PNF(idx) = algorithm_HTS_PNF_WOA(para_temp, H, user_r, user_theta);
    
    % Run WOA HTS robust
    SE_WOA_robust(idx) = algorithm_HTS_robust_WOA(para_temp, H, user_r, user_theta);
end

% Plot
figure;
plot(t_max_range * 1e9, SE_orig_PNF, '-b', 'LineWidth', 1.5);
hold on;
plot(t_max_range * 1e9, SE_orig_robust, '-.r', 'LineWidth', 1.5);
plot(t_max_range * 1e9, SE_WOA_PNF, '--g', 'LineWidth', 1.5);
plot(t_max_range * 1e9, SE_WOA_robust, ':m', 'LineWidth', 1.5);
xlabel('Maximum Time Delay (ns)', 'Interpreter', 'Latex');
ylabel('Spectral Efficiency (bit/s/Hz)', 'Interpreter', 'Latex');
legend('HTS PNF (Original)', 'HTS Robust (Original)', 'HTS PNF (WOA)', 'HTS Robust (WOA)', 'Interpreter', 'Latex');
title('Spectral Efficiency vs Maximum Time Delay', 'Interpreter', 'Latex');
grid on;
box on;
hold off;