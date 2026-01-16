function plot_SE_vs_N_T_WOA(para, H, user_r, user_theta, N_T_range)
%Plot spectral efficiency versus number of TTDs for WOA methods
%  plot_SE_vs_N_T_WOA(para, H, user_r, user_theta, N_T_range)
%Inputs:
%   para: structure of the initial parameters
%   H: channel for all users
%   user_r: distance for all users
%   user_theta: angle for all users
%   N_T_range: array of number of TTDs values (e.g., [4, 8, 12, 16])
%Date: 17/01/2026

if nargin < 1
    para = para_init();
end
if nargin < 2
    user_r = rand(para.K, 1) * 10 + 5;
end
if nargin < 3
    user_theta = sort(rand(para.K, 1) * pi);
end
if nargin < 4
    H = generate_channel(para, user_r, user_theta);
end
if nargin < 5
    N_T_range = [8, 16, 32, 64]; % Default N_T range, divisors of 512
end

SE_WOA_PNF = zeros(length(N_T_range), 1);
SE_WOA_robust = zeros(length(N_T_range), 1);
SE_WOA_fully_digital = zeros(length(N_T_range), 1);
SE_orig_PNF = zeros(length(N_T_range), 1);
SE_orig_robust = zeros(length(N_T_range), 1);
SE_orig_fully_digital = zeros(length(N_T_range), 1);

for idx = 1:length(N_T_range)
    N_T = N_T_range(idx);
    para_temp = para;
    para_temp.N_T = N_T;
    
    % Ensure N_RF <= N_T or adjust if needed, but assume N_RF is fixed
    
    % Run original HTS PNF
    SE_orig_PNF(idx) = algorithm_HTS_PNF(para_temp, H, user_r, user_theta);
    
    % Run original HTS robust
    SE_orig_robust(idx) = algorithm_HTS_robust(para_temp, H, user_r, user_theta);
    
    % Run original fully digital
    SE_orig_fully_digital(idx) = algorithm_fully_digital(para_temp, H);
    
    % Run WOA HTS PNF
    SE_WOA_PNF(idx) = algorithm_HTS_PNF_WOA(para_temp, H, user_r, user_theta);
    
    % Run WOA HTS robust
    SE_WOA_robust(idx) = algorithm_HTS_robust_WOA(para_temp, H, user_r, user_theta);
    
    % Run WOA fully digital
    SE_WOA_fully_digital(idx) = algorithm_fully_digital_WOA(para_temp, H);
end

% Plot
figure;
plot(N_T_range, SE_orig_PNF, '-b', 'LineWidth', 1.5);
hold on;
plot(N_T_range, SE_orig_robust, '-.r', 'LineWidth', 1.5);
plot(N_T_range, SE_orig_fully_digital, ':k', 'LineWidth', 1.5);
plot(N_T_range, SE_WOA_PNF, '--g', 'LineWidth', 1.5);
plot(N_T_range, SE_WOA_robust, ':m', 'LineWidth', 1.5);
plot(N_T_range, SE_WOA_fully_digital, '-c', 'LineWidth', 1.5);
xlabel('Number of TTDs', 'Interpreter', 'Latex');
ylabel('Spectral Efficiency (bit/s/Hz)', 'Interpreter', 'Latex');
legend('HTS PNF (Original)', 'HTS Robust (Original)', 'Fully Digital (Original)', 'HTS PNF (WOA)', 'HTS Robust (WOA)', 'Fully Digital (WOA)', 'Interpreter', 'Latex');
title('Spectral Efficiency vs Number of TTDs', 'Interpreter', 'Latex');
grid on;
box on;
hold off;