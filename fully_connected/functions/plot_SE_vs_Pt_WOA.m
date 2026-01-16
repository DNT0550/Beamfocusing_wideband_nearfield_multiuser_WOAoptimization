function plot_SE_vs_Pt_WOA(para, H, user_r, user_theta, Pt_range)
%Plot spectral efficiency versus transmit power for WOA methods
%  plot_SE_vs_Pt_WOA(para, H, user_r, user_theta, Pt_range)
%Inputs:
%   para: structure of the initial parameters
%   H: channel for all users
%   user_r: distance for all users
%   user_theta: angle for all users
%   Pt_range: array of transmit power values (e.g., [10, 20, 30])
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
    Pt_range = [10, 20, 30, 40, 50]; % Default Pt range
end

SE_WOA_PNF = zeros(length(Pt_range), 1);
SE_WOA_robust = zeros(length(Pt_range), 1);
SE_WOA_fully_digital = zeros(length(Pt_range), 1);
SE_WOA_FDA_penalty = zeros(length(Pt_range), 1);
SE_orig_PNF = zeros(length(Pt_range), 1);
SE_orig_robust = zeros(length(Pt_range), 1);
SE_orig_fully_digital = zeros(length(Pt_range), 1);
SE_orig_FDA_penalty = zeros(length(Pt_range), 1);

for idx = 1:length(Pt_range)
    Pt = Pt_range(idx);
    para_temp = para;
    para_temp.Pt = Pt;
    
    % Run original HTS PNF
    SE_orig_PNF(idx) = algorithm_HTS_PNF(para_temp, H, user_r, user_theta);
    
    % Run original HTS robust
    SE_orig_robust(idx) = algorithm_HTS_robust(para_temp, H, user_r, user_theta);
    
    % Run original fully digital
    P_initial = randn(para_temp.N, para_temp.K) + 1i * randn(para_temp.N, para_temp.K);
    SE_orig_fully_digital(idx) = algorithm_fully_digital(para_temp, H, P_initial);
    
    % Run original FDA penalty
    SE_orig_FDA_penalty(idx) = algorithm_FDA_penalty_new(para_temp, H, user_r, user_theta);
    
    % Run WOA HTS PNF
    SE_WOA_PNF(idx) = algorithm_HTS_PNF_WOA(para_temp, H, user_r, user_theta);
    
    % Run WOA HTS robust
    SE_WOA_robust(idx) = algorithm_HTS_robust_WOA(para_temp, H, user_r, user_theta);
    
    % Run WOA fully digital
    SE_WOA_fully_digital(idx) = algorithm_fully_digital_WOA(para_temp, H, P_initial);
    
    % Run WOA FDA penalty
    SE_WOA_FDA_penalty(idx) = algorithm_FDA_penalty_new_WOA(para_temp, H, user_r, user_theta);
end

% Plot
figure;
plot(Pt_range, SE_orig_PNF, '-b', 'LineWidth', 1.5);
hold on;
plot(Pt_range, SE_orig_robust, '-.r', 'LineWidth', 1.5);
plot(Pt_range, SE_orig_fully_digital, ':k', 'LineWidth', 1.5);
plot(Pt_range, SE_orig_FDA_penalty, '--y', 'LineWidth', 1.5);
plot(Pt_range, SE_WOA_PNF, '--g', 'LineWidth', 1.5);
plot(Pt_range, SE_WOA_robust, ':m', 'LineWidth', 1.5);
plot(Pt_range, SE_WOA_fully_digital, '-c', 'LineWidth', 1.5);
plot(Pt_range, SE_WOA_FDA_penalty, '-.b', 'LineWidth', 1.5);
xlabel('Transmit Power (W)', 'Interpreter', 'Latex');
ylabel('Spectral Efficiency (bit/s/Hz)', 'Interpreter', 'Latex');
legend('HTS PNF (Original)', 'HTS Robust (Original)', 'Fully Digital (Original)', 'FDA Penalty (Original)', 'HTS PNF (WOA)', 'HTS Robust (WOA)', 'Fully Digital (WOA)', 'FDA Penalty (WOA)', 'Interpreter', 'Latex');
title('Spectral Efficiency vs Transmit Power', 'Interpreter', 'Latex');
grid on;
box on;
hold off;