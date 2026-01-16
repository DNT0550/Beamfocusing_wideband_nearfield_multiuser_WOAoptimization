clc
clear all
close all

addpath("functions/");
para = para_init();

theta = 45*pi/180; % user direction
r = 10; % user distance

para.N_T = 16; % number of TTDs
para.M = 256; % number of subcarriers

% Bandwidths to plot
B_values = [1e10, 2e10, 3e10]; % 10 GHz, 20 GHz, 30 GHz

% For simplicity, assume single user
user_r = r * ones(para.K, 1);
user_theta = theta * ones(para.K, 1);

% Initialize P for WOA
P_WOA_PNF_all = cell(length(B_values), 1);
P_WOA_robust_all = cell(length(B_values), 1);

for idx = 1:length(B_values)
    B = B_values(idx);
    m = 1:para.M;
    para.fm_all = para.fc + B*(2*m-1-para.M) / (2*para.M); % subcarrier frequencies
    
    % Generate channel
    H = generate_channel(para, user_r, user_theta);
    
    % Run WOA HTS PNF
    [R_PNF_WOA, A_PNF_WOA, D_PNF_WOA, t_PNF_WOA] = algorithm_HTS_PNF_WOA(para, H, user_r, user_theta);
    
    % Run WOA HTS robust
    [R_robust_WOA, A_robust_WOA, D_robust_WOA, t_robust_WOA] = algorithm_HTS_robust_WOA(para, H, user_r, user_theta);
    
    % Calculate beampattern for WOA
    [P_WOA_PNF, P_WOA_robust] = beampattern_WOA(para, theta, r, A_PNF_WOA, t_PNF_WOA, A_robust_WOA, t_robust_WOA);
    
    P_WOA_PNF_all{idx} = P_WOA_PNF;
    P_WOA_robust_all{idx} = P_WOA_robust;
end

% Plot for each bandwidth
for idx = 1:length(B_values)
    B = B_values(idx);
    m = 1:para.M;
    para.fm_all = para.fc + B*(2*m-1-para.M) / (2*para.M);
    
    P_WOA_PNF = P_WOA_PNF_all{idx};
    P_WOA_robust = P_WOA_robust_all{idx};
    
    % Call plot function
    plot_array_gain_WOA(para, theta, r, P_WOA_PNF, P_WOA_robust, B_values(idx));
end