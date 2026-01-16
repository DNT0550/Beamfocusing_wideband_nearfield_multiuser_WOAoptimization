clc
clear all
close all

addpath("functions/");
para = para_init();

% Parameters for plotting
user_r = rand(para.K, 1) * 10 + 5; % user distances 5 ~ 15 m
user_theta = sort(rand(para.K, 1) * pi); % user directions 0 ~ 180 degree

%% Generate channel matrix
[H] = generate_channel(para, user_r, user_theta);

%% Plot Spectral Efficiency vs Transmit Power
figure;
plot_SE_vs_Pt_WOA(para, H, user_r, user_theta);
title('Spectral Efficiency vs Transmit Power for WOA Methods');

%% Plot Spectral Efficiency vs Number of TTDs
figure;
plot_SE_vs_N_T_WOA(para, H, user_r, user_theta);
title('Spectral Efficiency vs Number of TTDs for WOA Methods');

%% Plot Energy Efficiency vs Number of TTDs
figure;
plot_EE_vs_N_T_WOA(para, H, user_r, user_theta);
title('Energy Efficiency vs Number of TTDs for WOA Methods');

%% Plot Spectral Efficiency vs Bandwidth
figure;
plot_SE_vs_B_WOA(para, H, user_r, user_theta);
title('Spectral Efficiency vs Bandwidth for WOA Methods');

%% Plot Spectral Efficiency vs Maximum Time Delay
figure;
plot_SE_vs_t_max_WOA(para, H, user_r, user_theta);
title('Spectral Efficiency vs Maximum Time Delay for WOA Methods');

%% Plot Array Gain
figure;
plot_array_gain_WOA(para, 45*pi/180, 10); % Using default single user
title('Array Gain for WOA Methods');