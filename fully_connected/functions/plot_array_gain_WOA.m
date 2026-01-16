function plot_array_gain_WOA(para, theta, r, P_WOA_PNF, P_WOA_robust, B_values)
%Plot array gain for WOA methods compared to original methods
%  plot_array_gain_WOA(para, theta, r, P_WOA_PNF, P_WOA_robust, B_values)
%Inputs:
%   para: structure of the initial parameters
%   theta: user angle
%   r: user distance
%   P_WOA_PNF: beampattern from WOA PNF (optional)
%   P_WOA_robust: beampattern from WOA robust (optional)
%   B_values: array of bandwidths to plot (e.g., [1e10, 2e10, 3e10])
%Date: 17/01/2026

if nargin < 6
    B_values = [1e10, 2e10, 3e10]; % Default bandwidths
end

figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');

for idx = 1:length(B_values)
    B = B_values(idx);
    m = 1:para.M;
    para.fm_all = para.fc + B*(2*m-1-para.M) / (2*para.M); % subcarrier frequencies

    % Calculate original beampatterns
    [P_prop, P_prop_robust, P_conv_CF, P_con_MCCM, P_conv_MCM] = beampattern(para, theta, r);

    subplot(length(B_values), 1, idx); hold on; box on;
    plot(para.fm_all/1e9, 10*log10(P_prop/para.N), '-', 'LineWidth', 1.5, 'Color', 'b');
    plot(para.fm_all/1e9, 10*log10(P_prop_robust/para.N), '-.', 'LineWidth', 1.5, 'Color', 'r');
    if ~isempty(P_WOA_PNF)
        plot(para.fm_all/1e9, 10*log10(P_WOA_PNF/para.N), '--', 'LineWidth', 1.5, 'Color', 'g');
    end
    if ~isempty(P_WOA_robust)
        plot(para.fm_all/1e9, 10*log10(P_WOA_robust/para.N), ':', 'LineWidth', 1.5, 'Color', 'm');
    end
    plot(para.fm_all/1e9, 10*log10(P_con_MCCM/para.N), '--', 'LineWidth', 1.5, 'Color', 'k');
    plot(para.fm_all/1e9, 10*log10(P_conv_MCM/para.N), '--', 'LineWidth', 1.5, 'Color', 'c');
    plot(para.fm_all/1e9, 10*log10(P_conv_CF/para.N), ':', 'LineWidth', 1.5, 'Color', 'y');

    xlabel('Frequency (GHz)', 'Interpreter', 'Latex');
    ylabel('Array Gain (dB)', 'Interpreter', 'Latex');
    title(sprintf('$B = %d$ GHz', B/1e9), 'Interpreter', 'Latex');
end

% Legend for the last subplot
legend("TTD-BF, Proposed method","TTD-BF, Proposed robust method",...
    "TTD-BF, WOA PNF", "TTD-BF, WOA robust",...
    "Conventional BF, MCCM", "Conventional BF, MCM", "Conventional BF, CF", 'Interpreter', 'Latex');

end