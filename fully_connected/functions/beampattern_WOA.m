function [P_WOA_PNF, P_WOA_robust] = beampattern_WOA(para, theta, r, A_PNF, t_PNF, A_robust, t_robust)
%Calculate the beam pattern achieved by WOA-based methods
%  [P_WOA_PNF, P_WOA_robust] = beampattern_WOA(para, theta, r, A_PNF, t_PNF, A_robust, t_robust)
%Inputs:
%   para: structure of the initial parameters
%   theta: angle
%   r: distance
%   A_PNF: analog beamformer from HTS_PNF_WOA
%   t_PNF: time delays from HTS_PNF_WOA
%   A_robust: analog beamformer from HTS_robust_WOA
%   t_robust: time delays from HTS_robust_WOA
%Outputs:
%   P_WOA_PNF: beampattern achieved by the WOA PNF method
%   P_WOA_robust: beampattern achieved by the WOA robust method
%Date: 17/01/2026

c = 3e8; % speed of light
N_sub = para.N/para.N_T; % the number of antennas connected to each TTD
e = ones(N_sub, 1);

%% WOA PNF method
if nargin >= 5 && ~isempty(A_PNF) && ~isempty(t_PNF)
    P_WOA_PNF = zeros(para.M, 1);
    for m = 1:para.M
        fm = para.fm_all(m);
        a_m = A_PNF(:,1,m); % Assuming single RF chain for simplicity, or adjust for multi-user
        bm = array_response_vector(r, theta, para.N, para.d, fm);
        P_WOA_PNF(m) = abs(bm.' * a_m);
    end
else
    P_WOA_PNF = [];
end

%% WOA robust method
if nargin >= 7 && ~isempty(A_robust) && ~isempty(t_robust)
    P_WOA_robust = zeros(para.M, 1);
    for m = 1:para.M
        fm = para.fm_all(m);
        a_m = A_robust(:,1,m); % Assuming single RF chain
        bm = array_response_vector(r, theta, para.N, para.d, fm);
        P_WOA_robust(m) = abs(bm.' * a_m);
    end
else
    P_WOA_robust = [];
end

end