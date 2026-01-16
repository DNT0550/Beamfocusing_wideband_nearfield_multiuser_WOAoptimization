function [R, P] = algorithm_fully_digital_WOA(para, h, P_initial)
%WMMSE algorithm for fully digital systems using WOA
%  [R, P] = algorithm_fully_digital_WOA(para, h, P_initial)
%Inputs:
%   para: structure of the initial parameters
%   h: channel for all users
%   P_initial: initial beamforming vectors
%Outputs:
%   R: achievable rates
%   P: optimal fully digital beamformers
%Date: 17/01/2026

P = zeros(para.N, para.K, para.M);
for m = 1:para.M
    % WOA parameters
    SearchAgents_no = 30;
    Max_iter = 100;
    dim = 2 * para.N * para.K; % Real and imag parts
    lb = -10 * ones(1, dim);
    ub = 10 * ones(1, dim);
    
    % Objective function
    fobj = @(x) objective_function(x, para, h(:,:,m));
    
    % Run WOA
    [Best_score, Best_pos, ~] = WOA(SearchAgents_no, Max_iter, lb, ub, dim, fobj);
    
    % Extract best P
    Pm_real = Best_pos(1:dim/2);
    Pm_imag = Best_pos(dim/2+1:end);
    P(:,:,m) = reshape(Pm_real + 1i * Pm_imag, para.N, para.K);
end
[R] = rate_fully_digital(para, P, h);
R = R/(para.M+para.Lcp);

end

%% Objective function
function obj = objective_function(x, para, h)
    dim = length(x)/2;
    P_real = x(1:dim);
    P_imag = x(dim+1:end);
    P = reshape(P_real + 1i * P_imag, para.N, para.K);
    
    [R_sum] = rate_single(para, P, h);
    obj = -R_sum; % Minimize negative rate
end

%% achievable rate
function [R_sum, R] = rate_single(para, P, h)

R = zeros(para.K, 1);
for k = 1:para.K
    hk = h(:,k);
    pk = P(:,k); 
    P_I = P; P_I(:,k) = [];
    Ik = norm(hk'*P_I)^2 + norm(P, 'fro')^2/para.Pt; 
    R(k) = log2( 1 + abs(hk'*pk)^2/Ik );
end
R_sum = sum(R);

end