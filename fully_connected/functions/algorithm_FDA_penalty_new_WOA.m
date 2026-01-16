function [R_convergence, penalty_convergence, A, D] = algorithm_FDA_penalty_new_WOA(para, H, user_r, user_theta)
%The new version of the panalty-based fully-digital approximation (FDA) method using WOA
%Compared to the method in this paper, the new version uses Whale Optimization Algorithm for optimization.
%
%   [R_convergence, penalty_convergence, A, D] = algorithm_FDA_penalty_new_WOA(para, H, user_r, user_theta)
%Inputs:
%   para: structure of the initial parameters
%   H: channel for all users
%   user_r: distance for all users (for initialization)
%   user_theta: angle for all users (for initialization)
%Outputs:
%   R_convergence: achievable rates at each iteration
%   penalty_convergence: penalty value at each iteration
%   A: optimized analog beamforming matrix
%   D: optimized digital beamforming matrix
%Date: 17/01/2026

R_convergence = [];
penalty_convergence = [];

% Initialization using the fully-digital solution and the PNF-based HTS approach
W_initial = randn(para.N, para.K) + 1i * randn(para.N, para.K);
W_initial = W_initial / norm(W_initial, 'fro') * sqrt(para.Pt);
[~, W] = algorithm_fully_digital(para, H, W_initial);
[~, A, D, t] = algorithm_HTS_PNF(para, H, user_r, user_theta, W_initial);

% WOA parameters
SearchAgents_no = 30; % Number of search agents
Max_iter = 100; % Maximum number of iterations
dim = 2 * para.N * para.K * para.M; % Dimension: real and imag parts of W
lb = -10 * ones(1, dim); % Lower bound
ub = 10 * ones(1, dim); % Upper bound

% Objective function
fobj = @(x) objective_function(x, para, H, A, D);

% Run WOA
[Best_score, Best_pos, Convergence_curve] = WOA(SearchAgents_no, Max_iter, lb, ub, dim, fobj);

% Extract best W
W = reshape(Best_pos(1:dim/2) + 1i * Best_pos(dim/2+1:end), para.N, para.K, para.M);

% Update A and D based on best W
[A, t] = update_analog_beamformer(para, W, D, t);
[D] = update_digital_beamformer(para, W, A, D);

% Calculate final rate
W_hybrid = zeros(para.N, para.K, para.M);
for m = 1:para.M
    W_hybrid(:,:,m) = A(:,:,m)*D(:,:,m);
end
[R_sum] = rate_fully_digital(para, W_hybrid, H);
R_convergence = R_sum/(para.M+para.Lcp);
penalty_convergence = -Best_score; % Since objective is negative rate

end

%% Objective function for WOA
function obj = objective_function(x, para, H, A, D)
    dim = length(x)/2;
    W_real = x(1:dim);
    W_imag = x(dim+1:end);
    W = reshape(W_real + 1i * W_imag, para.N, para.K, para.M);
    
    % Calculate rate
    W_hybrid = zeros(para.N, para.K, para.M);
    for m = 1:para.M
        W_hybrid(:,:,m) = A(:,:,m)*D(:,:,m);
    end
    [R_sum] = rate_fully_digital(para, W_hybrid, H);
    obj = -R_sum; % Minimize negative rate
end

%% Update digital beamformer
function [D] = update_digital_beamformer(para, W, A, D)
    for m = 1:para.M
        D(:,:,m) = pinv(A(:,:,m))*W(:,:,m); % Equation (32)
    end
end

%% Update analog beamformer
function [A, t] = update_analog_beamformer(para, P, D, t)

    t_search = 0:para.t_max/1e3:para.t_max; % search space of TTDs' time delay

    PD = zeros(para.N, para.N_RF, para.M);
    for m = 1:para.M
        PD(:,:,m) = P(:,:,m)*pinv(D(:,:,m));
    end

    % update PS coefficients
    N_sub = para.N/para.N_T;
    A_PS = zeros(para.N, para.N_RF);
    for n = 1:para.N_RF
        for l = 1:para.N_T
            p_l = PD((l-1)*N_sub+1:l*N_sub, n, :);  
            p_l = squeeze(p_l);
            t_l = t(l,n);
            a_l = sum(p_l*diag(exp(1i*2*pi*para.fm_all*t_l)), 2);
            A_PS((l-1)*N_sub+1:l*N_sub, n) = a_l./abs(a_l);
        end
    end

    % update TTD coefficients
    t = zeros(para.N_T, para.N_RF);
    for n = 1:para.N_RF
        for l = 1:para.N_T
            p_l = PD((l-1)*N_sub+1:l*N_sub, n, :);  
            p_l = squeeze(p_l);
            a_l = A_PS((l-1)*N_sub+1:l*N_sub, n);
            psi_nq = p_l'*a_l;
            obj_value = real(psi_nq.'*exp(-1i*2*pi*para.fm_all'*t_search));
            [~,I] = max(obj_value); % one-dimensional search
            t(l,n) = t_search(I);
        end
    end

    % calculate overall analog beamformer
    A = zeros(para.N, para.N_RF, para.M);
    for m = 1:para.M
        A(:,:,m) = analog_bamformer(para, A_PS, t, para.fm_all(m));
    end
end

%% Calculate the overall TTD-based analog beamformer
function [Am] = analog_bamformer(para, A, t, fm)
    e = ones(para.N/para.N_T,1);
    Tm = exp(-1i*2*pi*fm*t);
    Am = A .* kron(Tm, e);
end