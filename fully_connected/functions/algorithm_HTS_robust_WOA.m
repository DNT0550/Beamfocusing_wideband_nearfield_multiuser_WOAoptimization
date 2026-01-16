function [R, A, D, t] = algorithm_HTS_robust_WOA(para, H, user_r, user_theta, W_initial, theta_plot, r_plot, plot_flag)
%The robust heuristic two-stage approach using WOA
%  [R, A, D, t] = algorithm_HTS_robust_WOA(para, H, user_r, user_theta, W_initial, theta_plot, r_plot, plot_flag)
%Inputs:
%   para: structure of the initial parameters
%   H: channel for all users
%   user_r: distance for all users
%   user_theta: angle for all users
%   W_initial: initialized digital beamformers (for FDA approach)
%   theta_plot: angle for plotting array gain (optional)
%   r_plot: distance for plotting array gain (optional)
%   plot_flag: flag to plot array gain (optional, default false)
%Outputs:
%   R: optimized spectral efficiency
%   A: optimized analog beamforming matrix
%   D: optimized digital beamforming matrix
%   t: optimized time delays of TTDs
%Date: 17/01/2026

%% Initialization
switch nargin
    case 4
    W_initial = randn(para.N, para.K) + 1i * randn(para.N, para.K);
    W_initial = W_initial / norm(W_initial, 'fro') * sqrt(para.Pt);
    plot_flag = false;
    case 5
    plot_flag = false;
    case 7
    plot_flag = false;
end

c = 3e8; % speed of light
t_search = 0:para.t_max/1e3:para.t_max; % search space of TTDs' time delay

N_sub = para.N/para.N_T; % number of antennas connected to each TTD
e = ones(N_sub, 1);

%% Analog beamformer design using WOA
% WOA parameters
SearchAgents_no = 30;
Max_iter = 100;
dim = para.N_T * para.N_RF + para.N * para.N_RF; % Time delays and PS coefficients (real and imag)
lb = [-10 * ones(1, para.N * para.N_RF), zeros(1, para.N_T * para.N_RF)];
ub = [10 * ones(1, para.N * para.N_RF), para.t_max * ones(1, para.N_T * para.N_RF)];

% Objective function
fobj = @(x) objective_analog(x, para, user_r, user_theta);

% Run WOA
[~, Best_pos, ~] = WOA(SearchAgents_no, Max_iter, lb, ub, dim, fobj);

% Extract t and A_PS
A_PS_real = Best_pos(1:para.N * para.N_RF);
A_PS_imag = Best_pos(para.N * para.N_RF + 1 : 2*para.N * para.N_RF);
t_vec = Best_pos(2*para.N * para.N_RF + 1 : end);

A_PS = reshape(A_PS_real + 1i * A_PS_imag, para.N, para.N_RF);
t = reshape(t_vec, para.N_T, para.N_RF);

A = zeros(para.N, para.N_RF, para.M);
H_equal = zeros(para.N_RF, para.K, para.M);
for m = 1:para.M
    A(:,:,m) = analog_bamformer(para, A_PS, t, para.fm_all(m)); % overall analog beamformer
    H_equal(:,:,m) = A(:,:,m)'*H(:,:,m); % equivalent channel
end

%% Digital beamformer design
D = zeros(para.N_RF, para.K, para.M);
for m = 1:para.M
    Dm = pinv(A(:,:,m))*W_initial;
    [~, Dm] = RWMMSE(para, H(:,:,m), H_equal(:,:,m), Dm, A(:,:,m));
    D(:,:,m) = Dm;
end

%% Calculate the spectral efficiency
W = zeros(para.N, para.K, para.M);
for m = 1:para.M
    W(:,:,m) = A(:,:,m)*D(:,:,m);
end

[R] = rate_fully_digital(para, W, H);
R = R/(para.M+para.Lcp);

%% Plot array gain if requested
if plot_flag && nargin >= 7
    [~, P_WOA_robust] = beampattern_WOA(para, theta_plot, r_plot, [], [], A, t);
    figure;
    plot(para.fm_all/1e9, 10*log10(P_WOA_robust/para.N), 'LineWidth', 1.5);
    xlabel('Frequency (GHz)', 'Interpreter', 'Latex');
    ylabel('Array Gain (dB)', 'Interpreter', 'Latex');
    title('Array Gain for WOA Robust Method', 'Interpreter', 'Latex');
    grid on;
end
end

%% Objective function for analog beamformer
function obj = objective_analog(x, para, user_r, user_theta)
    N_sub = para.N/para.N_T;
    e = ones(N_sub, 1);
    A_PS_real = x(1:para.N * para.N_RF);
    A_PS_imag = x(para.N * para.N_RF + 1 : 2*para.N * para.N_RF);
    t_vec = x(2*para.N * para.N_RF + 1 : end);
    
    A_PS = reshape(A_PS_real + 1i * A_PS_imag, para.N, para.N_RF);
    t = reshape(t_vec, para.N_T, para.N_RF);
    
    obj = 0;
    for n = 1:para.N_RF
        theta = user_theta(n); r = user_r(n);
        array_response = zeros(para.N, para.M);
        for m = 1:para.M
            fm = para.fm_all(m);
            bm = array_response_vector(r, theta, para.N, para.d, fm);
            array_response(:,m) = conj(bm);
        end
        
        gamma = zeros(para.N_T, para.M);
        for l = 1:para.N_T
            phi_n_l = A_PS(((l-1)*N_sub+1):l*N_sub, n);
            for m = 1:para.M
                gamma(l,m) = array_response( ((l-1)*N_sub+1):l*N_sub ,m)' * phi_n_l;
            end
        end
        
        for m = 1:para.M
            term = sum(gamma(:,m) .* exp(-1i*2*pi*para.fm_all(m)*t(:,n)));
            obj = obj - abs(term);
        end
    end
end

%% Calculate the overall analog beamformer at frequency f
function [A] = analog_bamformer(para, A_PS, t, f)
    e = ones(para.N/para.N_T,1);
    T = exp(-1i*2*pi*f*t);
    A = A_PS .* kron(T, e);
end

%% RWMMSE method for optimizing the digital beamformer
function [R, D] = RWMMSE(para, H, H_equal, D, A)
R_pre = 0;
for i = 1:20
    E = eye(para.K);
    Phi = 0; Upsilon = 0;
    for k = 1:para.K
        hk = H_equal(:,k);
        dk = D(:,k); 
        I = norm(hk'*D)^2 + norm(A*D, 'fro')^2/para.Pt; 
        w_k = 1 + abs(hk'*dk)^2 / (I - abs(hk'*dk)^2);
        v_k = hk'*dk / I;
    
        Phi = Phi + w_k*abs(v_k)^2 * ( hk*hk' + eye(para.N_RF)/para.Pt );
        Upsilon = Upsilon + w_k*conj(v_k)*E(:,k)*hk';
    
    end
    
    D = Phi\Upsilon';

    % check convergence
    [R] = rate_single_carrier(para, A*D, H);
    if abs(R - R_pre)/R <= 1e-4
        break;
    end
    R_pre = R;
end

end