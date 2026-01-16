function [R, A, D, t] = algorithm_HTS_PNF_WOA(para, H, user_r, user_theta, W_initial)
%The piecewise-near-field (PNF)-based heuristic two-stage approach using WOA
%  [R, A, D, t] = algorithm_HTS_PNF_WOA(para, H, user_r, user_theta, W_initial)
%Inputs:
%   para: structure of the initial parameters
%   H: channel for all users
%   user_r: distance for all users
%   user_theta: angle for all users
%   W_initial: initialized digital beamformers (for FDA approach)
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
end

c = 3e8; % speed of light
m = 1:para.M;
para.fm_all = para.fc + para.B*(2*m-1-para.M) / (2*para.M); % subcarrier frequencies
t_search = 0:para.t_max/1e3:para.t_max; % search space of TTDs' time delay

%% Analog beamformer design using WOA
% WOA parameters
SearchAgents_no = 30;
Max_iter = 100;
dim = para.N_T * para.N_RF; % Time delays
lb = zeros(1, dim);
ub = para.t_max * ones(1, dim);

% Objective function
fobj = @(x) objective_analog(x, para, user_r, user_theta);

% Run WOA
[~, Best_pos, ~] = WOA(SearchAgents_no, Max_iter, lb, ub, dim, fobj);

t = reshape(Best_pos, para.N_T, para.N_RF);

% Calculate PS coefficients and A
A_PS = zeros(para.N, para.N_RF);
for n = 1:para.N_RF
    theta = user_theta(n); r = user_r(n);
    N_sub = para.N/para.N_T;
    a_n = zeros(para.N, 1);
    for l = 1:para.N_T
        xi_l = (l-1-(para.N_T-1)/2)*N_sub;
        r_l = sqrt(r^2 + xi_l^2*para.d^2 - 2*r*xi_l*para.d*cos(theta));
        theta_l = acos( (r*cos(theta) - xi_l*para.d)/r_l );
        q = (0:(N_sub-1))';
        delta_q = (q-(N_sub-1)/2) * para.d;
        a_n((l-1)*N_sub+1 : l*N_sub) = exp( 1i * 2 * pi * para.fc/c...
            * (sqrt(r_l^2 + delta_q.^2 - 2*r_l*delta_q*cos(theta_l)) - r_l) );
    end
    A_PS(:,n) = a_n;
end

A = zeros(para.N, para.N_RF, para.M);
H_equal = zeros(para.N_RF, para.K, para.M);
for m = 1:para.M
    A(:,:,m) = analog_bamformer(para, A_PS, t, para.fm_all(m));
    H_equal(:,:,m) = A(:,:,m)'*H(:,:,m);
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
end

%% Objective function for analog beamformer
function obj = objective_analog(x, para, user_r, user_theta)
    t = reshape(x, para.N_T, para.N_RF);
    c = 3e8;
    obj = 0;
    for n = 1:para.N_RF
        theta = user_theta(n); r = user_r(n);
        N_sub = para.N/para.N_T;
        r_n = zeros(para.N_T, 1);
        for l = 1:para.N_T
            xi_l = (l-1-(para.N_T-1)/2)*N_sub;
            r_l = sqrt(r^2 + xi_l^2*para.d^2 - 2*r*xi_l*para.d*cos(theta));
            r_n(l) = r_l;
        end
        for m = 1:para.M
            fm = para.fm_all(m);
            term = sum(exp(-1i*2*pi*fm*( (r_n - r)/c + t(:,n) )));
            obj = obj - abs(term); % Maximize sum abs, so minimize negative
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