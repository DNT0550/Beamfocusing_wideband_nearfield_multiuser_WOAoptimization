function [Best_score, Best_pos, Convergence_curve] = WOA(SearchAgents_no, Max_iter, lb, ub, dim, fobj)
% Whale Optimization Algorithm (WOA) for optimization problems
% Inputs:
%   SearchAgents_no: Number of search agents
%   Max_iter: Maximum number of iterations
%   lb: Lower bound of variables
%   ub: Upper bound of variables
%   dim: Dimension of the problem
%   fobj: Objective function handle
% Outputs:
%   Best_score: Best objective value
%   Best_pos: Best solution
%   Convergence_curve: Convergence history

% Initialize the positions of search agents
Positions = initialization(SearchAgents_no, dim, ub, lb);

% Initialize convergence curve
Convergence_curve = zeros(1, Max_iter);

% Loop counter
iter = 0;

% Main loop
while iter < Max_iter
    for i = 1:size(Positions, 1)
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub = Positions(i, :) > ub;
        Flag4lb = Positions(i, :) < lb;
        Positions(i, :) = (Positions(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

        % Calculate objective function for each search agent
        fitness(i) = fobj(Positions(i, :));
    end

    % Update the leader
    [Best_score, index] = min(fitness);
    Best_pos = Positions(index, :);

    % Update convergence curve
    Convergence_curve(iter + 1) = Best_score;

    % a decreases linearly from 2 to 0
    a = 2 - iter * (2 / Max_iter);

    % a2 linearly decreases from -1 to -2
    a2 = -1 + iter * (-1 / Max_iter);

    for i = 1:size(Positions, 1)
        r1 = rand();
        r2 = rand();

        A = 2 * a * r1 - a;
        C = 2 * r2;

        b = 1;
        ll = (a2 - 1) * rand + 1;

        p = rand();

        for j = 1:size(Positions, 2)
            if p < 0.5
                if abs(A) >= 1
                    rand_leader_index = floor(SearchAgents_no * rand() + 1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand = abs(C * X_rand(j) - Positions(i, j));
                    Positions(i, j) = X_rand(j) - A * D_X_rand;
                elseif abs(A) < 1
                    D_Leader = abs(C * Best_pos(j) - Positions(i, j));
                    Positions(i, j) = Best_pos(j) - A * D_Leader;
                end
            elseif p >= 0.5
                distance2Leader = abs(Best_pos(j) - Positions(i, j));
                Positions(i, j) = distance2Leader * exp(b * ll) * cos(ll * 2 * pi) + Best_pos(j);
            end
        end
    end

    iter = iter + 1;
end
end

function Positions = initialization(SearchAgents_no, dim, ub, lb)
Boundary_no = size(ub, 2); % number of boundaries

% If the boundaries of all variables are equal and user enters a signle
% number for both ub and lb
if Boundary_no == 1
    Positions = rand(SearchAgents_no, dim) .* (ub - lb) + lb;
end

% If each variable has a different lb and ub
if Boundary_no > 1
    for i = 1:dim
        ub_i = ub(i);
        lb_i = lb(i);
        Positions(:, i) = rand(SearchAgents_no, 1) .* (ub_i - lb_i) + lb_i;
    end
end
end</content>
<parameter name="filePath">c:\Users\ADMIMN\Desktop\beamfocusing-optimization-for-near-field-wideband-multi-user-communications\fully_connected\functions\WOA.m