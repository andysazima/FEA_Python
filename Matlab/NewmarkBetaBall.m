function [d, v, a, algorithm_name] = NewmarkBetaBall(M,C,K,F,do,vo,beta,gamma,dt,t_total,Po,Ri)
% Algorithm Naming
switch beta
    case 0
        algorithm_name = 'Central Difference';
    case 1/12
        algorithm_name = 'Fox-Goodwin';
    case 1/6
        algorithm_name = 'Linear Acceleration';
    case 0.25
        algorithm_name = 'Average Acceleration';
    otherwise
        algorithm_name = 'Unknown Algorithm';
end

% Check for explicit vs implicit
explicit_check = 0;
M_check1 = 1 ./ M;
M_check2 = inv(M);
if beta == 0 && M_check1(end,1) == M_check2(end,1)
    explicit_check = 1;
    M_inv = 1 ./ M;
    % Explicit Initial Acceleration
    ao = M_inv * (F - C*vo - K*do);
else
    % Implicit Initial Acceleration
    ao = M \ (F - C*vo - K*do);
end
% Initialize Algorithm
t = 0:dt:t_total;
d = [do zeros(size(K,1),length(t)-1)];
v = [vo zeros(size(K,1),length(t)-1)];
a = [ao zeros(size(K,1),length(t)-1)];

if explicit_check == 1
    % EXPLICIT METHOD
    for n = 1:length(t)-1
        % Predictor Step
        d_bar = d(:,n) + v(:,n) * dt + 0.5*(1 - 2*beta) * dt^2 * a(:,n);
        v_bar = v(:,n) + (1 - gamma) * dt * a(:,n);
        % Force Calc Step
        Fi = 4*pi*(Ri + d(1,n))^2*(Po * Ri^3.75) / (Ri + d(1,n))^3.75;
        Force_BC = [1 Fi];
        F = ForceVector(K,Force_BC);
        % Solution Step
        F_star = F - C * v_bar - K * d_bar;
        M_star = M + gamma * dt * C;
        M_star_inv = 1 ./ M_star;
        a(:,n+1) = M_star_inv * F_star;
        % Corrector Step
        d(:,n+1) = d_bar + beta * dt^2 * a(:,n+1);
        v(:,n+1) = v_bar + gamma * dt * a(:,n+1);
    end
else
    % IMPLICIT METHOD
    for n = 1:length(t)-1
        % Predictor Step
        d_bar = d(:,n) + v(:,n) * dt + 0.5*(1 - 2*beta) * dt^2 * a(:,n);
        v_bar = v(:,n) + (1 - gamma) * dt * a(:,n);
        % Force Calc Step
        Fi = 4*pi*(Ri + d(1,n))^2*(Po * Ri^3.75) / (Ri + d(1,n))^3.75;
        Force_BC = [1 Fi];
        F = ForceVector(K,Force_BC);
        % Solution Step
        F_star = F - C * v_bar - K * d_bar;
        M_star = M + gamma * dt * C + beta * dt^2 * K;
        a(:,n+1) = M_star \ F_star;
        % Corrector Step
        d(:,n+1) = d_bar + beta * dt^2 * a(:,n+1);
        v(:,n+1) = v_bar + gamma * dt * a(:,n+1);
    end
end