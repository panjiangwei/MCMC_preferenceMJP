function [mus, sigmas] = sample_gaussians(paths, data, gprior, K)

% paths{i}: 2*ni
% data{i}.loc: 2*mi (lon, lat)
% data{i}.t: 1*mi

d = 2;
U = length(paths);

% sufficient statistic 1: sum of observations
sum_obs = zeros(d, K);

% sufficient statistic 2: sum of squares of observations
sum_obs2 = zeros(d, d, K);

% sufficient statistic 3: # of observations
num_obs = zeros(1, K);

for u = 1:U
    path = paths{u};
    states = path(2, :);
    times = path(1, :);

    obs_loc = data{u}.loc;
    obs_times = data{u}.t;
    
    j = 1;
    i = 1;
    while j < length(times) && i <= size(obs_loc,2)
        curr_st = states(j);
        tlimit = times(j+1);
        if obs_times(i) >= tlimit
            j = j + 1;
            continue;
        end

        % add the i-th observation to ss of curr_st
        num_obs(curr_st) = num_obs(curr_st) + 1;
        sum_obs(:, curr_st) = sum_obs(:, curr_st) + obs_loc(:, i);
        sum_obs2(:, :, curr_st) = sum_obs2(:, :, curr_st) + obs_loc(:, i) * obs_loc(:, i)';
        i = i + 1;
    end
end

mus = zeros(d, K);
sigmas = zeros(d, d, K);
for k = 1:K
    gpost = gprior;
    if num_obs(k) > 0
        obs_bar = sum_obs(:, k) / num_obs(k);
        
        gpost.kappa = gprior.kappa + num_obs(k);
        gpost.df = gprior.df + num_obs(k);
        gpost.mu = (gprior.kappa * gprior.mu + ...
            sum_obs(:, k)) / (gpost.kappa);
        gpost.psi = gprior.psi + sum_obs2(:, :, k) ...
            - num_obs(k) * (obs_bar * obs_bar') ...
            + (gprior.kappa * num_obs(k) / gpost.kappa) ...
            * (obs_bar - gprior.mu) * (obs_bar - gprior.mu)';
        gpost.psi = (gpost.psi + gpost.psi')/2;
    end

    % sample gaussian parameters from posterior
    sigmas(:, :, k) = iwishrnd(gpost.psi, gpost.df);
    mus(:, k) = mvnrnd(gpost.mu, sigmas(:, :, k) / gpost.kappa);
end

