function obs_rates = sample_obs_rates(paths, data, alpha_obs, beta_obs, K)

% sum_counts(k) = total # of observations with state k
sum_counts = zeros(1, K);

% sum_time(k) = total time in state k
sum_time = zeros(1, K);

U = length(data);
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
            sum_time(curr_st) = sum_time(curr_st) + tlimit - times(j);
            j = j + 1;
            continue;
        end

        sum_counts(curr_st) = sum_counts(curr_st) + 1;
        i = i + 1;
    end
end

obs_rates = gamrnd(alpha_obs + sum_counts, 1 ./ (beta_obs + sum_time))';
