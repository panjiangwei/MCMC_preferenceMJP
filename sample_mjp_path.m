function mjp_path = sample_mjp_path(samples_old, A, A_dom, obs_time, obs_ll, obs_rates)

K = length(A_dom);

A_diff = (A_dom + diag(A));

T = eye(K) + A ./ A_dom(:,ones(1,K));
T = T';                  % T(new, old)
log_T = log(T);


num_obs = length(obs_time);


p_init = ones(K,1);
p_init = p_init ./ sum(p_init);


% Sample candidate jump times given MJP trajectory
[R, r_indx] = get_candt_times(samples_old, A_diff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now sample MJP trajectory given candidate jump times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obs_indx = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg(K, r_indx) = 0;

for m = 1:r_indx-1
    
    poiss_prob =  -(R(m+1)-R(m)) .* A_dom;     % (\int_{R(m)}^{R(m+1)} u(t)dt
    
    if(m < r_indx-1)   % Because the last event is not a Poisson event
        poiss_prob = poiss_prob + log(A_dom);
    end
    
    if(m == 1)
        msg(:,m) = p_init;
    else
        msg(:,m) = T * exp(msg(:,m-1));
    end
    msg(:,m) = log(msg(:,m)) + poiss_prob; 
    
    %%%%%%%%%%%%%%%%%%%%%
    % Problem specific observation likelihood
    obs_count = 0;
    while(obs_indx <= num_obs && obs_time(obs_indx) < R(m+1))
        msg(:,m) = msg(:,m) + obs_ll(:,obs_indx);
        obs_indx = obs_indx + 1;
        obs_count = obs_count + 1;
    end
    msg(:,m) =  msg(:, m) - (R(m+1)-R(m)) .* obs_rates + ...
        obs_count .* log(obs_rates .* (R(m+1) - R(m)));
    %%%%%%%%%%%%%%%%%%%%%
    
    msg(:,m) = msg(:,m) - max(msg(:,m));
    
    if sum(isnan(msg(:, m))) > 0
        error('msg contains NaN');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Backward sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_curr   = r_indx-1;

p_curr   = msg(:,t_curr);
p_curr   = p_curr - logsumexp_v(p_curr);
s_curr   = sampleDiscrete(exp(p_curr),1);

s_old    = s_curr;
mjp_path = [R(t_curr); s_curr];

while(t_curr > 1)
    t_curr = t_curr - 1;
    
    p_curr   = msg(:,t_curr) + log_T(s_curr,:)';
    p_curr   = p_curr - logsumexp_v(p_curr);
    s_curr   = sampleDiscrete(exp(p_curr),1);
    
    if(s_curr == s_old)
        mjp_path(1,1) = R(t_curr);
    else
        mjp_path = [R(t_curr) mjp_path(1,:);
            s_curr mjp_path(2,:)];
        s_old    = s_curr;
    end
end

mjp_path = [mjp_path(1,:) R(end); mjp_path(2,:) -1];

