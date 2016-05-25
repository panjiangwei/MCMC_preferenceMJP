function A = sample_A_preferenceMJP(paths, theta, alpha, beta)

% number of users
U = length(paths);

% number of clusters
K = size(theta, 2);

% count(i, j): the # transitions from state i to j
count = zeros(K, K);

% sumtime(u, i): total time user u spent in state i
sumtime = zeros(U, K);

for u = 1:U
    states = paths{u}(2, :);
    times = paths{u}(1, :);
    for i = 1:length(states)-2
        curr_st = states(i);
        next_st = states(i+1);
        assert(curr_st ~= next_st);
        count(curr_st, next_st) = count(curr_st, next_st) + 1;
        sumtime(u, curr_st) = sumtime(u, curr_st) ...
                            + (times(i+1) - times(i));
    end
end

% sample A
A = gamrnd(count + alpha, 1 ./ (sumtime' * theta + beta));
A(1:K+1:end) = 0;
