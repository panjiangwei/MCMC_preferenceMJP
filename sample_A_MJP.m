function A = sample_A_MJP(paths, K, a, alpha, beta)

U = length(paths);

% count(i, j): the # transitions from state i to j
count = zeros(K, K);

% nvisits(i): # times state i is visited
nvisit = zeros(1, K); 

% sumtime(i): total time spent in state i
sumtime = zeros(1, K);

for u = 1:U
    path = paths{u};
    states = path(2, :);
    times = path(1, :);
    for i = 1:length(states)-2
        curr_st = states(i);
        next_st = states(i+1);
        assert(curr_st ~= next_st);
        count(curr_st, next_st) = count(curr_st, next_st) + 1;
        nvisit(curr_st) = nvisit(curr_st) + 1;
        sumtime(curr_st) = sumtime(curr_st) + times(i+1) - times(i);
    end
end

% sample normalized transition probability
normA = gamrnd(count + a, 1);
normA(1:(K+1):end) = 0;
normA = normA ./ repmat(sum(normA, 2), 1, K);

% sample diagonal rates
diag_rates = gamrnd(alpha + nvisit, 1 ./ (beta + sumtime));

A = normA .* repmat(diag_rates', 1, K);
A = A - diag(diag_rates);



