% visualize the MJP states by plotting observations assigned to each state

% derive the state of each observation in the last MCMC sample
N = length(all_loc);
q = zeros(1, N);

qid = 1;
for u = 1:U
    path = paths{u};
    states = path(2, :);
    times = path(1, :);
    
    obs = data{u}.loc;
    obs_times = data{u}.t;
    
    j = 1;
    i = 1;
    while j < length(times) && i <= size(obs,2)
        curr_st = states(j);
        tlimit = times(j+1);
        if obs_times(i) >= tlimit
            j = j + 1;
            continue;
        end
        
        q(qid) = curr_st;
        i = i + 1;
        qid = qid + 1;
    end
end

% plot states
h = figure;
nrow = 7;
ncol = 8;
assert(K <= nrow * ncol);
sampleInds = randi(N, 1, 1000);
for k = 1:K
    subplot(nrow, ncol, k);
    hold on
    
    lon = all_loc(:, 1);
    lat = all_loc(:, 2);
    
    plot(all_loc(sampleInds, 1), all_loc(sampleInds, 2), '.k', 'MarkerSize', 1)
    xlim([min(lon), max(lon)])
    ylim([min(lat), max(lat)])
    title(sprintf('%d - %d', k, sum(q==k)))
    J = find(q==k);
    plot(lon(J), lat(J), '.r', 'MarkerSize', 5)
end