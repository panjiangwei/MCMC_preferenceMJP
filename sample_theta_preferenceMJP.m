function theta = sample_theta_preferenceMJP(A, paths, a_theta, b_theta)

% number of states
K = size(A, 1);

% number of users
U = length(paths);

% numto(u, j) - num of times each user u transits to state j
%               from some state i
numto = zeros(U, K);

% sumtime(u, i) - total time each user u stays in state i
sumtime = zeros(U, K);

% count
for u = 1:U
    st_seq = paths{u}(2, :);
    t_seq = paths{u}(1, :);

    for i = 1:length(st_seq)-2
        st = st_seq(i);
        nst = st_seq(i+1);
        numto(u, nst) = numto(u, nst) + 1;
        sumtime(u, st) = sumtime(u, st) + t_seq(i+1) - t_seq(i);
    end
end

% sample posterior
A_theta = a_theta + numto;
B_theta = b_theta + sumtime * A;
theta = gamrnd(A_theta, 1 ./ B_theta);
