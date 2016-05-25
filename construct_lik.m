function obs_ll = construct_lik(mus, sigmas, traj)
% return log likelihoods of each state

K = size(mus, 2);
n = size(traj, 2);

obs_ll = zeros(K, n);

for k = 1:K
    %obs_ll(k, :) = log(mvnpdf(traj', mus(:, k)', sigmas(:, :, k))');
    obs_ll(k, :) = log_normpdf(traj, mus(:, k), sigmas(:,:,k))';
%     for i = 1:n
%         obs_ll(k, i) = logGaussPdf(traj(:, i), mus(:, k), sigmas(:, :, k));
%     end
end