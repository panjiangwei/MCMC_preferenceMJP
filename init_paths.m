function paths = init_paths(mus, sigmas, data, dt)

% every path has [tend; -1] as last dummy state
if nargin < 4
    dt = 0.1;
end

U = length(data);
paths = cell(1, U);

for u = 1:U
    tu = max(data{u}.t) + dt;
    obs_ll = construct_lik(mus, sigmas, data{u}.loc);
    n = length(data{u}.t);
    paths{u} = zeros(2, n);
    paths{u}(1, :) = data{u}.t;
    [~, paths{u}(2, :)] = max(obs_ll);

    % de-duplicate
    I = find(paths{u}(2, 1:end-1) == paths{u}(2, 2:end)) + 1;
    paths{u}(:, I) = [];
    paths{u} = [paths{u}, [tu; -1]];
end
