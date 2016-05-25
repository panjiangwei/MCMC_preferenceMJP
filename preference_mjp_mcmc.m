% MCMC for preference-MJP model

clear

% contains cell array data = cell(U, 1)
% data{u}.loc: 2 * n vector, locations (lon, lat) of observations
% data{u}.t: 1 * n vector, times of observations (normalized in [0,1])
load sample_checkin_data.mat

% number of states
K = 50;

% data dimension
d = 2;

% number of MCMC iterations
num_iters = 50;

% number of users
U = length(data);

% hyper-params for Gamma prior of rate matrix elements A(i,j)
alpha_A = 0.5;
beta_A = 2;

% hyper-params for Gamma prior of personal weight vector elements
alpha_theta = 0.5;
beta_theta = 0.5;

% Gamma prior hyper-parameters for obs_rates
alpha_obs = 10; % shape
beta_obs = 0.5; % rate

% hyper-params for gaussians - notation same as wikipedia
all_loc = [];
for u = 1:U
    all_loc = [all_loc; data{u}.loc'];
end

gprior.mu = mean(all_loc)';
gprior.kappa = 2;
gprior.df = 5;
gprior.psi = cov(all_loc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize personal preference vectors
theta = ones(U, K);

% initialize global transition rate matrix
A = rand(K)/200;
A(1:K+1:end) = 0;
omega_factor = 2;

% initialize Gaussian parameters of the states
sigmas = zeros(d, d, K);                            
mus = zeros(d, K);
for k = 1:K
    sigmas(:, :, k) = iwishrnd(gprior.psi, gprior.df);
    mus(:, k) = mvnrnd(gprior.mu, sigmas(:, :, k)/gprior.kappa);
end

% initialize mjp paths
paths = init_paths(mus, sigmas, data);

% initialize observation rates
obs_rates = ones(K, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gibbs sampling 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obs_ll = cell(U, 1);
tA = cell(U, 1);
tA_dom = cell(U, 1);

for iter = 1:num_iters
    fprintf(1, 'starting iteration %d\n', iter);

    % sample mjp paths
    parfor u = 1:U
        obs_ll{u} = construct_lik(mus, sigmas, data{u}.loc);
        tA{u} = A .* repmat(theta(u, :), K, 1);
        tA{u}(1:K+1:end) = - sum(tA{u}, 2);
        tA_dom{u} = -omega_factor .* diag(tA{u});
        paths{u} = sample_mjp_path(paths{u}, tA{u}, tA_dom{u}, data{u}.t, obs_ll{u}, obs_rates);
    end

    % sample global transition rate matrix
    A = sample_A_preferenceMJP(paths, theta, alpha_A, beta_A);

    % sample user weight vectors
    theta = sample_theta_preferenceMJP(A, paths, alpha_theta, beta_theta);

    % sample state parameters
    [mus, sigmas] = sample_gaussians(paths, data, gprior, K);
    
    % sample observation rates
    obs_rates = sample_obs_rates(paths, data, alpha_obs, beta_obs, K);
end
