% MCMC for Markov Jump Process (MJP)

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
num_iters = 1000;

% number of users
U = length(data);

% hyper-params for gaussians - notation same as wikipedia
all_loc = [];
for u = 1:U
    all_loc = [all_loc; data{u}.loc'];
end

gprior.mu = mean(all_loc)';
gprior.kappa = 2;
gprior.df = 5;
gprior.psi = cov(all_loc);

% Gamma prior hyper-parameters for observation rates
alpha_obs = 2;
beta_obs = alpha_obs / 10;

% Gamma prior hyper-parameters for |A(i,i)|
alpha_A = 1;
beta_A = 1;

% Dirichlet prior hyper-parameter for transition probabilities
a = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialization transition rate matrix
A = rand(K)/20;
A(1:K+1:end) = 0;
A(1:K+1:end) = -sum(A, 2);

omega_factor = 2;
A_dom = -omega_factor .* diag(A);

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

obs_ll = cell(1, U);

for iter = 1:num_iters
    fprintf(1, 'starting iteration %d\n', iter);
    
    % sample mjp paths
    parfor u = 1:U
        obs_ll{u} = construct_lik(mus, sigmas, data{u}.loc);
        paths{u} = sample_mjp_path(paths{u}, A, A_dom, data{u}.t, obs_ll{u}, obs_rates);
    end

    % sample transition rate matrix
    A = sample_A_MJP(paths, K, a, alpha_A, beta_A);
    A_dom = -omega_factor .* diag(A);

    % sample observation distribution
    [mus, sigmas] = sample_gaussians(paths, data, gprior, K);

    % sample observation rates
    obs_rates = sample_obs_rates(paths, data, alpha_obs, beta_obs, K);
end
