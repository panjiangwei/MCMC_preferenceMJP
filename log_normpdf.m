function y = log_normpdf(X, mu, sigma)
% X: d*n
% mu: d*1
% sigma: d*d

[d, n] = size(X);

X = X - repmat(mu,1,n);
[R,p] = chol(sigma);
if p ~= 0
    error('ERROR: sigma is not PD.');
end
Q = R'\X;
q = dot(Q, Q, 1);
c = d*log(2*pi) + 2*sum(log(diag(R)));
y = -(q+c)/2;