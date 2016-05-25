function y=logsumexp_v(x)
% y0 = log(sum(exp(x)))
% only preventing overflow/underflow

  xmax = max(x);
  if(xmax > -inf)
    x = x-xmax;
  end;
  y = xmax + log(sum(exp(x)));
  return

