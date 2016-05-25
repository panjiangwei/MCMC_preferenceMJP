function x = sampleDiscrete(Px,n)
% x = sampleDiscrete(Px,n)
%
% sampleDiscrete returns a set of sample from a discrete probability distribution
%
% Px :  Probability distribution from which to draw.
%       It can either be a vector if all samples should be drawn from the same distribution,
%       or a matrix of column vectors each containing a different distribution. All distributions
%       must have the same cardinallity
%
% n :   The number of samples to draw. If Px is a matrix, draw n samples from each distribution

% Added by Vinayak

  if (sum(sum(isnan(Px))))
    error "Probability vector has NaNs!"
  end;

  if (isvector(Px))
     Px = Px(:);
     if (abs(sum(Px) - 1) > 1E-10), error('Px is not a probability distribution'); end;
  else
     if(sum(abs(sum(Px,1) - ones(1,size(Px,2)))) > size(Px,2) * 10^-10), error('Px is not a probability distribution'); end;
  end;

  %Need to do a case split, as the octave permute function has a different semantics
  % it can't permute the 3rd dimension if there is no 3rd dimension, as matlab can (which
  % assumes that dimension has size 1)
  if (n == 1),
    i = repmat(rand(size(Px,2),n),[1, size(Px,1)])';
    j = cumsum(Px,1);
    tmp = (i <= j);
    x = repmat(size(Px,1) + 1, [1,size(Px,2)]) - sum(tmp,1);
  else
    i = permute(repmat(rand(size(Px,2),n),[1, 1, size(Px,1)]), [3 1 2]);
    j = repmat(cumsum(Px,1),[1, 1, n]);
    tmp = (i <= j);
    x = permute(repmat(size(Px,1) + 1, [1,size(Px,2),n]) - sum(tmp,1), [3 2 1]);
  end;

