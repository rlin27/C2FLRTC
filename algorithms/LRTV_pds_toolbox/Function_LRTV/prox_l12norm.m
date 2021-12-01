% calculate Proximal Projection
% X_opt = argmin_X lam|| X ||_(1,2) + 0.5|| X - Y ||_2^2
%
function X = prox_l12norm(Y,lam)

  [N K] = size(Y);
  w = sqrt(sum(Y.^2,2));
  
  s = max(ones(N,1) - lam*(ones(N,1) ./ w) , 0);
  X = bsxfun(@times,Y,s);
