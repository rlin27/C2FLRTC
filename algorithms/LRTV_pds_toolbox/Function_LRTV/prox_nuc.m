% calculate Proximal Projection
% X_opt = argmin_X lam|| X ||_nuclear_norm + 0.5|| X - Y ||_2^2
%
function X = prox_nuc(Y,lam)

  [II JJ] = size(Y);
  if II < JJ
    [U D v] = svd(Y*Y');
    D = sqrt(D);
    V = Y'*(U*diag(diag(D).^(-1)));
  else
    [u D V] = svd(Y'*Y);
    D = sqrt(D);
    U = Y*(V*diag(diag(D).^(-1)));
  end
  Ds = max(D - lam,0);
  X = U*Ds*V';

