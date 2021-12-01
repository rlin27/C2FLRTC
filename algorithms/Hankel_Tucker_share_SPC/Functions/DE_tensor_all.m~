function [ V S ] = DE_tensor_all( X, tau, S)
  % Delay Embedding for all directions (modes)
  % 
  N  = ndims(X);
  II = size(X);

  if isempty(S)
    for n = 1:N
      S{n} = make_duplication_matrix(II(n),tau(n));
    end
  end
  
  I2= II - tau + 1;
  JJ= [tau; I2]; JJ = JJ(:);
  V = tensor_allprod(X,S,0,size(X));
  V = reshape(full(V),[JJ']);

end
