function [ X ] = DE_tensor_all_adjoint( V, S )
  % Adjoint Operation of Delay Embedding for all directions (modes)
  % 
  N2  = ndims(V);
  N   = N2/2;
  JJ  = size(V);
  if mod(length(JJ),2)
      JJ = [JJ 1];
  end
  JJ2 = JJ(1:2:end-1) .* JJ(2:2:end);
  if length(JJ2) == 1
      JJ2 = [JJ2 1];
  end
  V   = reshape(V, JJ2);
  X = tensor_allprod(V,S,1,size(V));
