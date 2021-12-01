function S = make_duplication_matrix(T,tau)
  H = hankel(1:tau,tau:T);
  T2= prod(size(H));
  index = [(1:T2)' H(:)]-1;
  ind1  = 1 + index(:,1) + T2*index(:,2);
  S = sparse(T2,T);
  S(ind1) = 1;