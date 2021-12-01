function Z = tensor_allprod_exc(G,U,tr,exc,Dg)

  N = length(Dg);

  Z = G;
  for n = 1:N
    if ~isempty(U{n}) && n~=exc
      if tr == 0
        Z = tmult(Z,U{n},n,Dg);
        Dg(n) = size(U{n},1);
      else
        Z = tmult(Z,U{n}',n,Dg);
        Dg(n) = size(U{n},2);
      end
    end
  end
