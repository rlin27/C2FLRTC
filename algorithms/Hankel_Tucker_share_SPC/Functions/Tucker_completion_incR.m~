%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X, G, U, histo] = Tucker_completion_incR(T,Q,G,U,delta,incR,maxiter,inloop,tol,verb)
%
% This program solves Tucker decomposition of input incomplete tensor by rank increment approach.
%
% min  rank_1(X) + rank_2(X) + ... + rank_N(X),
% s.t. X = G * U{1} * U{2} * ... * U{N},
%      || Q.*(T - X) ||_F^2 / |Q| <= delta 
%
% inputs:
%   -- T       : input incomplete tensor
%   -- Q       : mask tensor, 0:missing, 1:available
%   -- G       : initialization of core tensor of Tucker decomposition
%   -- U       : initialization of (1xN)-cell array consisting of factor matrices
%   -- delta   : noise threshold parameter
%   -- incR    : (1xN)-cell array consisting of rank sequences for increment
%   -- maxiter : maximum number of iterations
%   -- inloop  : number of iterations for inner loop
%   -- tol     : tolerance parameter for checking convergence
%   -- verb    : verbosity for visualizing process of algorithm
%
% outputs:
%   -- X     : output complete tensor
%   -- G     : result of core tensor of Tucker decomposition
%   -- U     : result of (1xN)-cell array consisting of factor matrices
%   -- histo : history of || Q.*(T - X) ||_F^2 / |Q| recorded for each iteration
%   -- histoR: history of multilinear tensor ranks recorded for each iteration
%
% This code was written by Tatsuya Yokota (2017.08.28)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, G, U, histo, histoR] = Tucker_completion_incR(T,Q,G,U,delta,incR,maxiter,inloop,tol,verb)

  N  = length(U);
  for n = 1:N
    II(n) = size(U{n},1);
    R(n) = size(U{n},2);
  end

  T = Q.*T;
  GU= tensor_allprod(G,U,0,R);

  obj = (1/sum(Q(:)))*norm(T(Q(:)==1) - GU(Q(:)==1))^2;

  Z  = T;
  Ri = ones(1,N);
  for n = 1:N
    R(n) = incR{n}(1);
  end
  for iter = 1:maxiter

    % update parameters
    Z(Q(:)~=1) = GU(Q(:)~=1);
    for iter2 = 1:inloop
    for n = 1:N
      Y{n} = unfold(tensor_allprod_exc(Z,U,1,n,II),n);
      [U{n},~,~] = svds(Y{n}*Y{n}',R(n));
    end
    end
    G = tensor_allprod(Z,U,1,II);
    GU= tensor_allprod(G,U,0,R);

    % calc. cost function
    obj2 = (1/sum(Q(:)))*norm(T(Q(:)==1) - GU(Q(:)==1))^2;
    histo(iter) = obj2;
    histoR(iter,:) = R;

    % show process
    if mod(iter,verb) == 0
      fprintf('iter %d:: cost = %e :: cost_diff = %e :: R=[',iter,obj2,abs(obj2-obj));
      for n = 1:N
        fprintf('%d ',R(n));
      end
      fprintf('] \n');
    end

    % check rank increment condition
    if iter > 3 && abs(obj2 - obj) < tol
      
      % calc. mode residual
      dX = Q.*(T-GU);
      for n = 1:N
        E{n} = tensor_allprod_exc(dX,U,1,n,II);
        mode_residual(n) = norm(E{n}(:));
      end
      [value idd] = max(mode_residual);
      Ri(idd)= min(length(incR{idd}),Ri(idd) + 1);
      R(idd) = min(II(idd),incR{idd}(Ri(idd)));
      
    else
      obj = obj2;
    end

    % check termination condition
    if obj2 < delta
      break;
    end

  end
  X = GU;

