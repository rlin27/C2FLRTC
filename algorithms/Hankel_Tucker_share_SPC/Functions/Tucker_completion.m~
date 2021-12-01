%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X, G, U, histo] = Tucker_completion(T,Q,G,U,maxiter,inloop,tol,verb)
%
% This program solves fixed rank Tucker decomposition of input incomplete tensor.
%
% min  || Q.*(T - X) ||_F^2
% s.t. X = G * U{1} * U{2} * ... * U{N},
%      size(G) = (R_1, R_2, ..., R_N),
%
% inputs:
%   -- T       : input incomplete tensor
%   -- Q       : mask tensor, 0:missing, 1:available
%   -- G       : initialization of core tensor of Tucker decomposition
%   -- U       : initialization of (1xN)-cell array consisting of factor matrices
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
%
% This code was written by Tatsuya Yokota (2017.08.28)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, G, U, histo] = Tucker_completion(T,Q,G,U,maxiter,inloop,tol,verb)

  N  = length(U);
  for n = 1:N
    II(n) = size(U{n},1);
    R(n) = size(U{n},2);
  end

  T = Q.*T;
  GU= tensor_allprod(G,U,0,R);

  obj = (1/sum(Q(:)))*norm(T(Q(:)==1) - GU(Q(:)==1))^2;

  Z = T;
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
    
    % show process
    if mod(iter,verb) == 0
      fprintf('iter %d:: cost = %e :: cost_diff = %e \n',iter,obj2,abs(obj2-obj));
    end

    % convergence check
    if abs(obj2 - obj) < tol
      break;
    else
      obj = obj2;
    end

  end
  X = GU;


