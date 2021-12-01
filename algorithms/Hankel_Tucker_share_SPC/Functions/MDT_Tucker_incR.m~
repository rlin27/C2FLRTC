%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Xest, histo, histoR, G, U, S, D] = MDT_Tucker_incR(T,Q,tau,param)
%
% This program solves tensor completion problem via Tucker decomposition in multi-way delayembedded space.
%
% Simple flow
%  - Step 1. multi-way delay embedding transform (MDT): V = MDT(T); (N-th-order tensor --> M-th-order tensor)
%  - Step 2. Tucker decomposition                     : Vest = Tucker(V); (Tucker decomp. of M-th-order tensor)
%  - Step 3. inverse MDT                              : Xest = invMDT(Vest); (M-th-order tensor --> N-th-order tensor)
%
%
% Inputs:
%   -- T       : input incomplete tensor: (N-th-order tensor)
%   -- Q       : mask tensor, 0:missing, 1:available (N-th-order tensor)
%   -- tau     : delay-embedding parameter: (N-dimensional vector)
%   -- param   : optional parameters
%       |--> param.delta     : noise threshold parameter
%       |--> param.incR      : (1xM)-cell array consisting of rank sequences for increment
%       |--> param.maxiter   : maximum number of iterations
%       |--> param.inloop    : number of iterations for inner loop
%       |--> param.tol       : tolerance parameter for checking convergence
%       |--> param.verb      : verbosity for visualizing process of algorithm
%
% Outputs:
%   -- X     : output complete tensor
%   -- histo : history of cost function recorded for each iteration
%   -- histoR: history of multilinear tensor ranks recorded for each iteration
%   -- G     : result of core tensor of Tucker decomposition
%   -- U     : result of (1xM)-cell array consisting of factor matrices
%   -- S     : (1xM)-cell array consisting of duplication matrices used for MDT
%   -- D     : tensor of number of duplications used for inverse MDT
%
% This code was written by Tatsuya Yokota (2017.08.28)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xest, histo, histoR, G, U, S, D] = MDT_Tucker_incR(T,Q,tau,param)

  % check optional parameters
  if nargin == 3
    delta   = 1e-8;
    incR    = [];
    maxiter = 10000;
    inloop  = 1;
    tol     = 1e-4;
    verb    = 1;
  else
    if isfield(param,'delta'),   delta   = param.delta;   else, delta   = 1e-8;  end
    if isfield(param,'incR') ,   incR    = param.incR;    else, incR    = [];    end
    if isfield(param,'maxiter'), maxiter = param.maxiter; else, maxiter = 10000; end
    if isfield(param,'inloop'),  inloop  = param.inloop;  else, inloop  = 1;     end
    if isfield(param,'tol'),     tol     = param.tol;     else, tol     = 1e-4;  end
    if isfield(param,'verb'),    verb    = param.verb;    else, verb    = 1;     end
  end

  N   = ndims(T);
  Ni  = 2*N;
  II  = size(T);

  if length(tau) == 1
    tau = [tau 1];
    flag = 1;
  else
    flag = 0;
  end
  
  [v0, S] = DE_tensor_all(ones(II),tau,[]);
  II2 = size(v0);
  if mod(length(II2),2)
    II2 = [II2 1];
  end
  if flag
    II2 = [II2 1 1];
  end
  D  = DE_tensor_all_adjoint(v0,S);

  % Step1 : Multi-way delay embedding transform (MDT)
  Tms  = zeros(II);
  Tms(Q(:)==1) = T(Q(:)==1);
  QQ   = DE_tensor_all(Q,tau,S);
  V    = DE_tensor_all(Tms,tau,S);

  % Step2 : Tucker decomposition & its initialization
  if isempty(incR)
    for n = 1:Ni
      incR{n} = [1:II2(n)];
    end
  end
  for n = 1:Ni
    R(n) = incR{n}(1);
    U{n} = orth(randn(II2(n),R(n)));
  end
  G = randn(R);

  [Vest, G, U, histo, histoR] = Tucker_completion_incR(V,QQ,G,U,delta,incR,maxiter,inloop,tol,verb);

  % Step3 : inverse MDT
  Xest = DE_tensor_all_adjoint(Vest,S) ./ D;

