% calculate Proximal Projection
% x = argmin_x || x - a ||_2^2, s.t. || q .* (f - x) ||_2^2 <= delta
%
function x = prox_indicator(a,f,q,delta)

  II = size(a);
  x  = zeros(II);

  r  = sqrt(delta)/norm((f(q(:)==1)-a(q(:)==1)));
  qt = q*max(0,1-r);
  x  = qt.*f + (-qt+1).*a;

  % part I: case for q(i)==0
  %x(q==0) = a(q==0);

  % part II: case for q(i)==1
  %r = sqrt(delta)/norm((f(q(:)==1)-a(q(:)==1)));
  %x(q==1) = max(1-r,0)*f(q==1) + min(r,1)*a(q==1);


