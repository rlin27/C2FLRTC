function x = prox_dom(a,dom)

  x = a;
  x = max(dom(1),x);
  x = min(dom(2),x);
