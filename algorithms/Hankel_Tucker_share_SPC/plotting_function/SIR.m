function sir = SIR(Ori,Est)

  s = double(Ori(:));

  t = double(Est(:));

  sir = 10*log10( (s'*s) / ((s-t)'*(s-t)) );


