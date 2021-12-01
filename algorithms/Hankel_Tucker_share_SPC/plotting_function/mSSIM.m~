function [val vals ssimap] = mSSIM(X,Y);

  for t = 1:size(X,3)
    x = X(:,:,t);
    y = Y(:,:,t);
    [vals(t) ssimap{t}] = ssim_sub(x,y);
  end
  val = mean(vals);
