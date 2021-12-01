function [val ssimap] = SSIM(X,Y);

  if ndims(X) == 3 & size(X,3) == 3
    x = 0.299*X(:,:,1) + 0.587*X(:,:,2) + 0.114*X(:,:,3);
    y = 0.299*Y(:,:,1) + 0.587*Y(:,:,2) + 0.114*Y(:,:,3);
  end
  [val ssimap] = ssim(x,y);
