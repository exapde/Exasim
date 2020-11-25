writerObj = VideoWriter('RBm.avi');
open(writerObj);
for K = 0 : 199
  if K < 10
  filename = sprintf('movie/rbm.000%1d.png', K);
  elseif K < 100
  filename = sprintf('movie/rbm.00%2d.png', K);
  else
  filename = sprintf('movie/rbm.0%3d.png', K);
  end
  thisimage = imread(filename);
  writeVideo(writerObj, thisimage);
end
close(writerObj);