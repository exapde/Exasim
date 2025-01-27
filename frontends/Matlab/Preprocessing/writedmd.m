function writedmd(dmd, app, isd)

if app.modelnumber==0
    strn = "";
else
    strn = num2str(app.modelnumber);
end
filename = app.buildpath + "/datain" + strn + "/";

for i = 1:length(isd)
  d = isd(i);  
        
  fileID = fopen(filename + "mesh" + string(d) + ".bin",'r');
  tmp = fread(fileID,'double');
  fclose(fileID);
  
  if tmp(41)==0
    tmp(41) = length(dmd{d}.faceperm(:));
    tmp = [tmp(:); dmd{d}.faceperm(:)-1];
  else
    error("nsize(41) is not zero");
  end
  if tmp(42)==0
    tmp(42) = length(dmd{d}.nbintf(:));
    tmp = [tmp(:); dmd{d}.nbintf(:)-1];
  else
    error("nsize(42) is not zero");
  end  
  if tmp(43)==0
    tmp(43) = length(dmd{d}.facesend(:,2));
    tmp = [tmp(:); dmd{d}.facesend(:,2)-1];
  else
    error("nsize(43) is not zero");
  end
  if tmp(44)==0
    tmp(44) = length(dmd{d}.facesendpts(:));
    tmp = [tmp(:); dmd{d}.facesendpts(:)];
  else
    error("nsize(44) is not zero");
  end
  if tmp(45)==0
    tmp(45) = length(dmd{d}.facerecv(:,2));
    tmp = [tmp(:); dmd{d}.facerecv(:,2)-1];
  else
    error("nsize(45) is not zero");
  end
  if tmp(46)==0
    tmp(46) = length(dmd{d}.facerecvpts(:));
    tmp = [tmp(:); dmd{d}.facerecvpts(:)];
  else
    error("nsize(46) is not zero");
  end
  
  fileID2 = fopen(filename + "mesh" + string(d) + ".bin",'w');
  fwrite(fileID2,tmp(:),'double','native');
  fclose(fileID2);  
end
