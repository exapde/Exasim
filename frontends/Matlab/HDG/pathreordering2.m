function [fpath, lpath, fintf, lintf] = pathreordering2(epath, e2f)

[npaths, nep] = size(epath);
nfe = size(e2f,1);

fpath = zeros(2,npaths,nep);
lpath = zeros(2,npaths,nep);
fintf = zeros(nfe-2,npaths,nep);
lintf = zeros(nfe-2,npaths,nep);

for i = 1:npaths  
  [fpath(:,i,:), lpath(:,i,:), fintf(:,i,:), lintf(:,i,:)] = pathreorder(epath(i,:), e2f);
end

ne = npaths*nep;
fpath = reshape(fpath, [2 ne]);
lpath = reshape(lpath, [2 ne]);
fintf = reshape(fintf, [nfe-2 ne]);
lintf = reshape(lintf, [nfe-2 ne]);

end

function [fpath, lpath, fintf, lintf] = pathreorder(epath, e2f)

nfe = size(e2f,1);
nep = length(epath);

fpath = zeros(2, nep);
lpath = zeros(2, nep);
for i = 1:(nep-1)    
  e1 = epath(i);   % left element 
  e2 = epath(i+1); % right element 
  f1 = e2f(:,e1);  % faces on left element
  f2 = e2f(:,e2);  % faces on right element
  match = 0;
  for j = 1:nfe
    for k = 1:nfe
      if f1(j) == f2(k)
        match = 1;
        break;
      end
    end
    if match == 1
      break;
    end
  end

  % (j, f1(j)) is the face on left element
  % (k, f2(k)) is the face on right element
  
  if (i==1)
    jk = j;
    if nfe == 8
      if mod(jk, 2) == 1 % jk is odd 
          m = jk + 1;
      else  % jk is even 
          m = jk - 1;
      end    
    elseif nfe == 4
      if jk == 1
          m = 3;
      elseif jk==3
          m = 1;
      elseif jk == 2
          m = 4;
      elseif jk==4
          m = 2;    
      end      
    elseif nfe == 3
      if jk == 1
          m = 3;
      elseif jk ==2
          m = 1;
      elseif jk == 3
          m = 2;    
      end        
    end        
    fpath(1,i) = f1(m);   % start face
    fpath(2,i) = f1(j);   % left face           
    fpath(1,i+1) = f2(k); % right face      
    lpath(1,i) = m;
    lpath(2,i) = j;      
    lpath(1,i+1) = k;    
  end
    
  if (i==nep-1)
    jk = k;
    if nfe == 8
      if mod(jk, 2) == 1 % jk is odd 
          m = jk + 1;
      else  % jk is even 
          m = jk - 1;
      end    
    elseif nfe == 4
      if jk == 1
          m = 3;
      elseif jk==3
          m = 1;
      elseif jk == 2
          m = 4;
      elseif jk==4
          m = 2;    
      end      
    elseif nfe == 3
      if jk == 1
          m = 3;
      elseif jk ==2
          m = 1;
      elseif jk == 3
          m = 2;    
      end        
    end    
    fpath(1,nep) = f2(k);   % right face
    fpath(2,nep) = f2(m);   % end face
    fpath(2,nep-1) = f1(j); % left face
    lpath(1,nep) = k;
    lpath(2,nep) = m;  
    lpath(2,nep-1) = j;    
  end

  if (i ~= 1) && (i ~= nep-1)
    fpath(2,i) = f1(j);     % left face
    fpath(1,i+1) = f2(k);   % right face
    lpath(2,i) = j;
    lpath(1,i+1) = k;
  end
end

fintf = zeros(nfe-2, nep);
lintf = zeros(nfe-2, nep);
for i = 1:(nep)    
  j = lpath(1,i);
  m = lpath(2,i);
  fi = e2f(:,epath(i));
  n = 0; 
  for l = 1:nfe
    if (l ~= m) && (l ~= j)
      n = n + 1;
      fintf(n,i) = fi(l);          
      lintf(n,i) = l;      
    end
  end      
end

end




