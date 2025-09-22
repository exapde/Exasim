function p=mapp(p,x)
% map points in a reference element to points in physical element
% p : points in reference element
% x : coordinates of physical element

nx = size(x,1);
nd = size(x,2);

if nd==1
     r=[0,1];
     
     C=x(:,1)'/[1,1; r];
     
     p=[0*p+1,p]*C';
elseif (nd==2) && (nx==4) % 2-D bilinear map    
     r=[0,1,0,1];
     s=[0,0,1,1];
     
     C=x(:,1)'/[1,1,1,1; r; s; r.*s];
     D=x(:,2)'/[1,1,1,1; r; s; r.*s];

     px=p(:,1); py=p(:,2);
     p=[0*px+1, px, py, px.*py]*[C;D]';
elseif (nd==3) && (nx==8) % 3-D bilinear map    
     r=[0,1,0,1,0,1,0,1];
     s=[0,0,1,1,0,0,1,1];
     t=[0,0,0,0,1,1,1,1];
     
     C=x(:,1)'/[1,1,1,1,1,1,1,1; r; s; t; r.*s; r.*t; s.*t; r.*s.*t];
     D=x(:,2)'/[1,1,1,1,1,1,1,1; r; s; t; r.*s; r.*t; s.*t; r.*s.*t];
     E=x(:,3)'/[1,1,1,1,1,1,1,1; r; s; t; r.*s; r.*t; s.*t; r.*s.*t];

     px=p(:,1); py=p(:,2); pz=p(:,3);
     p=[0*px+1, px, py, pz, px.*py, px.*pz, py.*pz, px.*py.*pz]*[C;D;E]';
elseif (nd==2) && (nx==3) % 2-D affine map    
     r=[0,1,0];
     s=[0,0,1];
     
     C=x(:,1)'/[1,1,1; r; s];
     D=x(:,2)'/[1,1,1; r; s];

     px=p(:,1); py=p(:,2);
     p=[0*px+1, px, py]*[C;D]'; 
elseif (nd==3) && (nx==4) % 3-D affine map    
     r=[0,1,0,0];
     s=[0,0,1,0];
     t=[0,0,0,1];
     
     C=x(:,1)'/[1,1,1,1; r; s; t];
     D=x(:,2)'/[1,1,1,1; r; s; t];
     E=x(:,3)'/[1,1,1,1; r; s; t];

     px=p(:,1); py=p(:,2); pz=p(:,3);
     p=[0*px+1, px, py, pz]*[C;D;E]';
end
