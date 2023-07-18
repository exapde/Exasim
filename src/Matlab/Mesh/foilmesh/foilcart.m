function [X,Y,X3,Y3] = foilcart(xfl,yfl,xfu,yfu,porder,nl,nu,nr1,nr2,nw,nt,a,b1,b2,c,R1,R2,L)

[X,Y] = foilcart1(xfl,yfl,xfu,yfu,nl,nu,nr1,porder,a,b1,R1);

[X,Y] = foilcart2(X,Y,nr2,porder,b2,R1,R2);

[X3, Y3] = foilcart3(X,Y,nw,nt,porder,c,L);



