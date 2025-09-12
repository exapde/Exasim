function [pc, hc, dgnodes, XA, YA, XB, YB] = mshOrion( order, Nx1, Nx2, Nx3, Nbl)
% mshOrion - generates a FEM mesh for the axisymmetric Orion capsule geometry
%
% Input parameters:
%       order - order of FEM approximation (element nodes are equally spaced)
%       Nx1 - number of elements along streamwise direction for block1 
%       Nx2 - number of elements along streamwise direction for block2 
%       Nx3 - number of elements along streamwise direction for block3 
%       Nbl - number of elements in the normal direction away from the body 
%
% Output parameters:
%       p(NP, 3) - x,y,z coordinates of the nodes
%       h(NE, (order+1)^2) - hex element connectivitites
%       dgnodes((order+1)^2, 2, NE)
%

    %Nx1 = 40;
    %Nx2 = 50;
    %Nx3 = 20;
    %Nbl = 50;    % Normal Direction - Same for all block

    %order = 2;

    lambda_bl = 6;  % Boundary layer stretching

    XA = -0.1;   % Outer boudnary location
    YA =  4.0;
    
    XB = 7.0;
    YB = 5.0;
    
    snap = 0.5e-10;

    Nx_1 = Nx1*order + 1;
    Nx_2 = Nx2*order + 1;
    Nx_3 = Nx3*order + 1;
    Ny   = Nbl*order + 1;
    
    x0 = 0.0000;
    y0 = 0.0000;
    
    x1 = 0.4812;
    y1 = 2.3615;
    
    x2 = 0.7126;
    y2 = 2.5146;
    
    x3 = 0.8477;
    y3 = 2.4752;
    
    x4 = 3.2785;
    y4 = 0.9266;
    
    x5 = 3.3020;
    y5 = 0.8838;
    
    x6 = 3.3020;
    y6 = 0.0000;
    
    Rn = 6.0350;
    Rs = 0.2515;
    Rb = 2.5146;
    Ras = 0.0508;
    
    % Arc 0
    ThN = 23.0353*pi/180;
    th = 0:ThN/40:ThN;
    xc0 = Rn*(1 - cos(th));
    yc0 = Rn*sin(th);
    
    % Arc 1
    xc1 = Rs*cos(th);
    yc1 = Rs*sin(th);
    
    Thi = pi - ThN;
    c1_x =  xc0(end) - Rs*cos(Thi);
    c1_y =  yc0(end) - Rs*sin(Thi);
    Thf = atan2(YA - c1_y, XA - c1_x);
    
    th = Thi:(Thf-Thi)/10:Thf;
    xc1_1 = Rs*cos(th);
    yc1_1 = Rs*sin(th);
    
    xc1_1 = xc1_1 + c1_x;
    yc1_1 = yc1_1 + c1_y;
    
    l = sqrt((XA-c1_x)^2 + (YA-c1_y)^2);
    XC = xc1_1(end);
    YC = yc1_1(end);
    
    Thi = Thf;
    ThA = 32.5*pi/180;
    Thf = pi/2 - ThA;
    
    th = Thi:(Thf-Thi)/10:Thf;
    xc1_2 = Rs*cos(th);
    yc1_2 = Rs*sin(th);
    
    xc1_2 = xc1_2 + c1_x;
    yc1_2 = yc1_2 + c1_y;
    
    xl1 = xc1_2(end):(x4-xc1_2(end))/20:x4;
    yl1 = yc1_2(end):(y4-yc1_2(end))/20:y4;
    
    c2_x = x5 - Ras;
    c2_y = y5;
    
    l = sqrt((XB-c2_x)^2 + (YB-c2_y)^2);
    xpc = c2_x + (XB - c2_x)*Ras/l;
    ypc = c2_y + (YB - c2_y)*Ras/l;
    
    xb3 = xpc:(XB-xpc)/20:XB;
    yb3 = ypc:(YB-ypc)/20:YB;
    
    Thi = atan2(y4-c2_y, x4-c2_x);
    Thf = atan2(YB - c2_y, XB - c2_x);
    th = Thi:(Thf-Thi)/5:Thf;
    
    xc2_1 = c2_x + Ras*cos(th);
    yc2_1 = c2_y + Ras*sin(th);
    
    Thi = Thf;
    Thf = 0;
    th = Thi:(Thf-Thi)/10:Thf;
    
    xc2_2 = c2_x + Ras*cos(th);
    yc2_2 = c2_y + Ras*sin(th);
    
    yl2 = yc2_2(end):(y6-yc2_2(end))/20:y6;
    xl2 = xc2_2(end) + yl2*0;
    
    %plot(xc0,yc0,'k-',xc1_1,yc1_1,'r-',xc1_2,yc1_2,'g-',xl1,yl1,'b-'); hold on;
    %plot(xc2_1,yc2_1,'k-',xc2_2,yc2_2,'g-', xl2,yl2,'b-'); 
    %axis equal;
    
    c0_x = Rn;
    c0_y = 0;
    
    % plot(XA, YA, 'r*', c0_x, c0_y, 'r*', c1_x, c1_y, 'r*');
    
    ang0 = atan2(YA,XA-c0_x);
    R = sqrt(YA^2 + (XA-c0_x)^2);
    th = pi:(ang0-pi)/20:ang0;
    
    xb0 = c0_x + R*cos(th);
    yb0 = c0_y + R*sin(th);
    
    xb1 = XC:(XA-XC)/20:XA;
    yb1 = YC:(YA-YC)/20:YA;
    
    xb2 = XA:(XB-XA)/20:XB;
    yb2 = YA:(YB-YA)/20:YB;
    
    l = sqrt((XB-c2_x)^2 + (YB-c2_y)^2);
    xpc = c2_x + (XB - c2_x)*Ras/l;
    ypc = c2_y + (YB - c2_y)*Ras/l;
    
    xb3 = xpc:(XB-xpc)/20:XB;
    yb3 = ypc:(YB-ypc)/20:YB;
    
    yb4 = YB:-YB/20:0;
    xb4 = XB + yb4*0;
    
    % plot(xb0, yb0, 'k', xb1, yb1, 'k', xb2, yb2, 'k', xb3, yb3, 'k', xb4, yb4, 'k', 'LineWidth', 2);
    
    
    % Block 1
    [xl, yl] = distribute_points_exponential([xc0, xc1_1(1,2:end)], [yc0, yc1_1(1,2:end)], Nx_1, 2, 3);
    [xu, yu] = distribute_points_exponential(xb0, yb0, Nx_1, 0.01, 0.01);
    
    b1 = zeros(Nx_1, Ny, 2);
    for i = 1:Nx_1
        %plot([xl(1,i); xu(1,i)], [yl(1,i); yu(1,i)], 'k-');
        [b1(i,:,1), b1(i,:,2)] = distribute_points_exponential([xl(1,i); xu(1,i)], [yl(1,i); yu(1,i)], Ny, lambda_bl, 0.01);
    end
    
    %for j = 1:Ny
    %    plot(b1(:,j,1), b1(:,j,2), 'k-');
    %end
    
    % Block 2
    [xl, yl] = distribute_points_exponential([xc1_2, xl1(1,2:end), xc2_1(1,2:end)], [yc1_2, yl1(1,2:end), yc2_1(1,2:end)], Nx_2, 4, 4);
    [xu, yu] = distribute_points_exponential(xb2, yb2, Nx_2, 0.01, 0.01);
    
    b2 = zeros(Nx_2, Ny, 2);
    for i = 1:Nx_2
        %plot([xl(1,i); xu(1,i)], [yl(1,i); yu(1,i)], 'k-');
        [b2(i,:,1), b2(i,:,2)] = distribute_points_exponential([xl(1,i); xu(1,i)], [yl(1,i); yu(1,i)], Ny, lambda_bl, 0.01);
    end
    
    %for j = 1:Ny
    %    plot(b2(:,j,1), b2(:,j,2), 'k-');
    %end
    
    % Block 3
    [xl, yl] = distribute_points_exponential([xc2_2, xl2(1,2:end)], [yc2_2, yl2(1,2:end)], Nx_3, 3, 1);
    [xu, yu] = distribute_points_exponential(xb4, yb4, Nx_3, 0.01, 0.01);
    
    b3 = zeros(Nx_3, Ny, 2);
    for i = 1:Nx_3
        %plot([xl(1,i); xu(1,i)], [yl(1,i); yu(1,i)], 'k-');
        [b3(i,:,1), b3(i,:,2)] = distribute_points_exponential([xl(1,i); xu(1,i)], [yl(1,i); yu(1,i)], Ny, lambda_bl, 0.01);
    end

    %for j = 1:Ny
    %    plot(b3(:,j,1), b3(:,j,2), 'k-');
    %end

    [p1, h1] = extract_elements(b1, order);
    [p2, h2] = extract_elements(b2, order);
    h2 = h2 + size(p1,1);
    [p3, h3] = extract_elements(b3, order);
    h3 = h3 + size(p1,1) + size(p2,1);

    hg = cat(1, h1, h2, h3);
    pg = cat(1, p1, p2, p3);

    [foo,ix,jx] = unique(round(pg/snap)*snap,'rows');
    p = pg(ix,:);
    h = jx(hg);

    order1 = order + 1;
    aux1 = [ 1, order1];
    localnodes = cat(2, aux1, aux1 + order1*order);

    npt = size(p,1);
    lp = zeros(npt,1);
    for ie = 1:size(h,1)
       lp(h(ie,localnodes(:))) = 1;
    end
    
    acc = 0;
    for ip = 1:size(lp,1)
       acc = acc + lp(ip);
       lp(ip) = lp(ip)*acc;
    end
    
    pc = p(lp>0,:);
    hc = lp(h(:,localnodes([1,2,4,3])));

    dgnodes = zeros(size(h,2), 2, size(h,1));
    for ie = 1:size(h,1)
        dgnodes(:,:,ie) = p(h(ie,:),:);
    end

    x = dgnodes(:,1,:);
    y = dgnodes(:,2,:);
    plot(x(:), y(:), 'k.'); axis equal; hold on;

    for ie = 1:size(hc,1)
       plot(pc([hc(ie,:),hc(ie,1)],1), pc([hc(ie,:),hc(ie,1)],2),  'k-', 'LineWidth', 2);
    end


end

function [p, h] = extract_elements(mesh, order)
    [mx, my, d] = size(mesh);
    nx = (mx-1)/order;
    ny = (my-1)/order;

    node = zeros(mx, my);
    node(:) = 1:(mx*my);

    p = reshape(mesh, mx*my, 2);

    h = zeros(nx*ny,(order+1)^2);

    order1 = order + 1;

    ie = 0;
    for j = 1:ny
         for i = 1:nx
             ie = ie + 1;
             for jl = 1:order1
                    jn = (j-1)*order + jl;
                    for il = 1:order1
                        in = (i-1)*order +il;
                        h(ie, il + (jl-1)*order1) = node(in, jn);
                    end
             end
        end
    end
end