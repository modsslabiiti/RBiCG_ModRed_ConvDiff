function [a,rhs] = disnsy(nx, ny, xl, xu, yl, yu, hx, hy)
x2y = hx/(2*hy);
hx2 = hx/2;
y2x = hy/(2*hx);
hy2 = hy/2;
yy = yl;

% discretisation see tr-11, page -3-
j = 1;
for iy = 1:ny
    xx = xl;
    for ix = 1:nx
        a1 = -x2y*(disdcc(xx-hx2,yy-hy2,xl,xu,yl,yu,@ec) + ...
                disdcc(xx+hx2,yy-hy2,xl,xu,yl,yu,@ec));
        aj1 = -disbcc(xx,yy,iy,ny,@dery)*hx2;
        a2 = -y2x*(disdcc(xx-hx2,yy-hy2,xl,xu,yl,yu,@dc) + ...
                disdcc(xx-hx2,yy+hy2,xl,xu,yl,yu,@dc));
        aj2 = -disbcc(xx,yy,ix,nx,@derx)*hy2;
        a4 = -y2x*(disdcc(xx+hx2,yy+hy2,xl,xu,yl,yu,@dc) + ...
                disdcc(xx+hx2,yy-hy2,xl,xu,yl,yu,@dc));
        a5 = -x2y*(disdcc(xx-hx2,yy+hy2,xl,xu,yl,yu,@ec) + ...
                disdcc(xx+hx2,yy+hy2,xl,xu,yl,yu,@ec));
        a(j,3) = -(a1+a2+a4+a5) + hx*hy/4* ...
              (disdcc(xx-hx2,yy-hy2,xl,xu,yl,yu,@cc)+ ...
               disdcc(xx-hx2,yy+hy2,xl,xu,yl,yu,@cc)+ ...
               disdcc(xx+hx2,yy-hy2,xl,xu,yl,yu,@cc)+ ...
               disdcc(xx+hx2,yy+hy2,xl,xu,yl,yu,@cc));
        rhs(j) = hx*hy/4*(disdcc(xx-hx2,yy-hy2,xl,xu,yl,yu,@fc)+ ...
              disdcc(xx-hx2,yy+hy2,xl,xu,yl,yu,@fc)+ ... 
              disdcc(xx+hx2,yy-hy2,xl,xu,yl,yu,@fc)+ ...
              disdcc(xx+hx2,yy+hy2,xl,xu,yl,yu,@fc));
        a(j,1) = a1 + aj1;
        a(j,2) = a2 + aj2;
        a(j,4) = a4 - aj2;
        a(j,5) = a5 - aj1;
        j = j + 1;
        xx = xx + hx;
    end
    yy = yy + hy;
end

% dirichlet boundary condition, south boundary
ij = nx + 1;
for i = 1:nx
    [alp, bet, gam] = mbnds;
    rhs(ij) = rhs(ij) - a(ij,1)*gam/alp;
    ij = ij + 1;
end

% dirichlet condition on west boundary
ij = 2;
for i=1:ny
    [alp, bet, gam] = mbndw;
    rhs(ij) = rhs(ij) - a(ij,2)*gam/alp;
    ij = ij + nx;
end

% dirichlet boundary condition on north boundary
ij = nx*(ny-2) + 1;
for i=1:nx
    [alp, bet, gam] = mbndn;
    rhs(ij) = rhs(ij) - a(ij,5)*gam/alp;
    ij = ij + 1;
end

% dirichlet condition on east boundary
ij = nx - 1;
for i=1:ny
    [alp, bet, gam] = mbnde;
    rhs(ij) = rhs(ij) - a(ij,4)*gam/alp;
    ij = ij + nx;
end

% the linear system will be comprimised for dirichlet boundary condition
ny = ny - 1;
n = nx*ny;
for i=1:n
    a(i,1) = a(i+nx,1);
    a(i,2) = a(i+nx,2);
    a(i,3) = a(i+nx,3);
    a(i,4) = a(i+nx,4);
    a(i,5) = a(i+nx,5);
    rhs(i) = rhs(i+nx);
end

is = 1;
nx = nx - 1;
ij = 1;
for i=1:ny
    for j=1:nx
        a(ij,1) = a(ij+is,1);
        a(ij,2) = a(ij+is,2);
        a(ij,3) = a(ij+is,3);
        a(ij,4) = a(ij+is,4);
        a(ij,5) = a(ij+is,5);
        rhs(ij) = rhs(ij+is);
        ij = ij + 1;
    end
    is = is + 1;
end

ny = ny - 1;
nx = nx - 1;
n = nx*ny;
is = 1;
ij = nx + 1;
for i=2:ny
    for j=1:nx
        a(ij,1) = a(ij+is,1);
        a(ij,2) = a(ij+is,2);
        a(ij,3) = a(ij+is,3);
        a(ij,4) = a(ij+is,4);
        a(ij,5) = a(ij+is,5);
        rhs(ij) = rhs(ij+is);
        ij = ij + 1;
    end
    is = is + 1;
end

ij = n - nx + 1;
for i=1:nx
    a(i,1) = 0;
    a(ij,5) = 0;
    ij = ij + 1;
end
i1 = 1;
i2 = nx;
for i=1:ny
    a(i1,2) = 0;
    a(i2,4) = 0;
    i1 = i1 + nx;
    i2 = i2 + nx;
end

function val = disdcc(x, y, xl, xu, yl, yu, xc)
if (x < xl || x > xu || y < yl || y > yu)
    val = 0;
else
    val = xc(x,y);
end

function val  = disbcc(x, y, i, ni, der)
if (i == 1 || i == ni)
  val = 0;
else
  val = der(x,y);
end
