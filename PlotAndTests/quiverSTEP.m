x=mesh.xgrid*1e4; y=mesh.ygrid*1e4;

[X,Y]=meshgrid(x,y);

X0=x;
Y0=y;
dx=diff(X0);
sogx=.1;
su=0;
sX=[];
for k=1:length(dx)
    su=su+dx(k);
    if su>sogx
        su=0;
        sX=[sX k];
    end
end
sX0=sX;

dx=diff(Y0);
sogx=.1;
su=0;
sY=[];
for k=1:length(dx)
    su=su+dx(k);
    if su>sogx
        su=0;
        sY=[sY k];
    end
end
sY0=sY;

X0=x;
Y0=y;
dx=diff(X0);
sogx=SZ;
su=0;
sX=[];
for k=1:length(dx)
    su=su+dx(k);
    if su>sogx
        su=0;
        sX=[sX k];
    end
end
sX1=sX;

dx=diff(Y0);
sogx=SZ;
su=0;
sY=[];
for k=1:length(dx)
    su=su+dx(k);
    if su>sogx
        su=0;
        sY=[sY k];
    end
end
sY1=sY;
  
SX=X(sY0,sX0);
SY=Y(sY0,sX0);