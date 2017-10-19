%O. Hellwig, 22.07.2010
%last update 04.05.2012

%close all
clear all

% This m-file is for making movies of wave propagation.

%%%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%%
% Here the basic name of the binary snapshot files must be given:
% (The default extension is *.bin. The increasing number of the snapshot
% is added to the basic name, e.g. if ten snapshots were computed
% and the basic name is snap, then the filenames for the snapshots
% are snap1.bin, snap2.bin, ..., snap10.bin.

%model name
name = 'Test';

%parameter input file
%file = '../model/fd3dp.inp';
file = '../model/fd3dps.inp';

%path to par-directoryl
%par = '/data/hellwig/software/fd3d_struct/par/';
par = '~/software/fd3d_struct5/par/';

%plane to display
%p=1 ... y-z-plane at x=x0 in m
%p=2 ... x-z-plane at y=y0 in m
%p=3 ... x-y-plane at z=z0 in m
p = 2;

x0 = 110;
y0 = 20;
z0 = 30;
%x0 = 1000;
%y0 = 1000;
%z0 = 1000;

%range to plot wavefield
%p=1 ... [ymin ymax zmin zmax]
%p=2 ... [xmin xmax zmin zmax]
%p=3 ... [xmin xmax ymin ymax]
%ar = [0 40  0 60];
ar = [0 220 0 60];
%ar = [0 220 0 40];
%ar = [0 2000  0 2000];
%ar = [0 2000 0 2000];
%ar = [0 2000  0 2000];

%display source and receiver positions (yes...k=1, no...else)
k = 1;
%receiver numbers to display
%recnum = [1:60];
recnum = [1:54];
%source numbers to display
srcnum = [1];

%l=1 ... plot vx, vy, vz
%l=2 ... plot div, curlx, curly, curlz
%l=3 ... plot p
%l=4 ... plot txx, txy, txz, tyy, tyz, tzz
%l=5 ... plot vx, vy, vz, p
l = 1;

%model as background (yes...m=1, no...else)
m = 1;

%color axis
caxis_value1 = 1.0e-8;
caxis_value2 = 1.0e-8;

%save snapshots as
%n=1 ... eps
%n=2 ... jpg
%n=3 ... ppm
%n=4 ... png
%n=5 ... jpg and png
n = 4;
%save snapshots in
snapfile = '../pics/snap';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read parameter input file
disp(['reading parameter file ' file]);
disp(' ');

refmod = [0.0 0.0 0.0];
refsrc = [0.0 0.0 0.0];
refrec = [0.0 0.0 0.0];
FFID   = [];

fid = fopen(file,'r','ieee-le');
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break;
    end;
    
    if sum(findstr(tline,'(NPROCX)'))
        NPROCX = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(NPROCY)'))
        NPROCY = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(NPROCZ)'))
        NPROCZ = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(NX)'))
        NX = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(NY)'))
        NY = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(NZ)'))
        NZ = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(DX_FILE)'))
        DXFILE = tline(findstr(tline,'=')+1:end);
    end;
    if sum(findstr(tline,'(DY_FILE)'))
        DYFILE = tline(findstr(tline,'=')+1:end);
    end;
    if sum(findstr(tline,'(DZ_FILE)'))
        DZFILE = tline(findstr(tline,'=')+1:end);
    end;
    if sum(findstr(tline,'(REFMOD)'))
        refmod = eval(['[',tline(findstr(tline,'=')+1:end),']']);
    end;    
    if sum(findstr(tline,'(TIME)'))
        tmax = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(DT)'))
        DT = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(REFSRC)'))
        refsrc = eval(['[',tline(findstr(tline,'=')+1:end),']']);
    end;
    
    if sum(findstr(tline,'(MODEL_FILE)'))
        MFILE = tline(findstr(tline,'=')+1:end);
    end;
    if sum(findstr(tline,'(MODEL_IDX)'))
        MIDX = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(MODEL_IDY)'))
        MIDY = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(MODEL_IDZ)'))
        MIDZ = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(READREC)'))
        a = findstr(tline,'=');
        a = a(2);
        READREC = eval(tline(a+1:end));
    end;
    if sum(findstr(tline,'(REFREC)'))
        refrec = eval(['[',tline(findstr(tline,'=')+1:end),']']);
    end;    
    if sum(findstr(tline,'REC_FILE'))
        RECFILE = tline(findstr(tline,'=')+1:end);
    end;
    if sum(findstr(tline,'FFID'))
        FFID = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'SOURCE_FILE'))
        SRCFILE = tline(findstr(tline,'=')+1:end);
    end;
    if sum(findstr(tline,'(SNAP)'))
        snap = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(TSNAP1)'))
        TSNAP1 = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(TSNAP2)'))
        TSNAP2 = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(TSNAPINC)'))
        TSNAPINC = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(XSNAPMIN)'))
        XSNAPMIN = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(XSNAPMAX)'))
        XSNAPMAX = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(IDX)'))
        IDX = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(YSNAPMIN)'))
        YSNAPMIN = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(YSNAPMAX)'))
        YSNAPMAX = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(IDY)'))
        IDY = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(ZSNAPMIN)'))
        ZSNAPMIN = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(ZSNAPMAX)'))
        ZSNAPMAX = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(IDZ)'))
        IDZ = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(SNAP_FORMAT)'))
        format = eval(tline(findstr(tline,'=')+1:end));
    end;
    if sum(findstr(tline,'(SNAP_FILE)'))
        SNAPFILE = tline(findstr(tline,'=')+1:end);
    end;
end;
fclose(fid);
disp(' ');


if READREC~=1 %no sources and receivers displayed
    k=0;
end;
DXFILE = DXFILE(find(DXFILE~=' '));
if DXFILE(1:2) == './'
    DXFILE = DXFILE(3:end);
end;
DYFILE = DYFILE(find(DYFILE~=' '));
if DYFILE(1:2) == './'
    DYFILE = DYFILE(3:end);
end;
DZFILE = DZFILE(find(DZFILE~=' '));
if DZFILE(1:2) == './'
    DZFILE = DZFILE(3:end);
end;
MFILE = MFILE(find(MFILE~=' '));
if MFILE(1:2) == './'
    MFILE = MFILE(3:end);
end;
SNAPFILE = SNAPFILE(find(SNAPFILE~=' '));
if SNAPFILE(1:2) == './'
    SNAPFILE = SNAPFILE(3:end);
end;
if (~isempty(FFID))
    SNAPFILE = [SNAPFILE,num2str(FFID,'%.4i')];
end;
RECFILE = RECFILE(find(RECFILE~=' '));
if RECFILE(1:2) == './'
    RECFILE = RECFILE(3:end);
end;
SRCFILE = SRCFILE(find(SRCFILE~=' '));
if SRCFILE(1:2) == './'
    SRCFILE = SRCFILE(3:end);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load distances between grid points
dx = load([par,DXFILE]);
dy = load([par,DYFILE]);
dz = load([par,DZFILE]);

DX = [];
for ii = 1:size(dx,1)
    DX = [DX;repmat(dx(ii,2),[dx(ii,1),1])];
end;
DY = [];
for ii = 1:size(dy,1)
    DY = [DY;repmat(dy(ii,2),[dy(ii,1),1])];
end;
DZ = [];
for ii = 1:size(dz,1)
    DZ = [DZ;repmat(dz(ii,2),[dz(ii,1),1])];
end;

%initialize x-, y- and z-vectors
x = zeros(NX+1,1);
y = zeros(NY+1,1);
z = zeros(NZ+1,1);

x(1) = 0;
x(2:min(length(DX),NX)+1)    = DX(1:min(length(DX),NX));
x(min(length(DX),NX)+2:NX+1) = DX(end);
x = 0.5*(x(1:NX)+x(2:NX+1));

y(1) = 0;
y(2:min(length(DY),NY)+1)    = DY(1:min(length(DY),NY));
y(min(length(DY),NY)+2:NY+1) = DY(end);
y = 0.5*(y(1:NY)+y(2:NY+1));

z(1) = 0;
z(2:min(length(DZ),NZ)+1)    = DZ(1:min(length(DZ),NZ));
z(min(length(DZ),NZ)+2:NZ+1) = DZ(end);
z = 0.5*(z(1:NZ)+z(2:NZ+1));

%compute x-, y- and z-vectors
x = cumsum(x) + refmod(1);
y = cumsum(y) + refmod(2);
z = cumsum(z) + refmod(3);

%compute x-, y- and z-vectors (model)
xm1  = x(1:MIDX:NX);
ym1  = y(1:MIDY:NY);
zm1  = z(1:MIDZ:NZ);

%size x-, y- and z-vectors (model)
nxm = length(xm1);
nym = length(ym1);
nzm = length(zm1);

%compute x-, y- and z-vectors (snapshot)
xs1  = x(1:IDX:NX);
ys1  = y(1:IDY:NY);
zs1  = z(1:IDZ:NZ);

%compute x-, y- and z-vectors (snapshots)
if XSNAPMIN == []
    XSNAPMIN = 0;
end;
if XSNAPMAX == []
    XSNAPMAX = x(end);
end;
if YSNAPMIN == []
    YSNAPMIN = 0;
end;
if YSNAPMAX == []
    YSNAPMAX = y(end);
end;
if ZSNAPMIN == []
    ZSNAPMIN = 0;
end;
if ZSNAPMAX == []
    ZSNAPMAX = z(end);
end;

xi  = find((xs1>=XSNAPMIN)&(xs1<=XSNAPMAX));
yi  = find((ys1>=YSNAPMIN)&(ys1<=YSNAPMAX));
zi  = find((zs1>=ZSNAPMIN)&(zs1<=ZSNAPMAX));

xs1  = xs1(xi);
ys1  = ys1(yi);
zs1  = zs1(zi);

%size x-, y- and z-vectors (snapshots)
nxs = length(xs1);
nys = length(ys1);
nzs = length(zs1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%range of the axes for display
dxm = min(diff(xs1));
dym = min(diff(ys1));
dzm = min(diff(zs1));

x2 = [xs1(1):dxm:xs1(end)];
y2 = [ys1(1):dym:ys1(end)];
z2 = [zs1(1):dzm:zs1(end)];

%determine plane index
xmi = find(min(abs(x0-xm1)) == (abs(x0-xm1)));
xmi = xmi(1);
xsi = find(min(abs(x0-xs1)) == (abs(x0-xs1)));
xsi = xsi(1);

ymi = find(min(abs(y0-ym1)) == (abs(y0-ym1)));
ymi = ymi(1);
ysi = find(min(abs(y0-ys1)) == (abs(y0-ys1)));
ysi = ysi(1);

zmi = find(min(abs(z0-zm1)) == (abs(z0-zm1)));
zmi = zmi(1);
zsi = find(min(abs(z0-zs1)) == (abs(z0-zs1)));
zsi = zsi(1);

%determine plane position
x3 = xm1(xmi);
y3 = ym1(ymi);
z3 = zm1(zmi);

%mesh for interpolation
switch p
    case 1
        disp(['Snapshots at x0=',num2str(x3),' m']);
        [ym1,zm1] = meshgrid(ym1,zm1);
        [ys1,zs1] = meshgrid(ys1,zs1);
        [y2,z2]   = meshgrid(y2,z2);
    case 2
        disp(['Snapshots at y0=',num2str(y3),' m']);
        [xm1,zm1] = meshgrid(xm1,zm1);
        [xs1,zs1] = meshgrid(xs1,zs1);
        [x2,z2]   = meshgrid(x2,z2);
    case 3
        disp(['Snapshots at z0=',num2str(z3),' m']);
        [xm1,ym1] = meshgrid(xm1,ym1);
        [xs1,ys1] = meshgrid(xs1,ys1);
        [x2,y2]   = meshgrid(x2,y2);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load source and receiver positions
if k
    recs = load([par,RECFILE]);
    recs(:,1) = recs(:,1) + refrec(1);
    recs(:,2) = recs(:,2) + refrec(2);
    recs(:,3) = recs(:,3) + refrec(3);
    switch p
        case 1
            recs = recs(recnum,2:3);
        case 2
            recs = recs(recnum,1:2:3);
        case 3
            recs = recs(recnum,1:2);
    end;
    
    fid = fopen([par,SRCFILE],'r');
    nsrc = eval(fgetl(fid));
    src  = fscanf(fid,'%f',[8,nsrc]);
    fclose(fid);
    src = src';
    src(:,1) = src(:,1) + refsrc(1);
    src(:,2) = src(:,2) + refsrc(2);
    src(:,3) = src(:,3) + refsrc(3);
    switch p
        case 1
            src  = src(srcnum,2:3);
        case 2
            src  = src(srcnum,1:2:3);
        case 3
            src  = src(srcnum,1:2);
    end;
end;

%load model
if m
    %read model file
    fidm = fopen([par,MFILE,'_rho.bin'],'r','ieee-le');
    disp(['opening and reading file ',MFILE]);
    
    switch p
        case 1
            form   = [nzm,nym];
            skip   = 0;
            offset = 4*nzm*nym*(xmi-1);
            
            fseek(fidm,offset,-1);
            model = fread(fidm,form,'float32',skip);
        case 2
            form   = [nzm,1];
            skip   = 0;
            
            model = zeros(nzm,nxm);
            for ik=1:nxm
                offset = 4*nzm*(nym*(ik-1) + ymi-1);
                
                fseek(fidm,offset,-1);
                model(:,ik) = fread(fidm,form,'float32',skip);
            end;
        case 3
            form   = [nym,nxm];
            skip   = 4*(nzm-1);
            offset = 4*(zmi-1);
            
            fseek(fidm,offset,-1);
            model = fread(fidm,form,'float32',skip);
    end;
    fclose(fidm);
    
    %interpolate model to grid for display
    switch p
        case 1
            model = interp2(ym1,zm1,model,y2,z2);
        case 2
            model = interp2(xm1,zm1,model,x2,z2);
        case 3
            model = interp2(xm1,ym1,model,x2,y2);
    end;
    
    %scale model values for display to [-0.2 0.2]
    model = model-min(min(model));
    if max(max(model) ~= 0.0)
        model = 0.4*model/max(max(model))-0.2;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%open snapshot files
if (isempty(FFID))
    id = [];
else
    id = num2str(FFID,'%.4d');
end;

switch l
    case 1
        num = 3;
        sub1 = 2;
        sub2 = 2;
        fid1 = fopen([par,SNAPFILE,'_vx.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_vx.bin']);
        fid2 = fopen([par,SNAPFILE,'_vy.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_vy.bin']);
        fid3 = fopen([par,SNAPFILE,'_vz.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_vz.bin']);
    case 2
        num = 4;
        sub1 = 2;
        sub2 = 2;
        fid1 = fopen([par,SNAPFILE,'_div.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_div.bin']);
        fid2 = fopen([par,SNAPFILE,'_curlx.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_curlx.bin']);
        fid3 = fopen([par,SNAPFILE,'_curly.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_curly.bin']);
        fid4 = fopen([par,SNAPFILE,'_curlz.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_curlz.bin']);
    case 3
        num = 1;
        sub1 = 1;
        sub2 = 1;
        fid1 = fopen([par,SNAPFILE,'_p.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_p.bin']);
    case 4
        num = 6;
        sub1 = 3;
        sub2 = 2;
        fid1 = fopen([par,SNAPFILE,'_txx.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_txx.bin']);
        fid2 = fopen([par,SNAPFILE,'_tyy.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_tyy.bin']);
        fid3 = fopen([par,SNAPFILE,'_tzz.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_tzz.bin']);
        fid4 = fopen([par,SNAPFILE,'_txy.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_txy.bin']);
        fid5 = fopen([par,SNAPFILE,'_txz.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_txz.bin']);
        fid6 = fopen([par,SNAPFILE,'_tyz.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_tyz.bin']);
    case 5
        num = 4;
        sub1 = 2;
        sub2 = 2;
        fid1 = fopen([par,SNAPFILE,'_vx.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_vx.bin']);
        fid2 = fopen([par,SNAPFILE,'_vy.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_vy.bin']);
        fid3 = fopen([par,SNAPFILE,'_vz.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_vy.bin']);
        fid4 = fopen([par,SNAPFILE,'_p.bin'],'r','ieee-le');
        disp(['opening file ',SNAPFILE,'_p.bin']);
end;

%parameters for reading snapshots
switch p
    case 1
        form   = [nzs,nys];
        skip   = 0;
        label1 = 'y [m]';
        label2 = 'z [m]';
    case 2
        form   = [nzs,1];
        skip   = 0;
        label1 = 'x [m]';
        label2 = 'z [m]';
    case 3
        form   = [nys,nxs];
        skip   = 4*(nzs-1);
        label1 = 'x [m]';
        label2 = 'y [m]';
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting in window:
set(gcf,'Position',[50 50 1100 1100]);

%frame number
firstframe = 1;
lastframe  = (TSNAP2-TSNAP1)/TSNAPINC+1;
snapno     = 0;

for ii = firstframe:lastframe,
    hold off;
    disp(['loading snapshot no ',int2str(ii)]);
    
    %loading data:
    tsnap  = (ii-1)*TSNAPINC+TSNAP1;
   
    for ij = 1:num
        switch p
            case 1
                offset = 4*nzs*nys*(xsi - 1 + nxs*(ii-1));
                    
                fseek(eval(['fid',num2str(ij)]),offset,-1);
                v = fread(eval(['fid',num2str(ij)]),form,'float32',skip);
            case 2
                v = zeros(nzs,nxs);
                for ik=1:nxs
                    offset = 4*nzs*(nys*(nxs*(ii-1) + ik-1) + ysi-1);
                    
                    fseek(eval(['fid',num2str(ij)]),offset,-1);
                    v(:,ik) = fread(eval(['fid',num2str(ij)]),form,'float32',skip);
                end;
            case 3
                offset = 4*(nzs*nys*nxs*(ii-1) + zsi-1);
            
                fseek(eval(['fid',num2str(ij)]),offset,-1);
                v = fread(eval(['fid',num2str(ij)]),form,'float32',skip);
        end;
          
        %maximum value
        vm = max(max(abs(v)));
        disp([' Maximum amplitude of snapshots: ', num2str(vm)]);  
        
        %interpolation
        switch p
            case 1
                v = interp2(ys1,zs1,v,y2,z2);
            case 2
                v = interp2(xs1,zs1,v,x2,z2);
            case 3
                v = interp2(xs1,ys1,v,x2,y2);
        end;
                
        subplot(sub1,sub2,ij);
        hold on;
        
        switch ij
            case 1
                switch l
                    case 1
                        caxis_value = caxis_value1;
                        title([name]);
                        xlabel([label1,'     v_x']);
                    case 2
                        caxis_value = caxis_value1;
                        title([name]);
                        xlabel([label1,'     divergence']);
                    case 3
                        caxis_value = caxis_value1;
                        switch p
                            case 1
                                title([name,'         x=',sprintf('%1.2f',x3),' m         t = ',sprintf('%1.5f',tsnap),' s']);
                            case 2
                                title([name,'         y=',sprintf('%1.2f',y3),' m         t = ',sprintf('%1.5f',tsnap),' s']);
                            case 3
                                title([name,'         z=',sprintf('%1.2f',z3),' m         t = ',sprintf('%1.5f',tsnap),' s']);
                        end;
                        
                        xlabel([label1,'     pressure']);
                    case 4
                        caxis_value = caxis_value1;
                        title([name]);
                        xlabel([label1,'     \sigma_{xx}']);
                    case 5
                        caxis_value = caxis_value1;
                        title([name]);
                        xlabel([label1,'     v_x']);
                end;
            case 2    
                switch l
                    case 1
                        caxis_value = caxis_value1;
                        switch p
                            case 1
                                title(['x=',sprintf('%1.2f',x3),' m         t = ',sprintf('%1.5f',tsnap),' s']);
                            case 2
                                title(['y=',sprintf('%1.2f',y3),' m         t = ',sprintf('%1.5f',tsnap),' s']);
                            case 3
                                title(['z=',sprintf('%1.2f',z3),' m         t = ',sprintf('%1.5f',tsnap),' s']);
                        end;
                        xlabel([label1,'     v_y']);
                    case 2
                        caxis_value = caxis_value2;
                        switch p
                            case 1
                                title(['x=',sprintf('%1.2f',x3),' m         t = ',sprintf('%1.5f',tsnap),' s']);
                            case 2
                                title(['y=',sprintf('%1.2f',y3),' m         t = ',sprintf('%1.5f',tsnap),' s']);
                            case 3
                                title(['z=',sprintf('%1.2f',z3),' m         t = ',sprintf('%1.5f',tsnap),' s']);
                        end;
                        xlabel([label1,'     curl_x']);
                    case 4
                        caxis_value = caxis_value1;
                        xlabel([label1,'     \sigma_{xy}']);
                    case 5
                        caxis_value = caxis_value1;
                        switch p
                            case 1
                                title(['x=',sprintf('%1.2f',x3),' m         t = ',sprintf('%1.5f',tsnap),' s']);
                            case 2
                                title(['y=',sprintf('%1.2f',y3),' m         t = ',sprintf('%1.5f',tsnap),' s']);
                            case 3
                                title(['z=',sprintf('%1.2f',z3),' m         t = ',sprintf('%1.5f',tsnap),' s']);
                        end;
                        xlabel([label1,'     v_y']);
                end;
            case 3
                switch l
                    case 1
                        caxis_value = caxis_value1;
                        xlabel([label1,'     v_z']);
                    case 2
                        caxis_value = caxis_value2;
                        xlabel([label1,'     curl_y']);
                    case 4   
                        caxis_value = caxis_value1;
                        xlabel([label1,'     \sigma_{yy}']);
                    case 5
                        caxis_value = caxis_value1;
                        xlabel([label1,'     v_z']);
                        
                end;
            case 4
                switch l
                    case 2
                        caxis_value = caxis_value2;
                        xlabel([label1,'     curl_z']);
                    case 4
                        caxis_value = caxis_value1;
                        xlabel([label1,'     \sigma_{xz}']);
                    case 5   
                        caxis_value = caxis_value2;
                        xlabel([label1,'     pressure']);
                end;
            case 5
                switch l
                    case 4
                        caxis_value = caxis_value1;
                        xlabel([label1,'     \sigma_{zz}']);
                end;
            case 6
                switch l
                    case 4
                        caxis_value = caxis_value1;
                        xlabel([label1,'     \sigma_{yz}']);
                end;                
        end;
        
        %color map
        if m
            v(find(abs(v) < 0.2*caxis_value)) = caxis_value*model(find(abs(v) < 0.2*caxis_value));
            colormap(load('seismic_mod.map'));
        else
            colormap(load('seismic.map'));
        end;
        
        %model display
        switch p
            case 1
                imagesc(y2(1,:),z2(:,1),v);
            case 2
                imagesc(x2(1,:),z2(:,1),v);
            case 3
                imagesc(x2(1,:),y2(:,1),v);
        end;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %additional features 
%         line([400 400],[0 2000],'Color',[.5 .5 .5]);
%         line([800 800],[0 2000],'Color',[.5 .5 .5]);
%         line([1200 1200],[0 2000],'Color',[.5 .5 .5]);
%         line([1600 1600],[0 2000],'Color',[.5 .5 .5]);
%         line([0 2000],[400 400],'Color',[.5 .5 .5]);
%         line([0 2000],[800 800],'Color',[.5 .5 .5]);
%         line([0 2000],[1200 1200],'Color',[.5 .5 .5]);
%         line([0 2000],[1600 1600],'Color',[.5 .5 .5]);
%         line([80 80 1920 1920 80],[1920 80 80 1920 1920],'Color','k');
%         line([0 2000],[1000 1000],'Color','k');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %plot sources and receivers
        if k
            plot(src(:,1),src(:,2),'hk','MarkerFaceColor','k','MarkerSize',12);
            plot(recs(:,1),recs(:,2),'vk','MarkerFaceColor','k');
        end;
        hold off;
        caxis([-caxis_value caxis_value]);
              
        %end model display
        
        ylabel([label2]);
        set(gca,'DataAspectRatio',[1 1 1]);
        set(get(gca,'title'),'FontSize',14);
        set(get(gca,'title'),'FontWeight','bold');
        set(get(gca,'Ylabel'),'FontSize',14);
        set(get(gca,'Ylabel'),'FontWeight','bold');
        set(get(gca,'Xlabel'),'FontSize',14);
        set(get(gca,'Xlabel'),'FontWeight','bold');
        set(gca,'FontSize',12);
        set(gca,'FontWeight','bold');
        set(gca,'Linewidth',1.0);
        set(gca,'Box','on');
        axis(ar);
	    axis ij
        %set(gca,'YTick',[400 800 1200]) 
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    pause(0.1);
	%set(gcf,'Renderer','zbuffer')
   	%M(i-firstframe+1)=getframe(gcf);

    %brighten(0.5)

	snapno=snapno+1;
    %saving snapshot
    switch n
        case 1
            epsfile=[snapfile,'_',int2str(ii),'.eps'];
            eval(['print -deps ' epsfile]);
        case 2
            jpgfile=[snapfile,'_',int2str(ii),'.jpg'];
            eval(['print -djpeg100 ' jpgfile]);
        case 3
            ppmfile=[snapfile,'_',int2str(ii),'.ppm'];
            eval(['print -dppmraw ' ppmfile]);
        case 4
            pngfile=[snapfile,'_',int2str(ii),'.png'];
            eval(['print -dpng ' pngfile]);
        case 5
            jpgfile=[snapfile,'_',int2str(ii),'.jpg'];
            eval(['print -djpeg100 ' jpgfile]);
            pngfile=[snapfile,'_',int2str(ii),'.png'];
            eval(['print -dpng ' pngfile]);   
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch l
    case 1
        fclose(fid1);
        fclose(fid2);
        fclose(fid3);
    case 2
        fclose(fid1);
        fclose(fid2);
        fclose(fid3);
        fclose(fid4);
    case 3
        fclose(fid1);
    case 4
        fclose(fid1);
        fclose(fid2);
        fclose(fid3);
        fclose(fid4);
        fclose(fid5);
        fclose(fid6);
    case 5
        fclose(fid1);
        fclose(fid2);
        fclose(fid3);
        fclose(fid4);
end;
