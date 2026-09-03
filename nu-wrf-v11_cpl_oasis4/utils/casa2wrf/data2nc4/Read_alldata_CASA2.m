%===================================================================
% This code reads the txt data files and write them as binary files. 
% Jossy Jacob March 20, 2015
%===================================================================
%
%
function Read_alldata_CASA2()
infile = 'Original_data/FUE2010.txt'
TT = csvread(infile);
NLAT = 360; NLON = 720; nmon = 12;
FUE(1:NLON,1:NLAT,1:nmon)= 0.0;
temp(1:NLAT,1:NLON) = 0.0;
t1(1:NLAT,1:NLON) = 0.0; 
istart = 1;
for imon = 1:nmon
    i1 = istart;
    i2 = istart + NLAT-1;
    temp = TT(i1:i2,:);
    for j = 1: NLAT
        j1 = NLAT - j+1;
        t1(j,:) = temp(j1,:);
    end
    FUE(:,:,imon)=t1'; 
    istart = i2+1;
end
figure; 
for i = 1:nmon
    subplot(4,3,i) 
    tp(:,:)=FUE(:,:,i);tp(tp>12)=12;
    ttext = ['FUE(',num2str(i),')'];
    pcolor(tp');shading flat; caxis([0 12]); title(ttext);
   %  hc=colorbar('vert');set(hc,'position',[.88 .51 .05 .18])
end
 hc=colorbar('horiz');set(hc,'position',[.15 .075 .73 .02]);
pname = 'FUE.png'; printstr = ['print -dpng ',pname];eval(printstr);

filename1 = 'Binary/BFUEL_2010.dat'
write_bin3d(FUE,NLON, NLAT,imon,filename1); 
clear TT temp t1 tp FUE lon lat; 

infile = 'Original_data/NPP2010.txt'
TT = textread(infile);
NLAT = 360; NLON = 720; nmon = 12;
NPP(1:NLON,1:NLAT,1:nmon)= 0.0;
temp(1:NLAT,1:NLON) = 0.0;
t1(1:NLAT,1:NLON) = 0.0; 
istart = 1;
for imon = 1:nmon
    i1 = istart;
    i2 = istart + NLAT-1;
    temp = TT(i1:i2,:);
    for j = 1: NLAT
        j1 = NLAT - j+1 ;
        t1(j,:) = temp(j1,:);
    end
    NPP(:,:,imon)=t1'; 
    istart = i2+1; 
end
figure; 
for i = 1:nmon
    subplot(4,3,i) 
    tp(:,:)=NPP(:,:,i);
    ttext = ['NPP(',num2str(i),')'];
    pcolor(tp');shading flat; caxis([0 300]); title(ttext);
   %  hc=colorbar('vert');set(hc,'position',[.88 .51 .05 .18])
end
 hc=colorbar('horiz');set(hc,'position',[.15 .075 .73 .02]);
pname = 'NPP.png'; printstr = ['print -dpng ',pname];eval(printstr);

filename2 = 'Binary/NPP_2010.dat'
write_bin3d(NPP,NLON, NLAT,imon,filename2); 
clear TT  temp t1 tp NPP lon lat;

infile = 'Original_data/resp2010.txt'
TT = textread(infile);
NLAT = 360; NLON = 720; nmon = 12;
RESP(1:NLON,1:NLAT,1:nmon)= 0.0;
temp(1:NLAT,1:NLON) = 0.0;
t1(1:NLAT,1:NLON) = 0.0; 
istart = 1;
for imon = 1:12
    i1 = istart;
    i2 = istart + NLAT-1;
    temp = TT(i1:i2,:);
    for j = 1: NLAT
        j1 = NLAT - j+1;
        t1(j,:) = temp(j1,:);
    end
    RESP(:,:,imon)=t1'; 
    istart = i2+1; 
    [i1,i2, istart, imon];
end
%for i = 1:12
%    figure; 
%    tp(:,:)=RESP(:,:,1);
%    pcolor(tp');shading flat; colorbar; title('RESP (month = 1)')
%end
figure; 
for i = 1:nmon
    subplot(4,3,i) 
    tp(:,:)=RESP(:,:,i);tp(tp>250)=250;
    ttext = ['RESP(',num2str(i),')'];
    pcolor(tp');shading flat; caxis([0 250]); title(ttext);
   %  hc=colorbar('vert');set(hc,'position',[.88 .51 .05 .18])
end
 hc=colorbar('horiz');set(hc,'position',[.15 .075 .73 .02]);
pname = 'RESP.png'; printstr = ['print -dpng ',pname];eval(printstr);
filename2 = 'Binary/RESP_2010.dat'
write_bin3d(RESP,NLON, NLAT,imon,filename2); 
clear TT temp t1 tp lat lon RESP ;


%fname = 'RESP.dat'
%[var,nx,ny,ntime]= read_bin3d(fname)

infile = 'Original_data/fire_emiss_daily_2010.txt'
NLAT = 360; NLON = 720; nmon = 365;
TT = csvread(infile);
FIRE(1:NLON,1:NLAT,1:nmon)= 0.0;
temp(1:NLAT,1:NLON) = 0.0;
t1(1:NLAT,1:NLON) = 0.0; 
istart = 1;
for imon = 1:nmon
    i1 = istart;
    i2 = istart + NLAT-1;
    temp = TT(i1:i2,:);
    for j = 1: NLAT
        j1 = NLAT - j+1;
        t1(j,:) = temp(j1,:);
    end
    FIRE(:,:,imon)=t1'; 
    istart = i2+1;
    %[i1,i2, istart, imon]
end
%for i = 1:12
    figure; 
    tp(:,:)=FIRE(:,:,1);
    pcolor(tp');shading flat; colorbar; title('FIRE (month = 1)');
%end
figure; 
for i = 1:12
    subplot(4,3,i) 
    tp(:,:)=FIRE(:,:,i);tp(tp>10)=10;
    ttext = ['FIRE(',num2str(i),')'];
    pcolor(tp');shading flat; caxis([0 10]); title(ttext);
   %  hc=colorbar('vert');set(hc,'position',[.88 .51 .05 .18])
end
hc=colorbar('horiz');set(hc,'position',[.15 .075 .73 .02]);
pname = 'FIRE.png'; printstr = ['print -dpng ',pname];eval(printstr);
filename2 = 'Binary/FIRE_2010.dat';
write_bin3d(FIRE,NLON, NLAT,nmon,filename2); 
clear TT temp t1 tp lat lon FIRE ;
 
%
%End
% Test reading the data files: 

%fname = 'FIRE.dat'
%[var,nx,ny,ntime]= read_bin3d(fname);
%Read the taka data
% This data need to be rotated from -180 to 180
%
infile2 ='Original_data/Taka02.1x1.25.data'; 
ntime=12 ; nx=288; ny=181;
NLAT = ny; NLON = nx; nmon = ntime;
%lon0=0; lat0=-89.5; dx = 1.25; dy = 1.0;
lon0=0; lat0=-90.; dx = 1.25; dy = 1.0;
[var]= read_bin3d_nohead(infile2,ntime,nx,ny);
OCO2 = var(:,1:NLAT,:);
for i = 1: nx
lon(i) = lon0 + dx* (i-1);
end
for j = 1:ny
  lat(j) = lat0 + dy*(j-1);
end
[x,y]=meshgrid(lon,lat');
figure;

tp(:,:)=var(:,1:NLAT,1)* 1.e9;%tp(tp>10)=10;
pcolor(x,y,tp');shading flat;caxis([-3 3]);
%ttext = ['Ocean CO2 monthly'];
%text(.3, 1.02, ttext, 'FontSize',12, 'FontWeight', 'bold', 'Color','b');

m1 = 4; n1 =3 ;
for i = 1:ntime
    subplot(m1,n1,i)     
    tp(:,:)=var(:,1:NLAT,i)* 1.e9;%tp(tp>10)=10;
    ttext = ['month(',num2str(i),')'];
     pcolor(x,y,tp');shading flat;caxis([-3 3]);
    %title(ttext)
end
ttex = 'Ocean CO2 monthly * 1.0e9';
text(-600,870,ttex,'FontSize',12, 'FontWeight', 'bold', 'Color','b');
hc=colorbar('horiz');set(hc,'position',[.2 .05 .7 .02]);
print -dpng taka.png;
filename2 = 'Binary/OceanCO2_2010.dat'
write_bin3d(OCO2,NLON, NLAT,nmon,filename2); 
clear tp

%% TO plot rotated the field 180 deg longitude: 
n1 = nx/2 + 1;
lon0 = lon(n1)-360 + dx/2;
for i = 1:nx
lon1(i) = lon0 + dx *(i-1);
end

[x,y]=meshgrid(lon1,lat');
for nt = 1:ntime
var2(1:n1-1,:)=var(n1:nx,:,nt);
var2(n1:nx,:)=var(1:n1-1,:,nt);
OC2(:,:,nt)=var2(:,:);
end
figure;
m1 = 4; n1 =3 ;
for i = 1:ntime
    subplot(m1,n1,i) ;    
    tp(:,:)=OC2(:,1:NLAT,i)* 1.e9;%tp(tp>10)=10;
    ttext = ['month(',num2str(i),')'];
     pcolor(x,y,tp');shading flat;caxis([-3 3]);
    %title(ttext)
end
ttex = 'Ocean CO2 monthly * 1.0e9';
text(-600,870,ttex,'FontSize',12, 'FontWeight', 'bold', 'Color','b');
hc=colorbar('horiz');set(hc,'position',[.2 .05 .7 .02]);
print -dpng taka_rotated.png;
clear lon lat
%% Plot Fossil Fuel climatology 
%
fname = 'Original_data/D_float.2010' 
%Not sure whether ny=181 or 180
ntime=1; nx=288; ny=181;NLAT = ny; NLON = nx; nmon = ntime;
lon0=0; lat0=-90.; dx = 1.25; dy = 1.0;
[var]= read_bin3d_nohead(fname,ntime,nx,ny);
for i = 1: nx
lon(i) = lon0 + dx* (i-1);
end
for j = 1:ny
  lat(j) = lat0 + dy*(j-1);
end
[x,y]=meshgrid(lon,lat');
tp(:,:)=var(:,1:NLAT)*1.e7;%tp(tp>10)=10;

figure; 
 ttext = ['Fossil Fuel yearly climatology *1e7'];
 pcolor(x,y,tp');shading flat;caxis([0 1]);
 title(ttext,'FontSize',12, 'FontWeight', 'bold', 'Color','b'); 
 colorbar; 
 print -dpng D_float.png ;

filename2 = 'Binary/FossilFuel_2010.dat'
write_bin3d(var(1:NLON,1:NLAT,1:nmon),NLON, NLAT,nmon,filename2); 
clear tp lat lon car
end

%====================================================
function [var,nx,ny,ntime] = read_bin3d(fname)
fid      = fopen(fname,'r','b');
dummy = fread(fid,1,'int');
head = fread(fid,[1,3],'int');
nx = head(1); ny = head(2); ntime = head(3);
dummy = fread(fid,1,'int');
for it = 1:ntime 
   dummy = fread(fid,1,'int');
   temp=fread(fid,[nx,ny],'float');
   dummy = fread(fid,1,'int');
   var(:,:,it)=temp;
end
fclose(fid)
%return 
end

%====================================================
function [var] = read_bin3d_nohead(fname,ntime,nx,ny)
fid      = fopen(fname,'r','b');
for it = 1:ntime 
   %dummy = fread(fid,1,'int');
   temp=fread(fid,[nx,ny],'float');
   %dummy = fread(fid,1,'int');
   var(:,:,it)=temp;
end
fclose(fid)
%return
end
%====================================================
function write_bin3d(var,nx, ny, ntime,fname)
fid      = fopen(fname,'w','b');
head = [nx ny ntime];
n = size(head,2)* 4;
fwrite(fid,n,'int');
fwrite(fid,head,'int');
fwrite(fid,n,'int');
d1 = nx * ny * 4;
for it = 1:ntime 
   dt = var(:,:,it) ;
   fwrite(fid,d1,'int');
   fwrite(fid,dt,'float');
   fwrite(fid,d1,'int');
end
fclose(fid)
%return
end
%====================================================
