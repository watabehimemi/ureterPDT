% lookmcxyz.m
%   Looks at myname_F.bin, created by mcxyz.c 
%   where myname is the name of the run: myname_T.bin, myname_H.mci
%
% lookmcxyz_alc2.m, July 23
%       Absorption array is OFF.
%
%   Simulates angiolight catheter within 4-mm-dia. vessel
% with the vessel wall at a particular musp value. 
%   For this run, musp = 200 cm^-1.
%
% Reads 8 mcxyz runs, for 8 catheter positions, zs = 0.2 to 0.6 cm.
% For each run:
%   alc#_H.mci --> header:
%       timin,Nx,Ny,Nz,dy,dx,dz,xs,ys,zx,Nt,muav(),musv(),gv()
%   alc#_F.bin --> F(y,x,z) = relative fluence rate [1/cm^2]
%   alc#_T.bin --> T(y,x,z) = tissue types
%
% Displays
%   Tzx = end-view of tissue structure
%   Fzx = end-view of vessel @ source, ys = -0.15 cm
%   Fzy = side-view along length of vessel
%
% Saves
%   Fzy_data4.mat = Fzy y z zzs Fdet
%       Fzy(400,400,8) = 8 z,y images
%       Fdet(8,1) = signal [1/cm^2] @ detector fiber
%

home
clear
format compact
commandwindow

cc = 'rbgm'; % color

%%%% USER CHOICES <---------- you must specify -----
myname = 'disc';
nm     = 800;
%%%%

Nx = 150;
Ny = 150;
Nz = 150;
binsize = 0.1; %[mm]
dx = binsize;
dy = binsize;
dz = binsize;
%% Absorption A(y,x,z) 
filename = sprintf('%s.mc2',myname);
disp(['loading ' filename])

tic
fid = fopen(filename, 'rb');
if fid == -1
    error('ファイルを開けませんでした。');
end

[Data, count] = fread(fid, Ny*Nx*Nz, 'float');
fclose(fid);
toc

fid2 = fopen(filename, 'rb');
dataCalc = fread(fid2, inf, 'float');
fclose(fid2);

dataCalc = reshape(dataCalc, [Nx, Ny, Nz]);
Fl = reshape(Data, Nx, Ny, Nz); % Absorption(y,x,z)
Ab = Fl * 784.8;
%%
x = (0:Nx) * dx;  % mm
y = (0:Ny) * dy;
z = (0:Nz) * dz;
ux = 2:Nx-1;
uy = 2:Ny-1;
uz = 2:Nz-1;
zmin = min(z);
zmax = max(z);
zdiff = zmax - zmin;
xmin = min(x);
xmax = max(x);
xdiff = xmax - xmin;

%% Look at Fluence, Temperature, Thermal damage
disp(size(Ab(Nx/2,:,:))); % yz
disp(size(Ab(:,Ny/2,:))); % xz
disp(size(Ab(:,:,Ny/2))); % xy

Flxy = reshape(Fl(:,:,Nz/2), Nx, Ny)';

Flyzx = shiftdim(Fl, 1);   % Txyz --> Tyzx 
disp(size(Flyzx(:,:,Nx/2)));
Flyz = reshape(Flyzx(:,:,Nx/2), Ny, Nz)';

Flzxy = shiftdim(Fl, 2);   % Txyz --> Tzxy 
disp(size(Flzxy(:,:,Ny/2)));
Flxz = reshape(Flzxy(:,:,Ny/2), Nz, Nx)';
%% Fz, Fx

fsz = 14;
figure(1); clf
imagesc(x, z, log10(Flxz), [-5 0])
hold on
colorbar
set(gca, 'fontsize', fsz)
set(gca, 'fontname', 'Arial')
xlabel('\itz \rm[mm]')
ylabel('\itx \rm[mm]')
ylim([-1,1])
colormap('jet')
axis equal image
grid minor

figure(2); clf
imagesc(y, z, log10(Flyz), [-5 0])
hold on
colorbar
set(gca, 'fontsize', fsz)
set(gca, 'fontname', 'Arial')
xlabel('\ity \rm[mm]')
ylabel('\itz \rm[mm]')
ylim([-1,1])
colormap('jet')
axis equal image
grid minor

figure(3); clf
imagesc(x, y, log10(Flxy), [-5 0])
hold on
colorbar
set(gca, 'fontsize', fsz)
set(gca, 'fontname', 'Arial')
xlabel('\itx \rm[mm]')
ylabel('\ity \rm[mm]')
ylim([-1,1]);
colormap('jet')
axis equal image
grid minor

% figure(4);
% xslice = Nx/2;   
% yslice = Ny/2;
% zslice = Nz/2;
% t = slice(Fl, xslice, yslice, zslice);
% set(t, 'EdgeColor', 'none');
% alpha('color');
% alphamap('vup');
% xlabel('y [mm]')
% ylabel('x [mm]')
% zlabel('z [mm]')

%%
% figure(5);
% cmap = parula(256);
% sliceViewer(Fl, 'Colormap', cmap, 'SliceDirection', [1 0 0]);
% 
% figure(6);
% cmap = parula(256);
% sliceViewer(Fl, 'Colormap', cmap, 'SliceDirection', [0 1 0]);
% 
% figure(7);
% cmap = parula(256);
% sliceViewer(Fl, 'Colormap', cmap, 'SliceDirection', [0 0 1]);