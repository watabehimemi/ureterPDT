format compact
clc
home

%% Specify Monte Carlo parameters    
myname = 'disc';
Nx = 150;           % Number of x-axis voxels
Ny = 150;           % Number of y-axis voxels
Nz = 150;           % Number of z-axis voxels
binsize = 0.1;      % [mm]

dx = binsize;
dy = binsize;
dz = binsize;

x = ([1:Nx]' - Nx/2)*dx;  % x-axis coordinate
y = ([1:Ny]' - Ny/2)*dy;  % y-axis coordinate
z = ([1:Nz]' - Nz/2)*dz;  % z-axis coordinate

zmin = min(z);
zmax = max(z);
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

%%%%%%
% CREATE TISSUE STRUCTURE T(x, y, z)
%   Create T(x, y, z) by specifying a tissue type (an integer)
%   for each voxel in T.

T = double(0 * ones(Nx, Ny, Nz));  % Initialize as air
Txy = T(:,:,7);
Tzx = T(:,7,:);


% Create cylindrical tissue layer
radius1 = 1.5;  % Radius of the air [mm]
radius2 = 3; % Radius of the mucosa  [mm]
radius3 = 7.5; % Radius of the muscle [mm]
center_x = 0; % x-coordinate of the cylinder center
center_y = 0; % y-coordinate of the cylinder center
center_z = 0; % z-coordinate of the cylinder center

for iz = 1:Nz
    for ix = 1:Nx
        for iy = 1:Ny
            r = sqrt((x(ix) - center_x)^2 + (y(iy) - center_y)^2);
            if r < radius1
                T(ix, iy, iz) = 1;  % Tissue type 0 (air)
            elseif r < radius2
                T(ix, iy, iz) = 2;  % Tissue type 1 (mucosa)
            elseif r < radius3
                T(ix, iy, iz) = 3;  % Tissue type 2 (muscle)
            end
        end
    end
end

%% convert T to linear array of integer values, v(i)
v = uint8(reshape(T, Nx * Ny * Nz, 1));

%% write myname_T.bin file
filename = sprintf('%sT.bin', myname);
disp(['create ' filename])
fid = fopen(filename, 'wb');
fwrite(fid, v, 'uint8');
fclose(fid);

%% Visualization
figure(1); clf 
imagesc(x, y, squeeze(T(:, :, Nz/2)))
hold on
fsz = 22;
set(gca, 'fontsize', fsz)
xlabel('x [mm]')
ylabel('y [mm]')
axis equal image
grid minor
map = [1 1 1
       1 1 1   % 白
       1 0.7 0.6   % オレンジ
       1 0.6 0.8];     % ピンク
colormap(map)
axis([xmin xmax ymin ymax])

figure(2); clf
imagesc(x,z,squeeze(T(:, Ny/2,:)))
hold on
set(gca,'fontsize',fsz)
xlabel('z [mm]')
ylabel('x [mm]')
axis equal image
grid minor
map = [1 1 1
       1 1 1   % 白
       1 0.7 0.6   % オレンジ
       1 0.6 0.8];     % ピンク
colormap(map)
axis([zmin zmax xmin xmax])

%% read json file
json = jsondecode(fileread('disc.json'));
disp(json);
%% Calculation

command = 'mcx -f disc.json -F mc2 -O F -w D -G 1101';
status = system(command, '-echo');

disp('done')


