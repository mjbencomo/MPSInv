% This function reads the binary file containing medium parameters for a 2D
% salt model.
clear

% Grid information
Ny = 1911;
Nx = 5395;
dx = 12.5; %units meters
dy = 6.26; % units meters

xData = (0:Nx-1)*dx;
yData = (0:Ny-1)*dy;

xDom = [0,dx*(Nx-1)];
yDom = [0,dy*(Ny-1)];

units.x = 'm';
units.y = 'm';
units.beta  = 'm^3/kg';
units.kappa = 'Pa';

%% Computing buoyancy

% Reading in density
fileID = fopen('salt2D_den.bin','rb'); 
if fileID == -1
    error('File could not be opened. Check the path.');
end
rho = fread(fileID, [Ny,Nx], 'single').'; 
fclose(fileID); 

% Computing beta = buoyancy
% Note, we transpose since we expect beta(i,j) = beta(x_i,y_i)
betaData = 1./rho;

%% Computing bulk modulus

% Reading in velocity
fileID = fopen('salt2D_vel.bin','rb'); 
if fileID == -1
    error('File could not be opened. Check the path.');
end
vel = fread(fileID, [Ny,Nx], 'single').'; 
fclose(fileID); 

% Computing kappa = bulk modulus
% Note, we transpose since we expect kappa(i,j) = kappa(x_i,y_i)
kappaData = (vel.^2)./betaData;

%% Saving data
projectRoot = currentProject().RootFolder;
matFile = fullfile(projectRoot, ...
    "src","utils","medium","mediumData",...
    "salt2D_data.mat");
save(matFile,'betaData','kappaData','xData','yData','xDom','yDom','units');