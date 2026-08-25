% Viewing salt2D_data.mat 

%% Reading in data

clear
load("salt2D_data.mat");
who

%% Plotting

% Plotting bulk modulus
figure
imagesc(xData,yData,kappaData)
cb = colorbar;
cb.Label.String = 'Bulk Modulus (Pa)';
xlabel('Distance (m)')
ylabel('Depth (m)')

% Plotting buoyancy
figure
imagesc(xData,yData,betaData)
cb = colorbar;
cb.Label.String = 'Buoyancy (m^3/kg)';
xlabel('Distance (m)')
ylabel('Depth (m)')

% Plotting wave speed
figure
imagesc(xData,yData,sqrt(betaData.*kappaData))
cb = colorbar;
cb.Label.String = 'Wave Speed (m/s)';
xlabel('Distance (m)')
ylabel('Depth (m)')