%% Init
clc
close all

% load data
load("EPN_data.mat");

pos = EPN_data(:,1:3); % m
vel = EPN_data(:,4:6); % mm/year
ellipsoid = referenceEllipsoid('wgs84');

% length
len = length(pos);

% calculate lat, lon
[lon,lat,~] = cart2ellip(pos(:,1),pos(:,2),pos(:,3)); % degree

vE = zeros(len,1);
vN = zeros(len,1);

% tramsformation for velocity
for i = 1:len
    R = rotation(-(90 + lon(i)),'z') * rotation(-(90 - lat(i)),'x');
    ENU = R \ [vel(i,1);vel(i,2);vel(i,3)];
    vE(i,1) = ENU(1);
    vN(i,1) = ENU(2);
end

%%
Y = [vE;vN];

A = [zeros(len,1),pos(:,3),-pos(:,2);
     -pos(:,3),zeros(len,1),pos(:,1),];

Qee = eye(len * 2) * 0.1;

%%
%
Qs1s1_e = zeros(len);
for i = 1:len
    for j = 1:len
        d = distance(lat(i),lon(i),lat(j),lon(j),ellipsoid)/1000; % km 
        Qs1s1_e(i,j) = signalCovariance(d);
    end
end

Qs1s1_n = Qs1s1_e;
Qs1s1 = [Qs1s1_n, zeros(len); zeros(len), Qs1s1_e];

Qyy = Qs1s1 + Qee;

%



% to interpolate
[lat_new, lon_new] = meshgrid(40:60,10:20);
lat_new = lat_new(:);
lon_new = lon_new(:);

% ecef coordinates of new points

[x_new, y_new, z_new] = ellip2cart(lon_new/180*pi,lat_new/180*pi);
len2 = length(x_new);

Qs1s2_e = zeros(len2,len);
for i = 1:len
    for j = 1:len2
        d = distance(lat(i),lon(i),lat_new(j),lon_new(j),ellipsoid)/1000; % degree
        Qs1s2_e(j,i) = signalCovariance(d);
    end
end
Qs1s2_n = Qs1s2_e;
Qs1s2 = [Qs1s2_n,zeros(len2,len);zeros(len2,len),Qs1s2_e];
Qsy = Qs1s2;

% 
Qs2s2_e = zeros(len2);
for i = 1:len2
    for j = 1:len2
        d = distance(lat_new(i),lon_new(i),lat_new(j),lon_new(j),ellipsoid)/1000; % km 
        Qs2s2_e(i,j) = signalCovariance(d);
    end
end
Qs2s2_n = Qs2s2_e;
Qs2s2 = [Qs2s2_n,zeros(len2);zeros(len2),Qs2s2_e];


x_dach = inv(A' * inv(Qyy) * A) * A' * inv(Qyy) * Y;
s_dach = Qsy * inv(Qyy) * (Y - A * x_dach);

A_new = [zeros(len2,1),z_new,-y_new;
     -z_new,zeros(len2,1),x_new];

Y_new = A_new * x_dach + s_dach;

vE_new = Y_new(1:len2);
vN_new = Y_new(len2+1:end);


%% residual
Y1 = Y - A * x_dach;
Y2 = s_dach;

vE_r1 = Y1(1:len);
vN_r1 = Y1(len+1:end);
vE_r2 = Y2(1:len2);
vN_r2 = Y2(len2+1:end);


%% plot 
figure
worldmap('Europe')
land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(land, 'FaceColor', [1 1 1])
hold on
quiverm(lat,lon,vN,vE)
quiverm(lat_new,lon_new,vN_new,vE_new,'r')

figure
worldmap('Europe')
land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(land, 'FaceColor', [1 1 1])
hold on
quiverm([lat;lat_new,],[lon;lon_new],[vN;vN_new],[vE;vE_new])

figure
worldmap('Europe')
land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(land, 'FaceColor', [1 1 1])
hold on
quiverm(lat,lon,vN_r1,vE_r1)
quiverm(lat_new,lon_new,vN_r2,vE_r2,'r')

figure
worldmap('Europe')
land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(land, 'FaceColor', [1 1 1])
hold on
quiverm([lat;lat_new,],[lon;lon_new],[vN_r1;vN_r2],[vE_r1;vE_r2])



%%
H = Qsy * inv(Qyy);
G = inv(A' * inv(Qyy) * A) * A' * inv(Qyy);
Q_x_dach = inv(A' * inv(Qyy) * A);
Q_s_dach = Qs2s2 - Qsy * inv(Qyy) * Qsy'+ H * A * Q_x_dach * A' * H';
Q_t = Qs2s2 - Qsy * inv(Qyy) * Qsy'+ (H * A - A_new) * Q_x_dach * (A'* H' - A_new');