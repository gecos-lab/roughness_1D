%% Processing and power spectrum analysis of quasi-linear roughness profiles
%  Andrea Bistacchi 16/11/2009
%  update 24/04/2014 -> SCALING PSD - AT LAST! - then histogram
%  update 23/0372019 -> removed old matlab funcions like princomp -> pca

function ROUGHNESS_1D_08

clear; close all; clc; hold off

set(0,'DefaultFigureWindowStyle','docked','DefaultFigureColor','w');
set(0,'DefaultAxesFontName','Times','DefaultAxesFontSize',16);
set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on');

disp(' ');
disp('Processing and power spectrum analysis of quasi-linear roughness profiles');

figure(1);
figure(2);
figure(3);
figure(4);
figure(5);
figure(6);
figure(7);
figure(8);
figure(9);
figure(10);
figure(11);
figure(12);
figure(13);
figure(14);
figure(15);
figure(16);
figure(17);
figure(18);
figure(19);
figure(20);

displ = NaN; % needed to inizialize displacement
i=0; % needed for storing  notes about plots
plotNotes =[{'plot No'}, {'monoDimension'},...
    {'plane normal D/D'},{'slip vector P/T'},...
    {'profileLenght'},...
    {'displ'},...
    {'profile Min'}, {'profile Max'},...
    {'trace vector'}, {'slip2trace'}, {'sampling'}, {'taper'},...
    {'spectrum method'}, {'lower freq'}, {'upper freq'}, {'H'}, {'P1'}, {'sigmaRMS'}, {'Rsquare'},...
    {'file name'}];


% Root menu
while 1
    
    disp(' ');
    disp('__________');
    disp('ROOT MENU:');
    disp('1- Load rougness data from xyz file');
    disp('2- Preprocess data (clean and choose monotonic dimension)');
    disp('3- Define interpolating plane with Principal Component Analysis');
    disp('4- Input displacement');
    disp('5- Select "almost straight" segment from profile');
    disp('6- Calculate distance to plane');
    disp('7- Power spectrum analysis');
    disp('8- Save processed data');
    disp('9- Load previously saved data');
    disp('10- Power spectrum summary plot (load last results to plot)');
    disp('11- Power spectrum summary plot of multiple saved files');
    disp('12- Quit');
    disp(' ');
    
    action = 0;
    while (action < 1) || (action > 12)
        action = round(input(' > '));
    end
    
    switch action
        case 1
            [data,filename] = LOAD;
        case 2
            [x,y,z,monoDimension] = PREPROCESS(data);
        case 3
            [planeNormal,slipVector] = PLANE(x,y,z);
        case 4
            displ = DISPL(planeNormal,slipVector);
        case 5
            [xi,yi,zi,Min,Max] = SELECT(x,y,z,monoDimension);
        case 6
            [t,dist,tLS,distLS,sample,traceVector,slip2trace,profileLenght] = LINE(xi,yi,zi,monoDimension,planeNormal,slipVector);
        case 7
            [Amplitude,Freq,method,lowerF,upperF,H,P1,sigmaRMS,Rsquare,taperPercent] = SPECTRUM(t,dist,tLS,distLS);
        case 8
            SAVE(Amplitude,Freq,H,P1,sigmaRMS,slip2trace,monoDimension,planeNormal,slipVector,Min,Max,traceVector,sample,taperPercent,method,lowerF,upperF,filename,displ,profileLenght,Rsquare);
        case 9
            disp(' ');
            disp('9- Load previously saved data from .mat file');
            
            clear monoDimension planeNormal slipVector profileLenght displ Min Max traceVector slip2trace sample taperPercent method lowerF upperF H P1 sigmaRMS filename;
            
            sample = 0; % needed because of conflict with sample.m function
            displ = 0;  % needed for joints where displ is not defined
            profileLenght = -1;  % needed for backward compatibility
            
            [file, path] = uigetfile('*.mat');
            load([path file]);
            
            if exist('P1','var')==0, % needed for backward compatibility with files where P1 wasn't saved
                P1 = A10/(10^-(1+2*H));
            end
            
            fMax = 1e2; fMin = 1e-2; %  -> note fMax and fMin
            sigmaRMS = sqrt(P1 * (fMax^-(2*H) - fMin^-(2*H)) / -(2*H)); % needed for backward compatibility with files where sigmaRMS was wrong
            
            disp(' ');
            disp([' -> file ' file ' successfully loaded.']);
            
        case 10
            i = i+1;
            thisPlotNotes = [{num2str(i)}, {num2str(monoDimension)},...
                {dipdirstr(planeNormal)},{plungetrendstr(slipVector)},...
                {num2str(profileLenght)},...
                {num2str(displ)},...
                {num2str(Min)}, {num2str(Max)},...
                {plungetrendstr(traceVector)}, {num2str(slip2trace)}, {num2str(sample)}, {num2str(taperPercent)},...
                {num2str(method)}, {num2str(lowerF)}, {num2str(upperF)}, {num2str(H)}, {num2str(P1)}, {num2str(sigmaRMS)}, {num2str(Rsquare)},...
                {filename}];
            
            plotNotes = [plotNotes; thisPlotNotes];
            
            SUMMARY(Amplitude,Freq,H,P1,sigmaRMS,slip2trace,plotNotes,displ,lowerF,upperF);
        case 11
            MULTIPLOT;
        case 12
            break
    end
end

% set(0,'DefaultFigureWindowStyle','normal');

end
%%

%% 1- import x y z data (e.g. generated by gOcad -> export as Excel table).
%     Data should be stored in an ASCII file with x, y, z columns, and no header.
function [data,filename] = LOAD;

disp(' ');
disp('1- Load rougness data from xyz file');
disp('   Data should be stored in an ASCII file with x, y, z columns, and no header.');

[filename, pathname] = uigetfile('*.*');
data = load([pathname filename]);

disp(' ');
disp([' -> ' num2str(max(size(data))) ' rows successfully loaded.']);

end
%%

%% 2- clean data and choose which dimension must be monotonic
%     Input data is a n x 3 martix with 1st, 2nd and 3rd columns corrensponding
%     to x, y, z coordinates. This can be extracted from gOcad with "export to
%     Excel" for point or line objects.
%     Input monoDimension specifies which dimension must be monotonic: 1 => x,
%     2 => y, 3 => z.
%     Output are three column vectors with x, y, z coordinates. The chosen
%     monotonic dimension is strictly monotonically increasing, with no
%     duplicate values.
function [x,y,z,monoDimension] = PREPROCESS(data);

disp(' ');
disp('2- Preprocess data (clean and choose monotonic dimension)');
disp(' ');
disp('   Choose the monotonic dimension');
disp('   1 -> x');
disp('   2 -> y');
disp('   3 -> z');

monoDimension = 0;
while (monoDimension < 1) | (monoDimension > 3);
    monoDimension = round(input(' > '));
end

disp(' ');
disp(['   Monotonic dimension is ' num2str(monoDimension)]);

sortedData = sortrows(data,monoDimension);

isDuplicate = sortedData(1:end-1,monoDimension)==sortedData(2:end,monoDimension);

duplicateRowIndex = find(isDuplicate==1);

numDuplicates = length(duplicateRowIndex);

sortedData(duplicateRowIndex,monoDimension) = max(sortedData(:,monoDimension))*10*ones(numDuplicates,1);

reSortedData = sortrows(sortedData,monoDimension);

uniqueData = reSortedData(1:end-numDuplicates,:);

x = uniqueData(:,1); y = uniqueData(:,2); z = uniqueData(:,3);

disp(' ');
disp(' -> Preprocessing successfully completed.');
disp([' -> ' num2str(numDuplicates) ' duplicate rows deleted.']);

end
%%

%% 3- if suitable, find interpolating plane with Principal Component Analysis
%     otherwise, find interpolating line with PCA and define plane with an
%     arbitrary dip angle or direction
function [planeNormal,slipVector] = PLANE(x,y,z);

disp(' ');
disp('3a- Define interpolating plane');

xyz = [x y z];
%[coeff,~,roots] = princomp(xyz);
[coeff,~,roots] = pca(xyz);

firstEigenVector = coeff(:,1);
secondEigenVector = coeff(:,2);
thirdEigenVector = coeff(:,3);

firstEigenValue = roots(1);
secondEigenValue = roots(2);
thirdEigenValue = roots(3);

firstEigenValuePercent = 100*firstEigenValue/sum(roots);
secondEigenValuePercent = 100*secondEigenValue/sum(roots);
thirdEigenValuePercent = 100*thirdEigenValue/sum(roots);

disp(' ');
disp('    Interpolating plane parameters:');
disp('    Eigenvectors [plunge / trend]');
disp(['    1st ' plungetrendstr(firstEigenVector)]);
disp(['    2nd ' plungetrendstr(secondEigenVector)]);
disp(['    3rd ' plungetrendstr(thirdEigenVector)]);
disp(' ');
disp('    Eigenvalues [value / percent]');
disp(['    1st ' num2str(firstEigenValue)  ' / ' num2str(firstEigenValuePercent) ]);
disp(['    2nd ' num2str(secondEigenValue) ' / ' num2str(secondEigenValuePercent)]);
disp(['    3rd ' num2str(thirdEigenValue)  ' / ' num2str(thirdEigenValuePercent) ]);

if thirdEigenValuePercent >= 50,
    disp(' ');
    disp('    Data distribution almost spherical.');
elseif secondEigenValuePercent >= 50,
    disp(' ');
    disp('    Data distribution = oblate ellipsoid.');
else
    disp(' ');
    disp('    Data distribution = prolate ellipsoid.');
end

disp(' ');
disp('    Interpolating plane can be defined as:');
disp('    1- perpendicular to third eigenvector (good for oblate ellipsoid)');
disp('    2- containing first eigenvector, with arbitrary dip angle (good for horizontal prolate ellipsoid)');
disp('    3- containing first eigenvector, with arbitrary direction (good for vertical prolate ellipsoid)');
disp('    4- quit (gives error if plane has not been calculated yet)');
disp(' ');

method = 0;
while (method < 1) | (method > 4);
    method = round(input(' > '));
end

if     method == 1,
    planeNormal = thirdEigenVector;
    
elseif method == 2,
    L1 = firstEigenVector(1);  % direction cosines of unit vector from first eigenvector
    L2 = firstEigenVector(2);
    L3 = firstEigenVector(3);
    
    if L1 == 0 & L2 == 0, disp(' -> FIRST EIGENVECTOR IS VERTICAL'); end
    
    disp(' ');
    disp('    Input dip angle of interpolating plane');
    planeDip = input('    > ');
    
    N3 = -cos(planeDip*pi/180); % direction cosine of downward normal to plane defined by dip angle
    
    A = L1^2 + L2^2;        % 2nd degree equation
    B = 2*L1*L3*N3;
    C = (L2^2 + L3^2)*N3^2 - L2^2;
    
    N1a = ( -B + sqrt(B^2 - 4*A*C) ) / (2*A);    % 1st solution
    N2aa =  sqrt(1 - N1a^2 - N3^2);
    N2ab = -sqrt(1 - N1a^2 - N3^2);
    N2ac = (-L1*N1a - L3*N3) / L2;
    planeNormalAA = [N1a; N2aa; N3];
    planeNormalAB = [N1a; N2ab; N3];
    planeNormalAC = [N1a; N2ac; N3];
    
    N1b = ( -B - sqrt(B^2 - 4*A*C) ) / (2*A);    % 2nd solution
    N2ba =  sqrt(1 - N1b^2 - N3^2);
    N2bb = -sqrt(1 - N1b^2 - N3^2);
    N2bc = (-L1*N1b - L3*N3) / L2;
    planeNormalBA = [N1b; N2ba; N3];
    planeNormalBB = [N1b; N2bb; N3];
    planeNormalBC = [N1b; N2bc; N3];
    
    disp(' ');
    disp('    Choose interpolating plane [dip / direction / dot product]');  % THIS IS NECESSARY TO AVOID DIFFERENT PROBMLEMS WITH +o- ROOTS AND DIVISION BY 0
    disp(['    1st ' dipdirstr(planeNormalAA) ' / ' num2str(planeNormalAA'*firstEigenVector)]);
    disp(['    2nd ' dipdirstr(planeNormalAB) ' / ' num2str(planeNormalAB'*firstEigenVector)]);
    disp(['    3nd ' dipdirstr(planeNormalAC) ' / ' num2str(planeNormalAC'*firstEigenVector)]);
    disp(['    4rd ' dipdirstr(planeNormalBA) ' / ' num2str(planeNormalBA'*firstEigenVector)]);
    disp(['    5th ' dipdirstr(planeNormalBB) ' / ' num2str(planeNormalBB'*firstEigenVector)]);
    disp(['    6th ' dipdirstr(planeNormalBC) ' / ' num2str(planeNormalBC'*firstEigenVector)]);
    plane = 0;
    while (plane < 1) | (plane > 6);
        plane = round(input('    > '));
    end
    
    if     plane == 1, planeNormal = planeNormalAA;
    elseif plane == 2, planeNormal = planeNormalAB;
    elseif plane == 3, planeNormal = planeNormalAC;
    elseif plane == 4, planeNormal = planeNormalBA;
    elseif plane == 5, planeNormal = planeNormalBB;
    elseif plane == 6, planeNormal = planeNormalBC;
    end
    
elseif method == 3,
    
end

disp(' ');
disp([' -> Plane [dip / direction]: ' dipdirstr(planeNormal)]);

% plane baricenter

meanX = mean(x);
meanY = mean(y);
meanZ = mean(z);

disp(' ');
disp([' -> Plane baricenter x, y, z [m]: ' num2str(meanX) ' , ' num2str(meanY) ' , ' num2str(meanZ)]);

% slip vector

disp(' ');
disp('3b- Define slip vector');
disp('    1- as plunge / trend');
disp('    2- as -180°/+180° rake (Aki & Richards, 1980)');
disp(' ');

method = 0;
while (method < 1) | (method > 2);
    method = round(input(' > '));
end

if method == 1
    disp(' ');
    disp('    Input slip vector plunge');
    slipPlunge = input('    > ');
    disp('    Input slip vector trend');
    slipTrend = input('    > ');
else
    disp(' ');
    disp('             90 = thrust');
    disp('                |');
    disp('     0 = left --+-- 180 = right');
    disp('                |');
    disp('            -90 = normal');
    disp(' ');
    disp('    Input slip vector -180°/+180° rake (Aki & Richards)');
    rake = input('    > ');
    [planeDip,planeDir] = dipdir(planeNormal);
    [slipPlunge,slipTrend] = rake2slip(planeDip,planeDir,rake);
end

slipVector = [ sin(slipTrend*pi/180)*cos(slipPlunge*pi/180) ; cos(slipTrend*pi/180)*cos(slipPlunge*pi/180) ; -sin(slipPlunge*pi/180) ];

disp(' ');
disp([' -> Slip vector [plunge / trend]: ' plungetrendstr(slipVector)]);

end
%%

%% 4- input displacement
function displ = DISPL(planeNormal,slipVector);

disp(' ');
disp('4- Input displacement [m]');
disp(' ');
disp([' -> Plane [dip / direction]: ' dipdirstr(planeNormal)]);
disp([' -> Slip vector [plunge / trend]: ' plungetrendstr(slipVector)]);
disp(' ');

displ = input('Displacement [m] > ');

end
%%

%% 5- select "almost straight" segment(s) from profile -> selectsegment.m
function [xi,yi,zi,Min,Max] = SELECT(x,y,z,monoDimension);

disp(' ');
disp('5- Select "almost straight" segment from profile');
close(1); figure(1);

if monoDimension == 1,
    subplot(2,1,1); plot(x,y); xlabel('x [m]'); ylabel('y [m]'); axis equal;
    subplot(2,1,2); plot(x,z); xlabel('x [m]'); ylabel('z [m]'); axis equal;
    disp(' ');
    disp('   Input lower boundary (minimum x)');
    xMin = input('   > ');
    disp('   Input upper boundary (maximum x)');
    xMax = input('   > ');
    
    [dunmmy,minIndex] = min(abs(x-xMin));
    [dunmmy,maxIndex] = min(abs(x-xMax));
    
    xi = x(minIndex:maxIndex,1);
    yi = y(minIndex:maxIndex,1);
    zi = z(minIndex:maxIndex,1);
    
    subplot(2,1,1); hold on; plot(xi(1),yi(1),'or'); hold on; plot(xi(end),yi(end),'or'); axis equal;
    subplot(2,1,2); hold on; plot(xi(1),zi(1),'or'); hold on; plot(xi(end),zi(end),'or'); axis equal;
    
    Min = min(xi); Max = max(xi);
    
elseif monoDimension == 2,
    subplot(2,1,1); plot(y,x); xlabel('y [m]'); ylabel('x [m]'); axis equal;
    subplot(2,1,2); plot(y,z); xlabel('y [m]'); ylabel('z [m]'); axis equal;
    disp(' ');
    disp('   Input lower boundary (minimum y)');
    yMin = input('   > ');
    disp('   Input upper boundary (maximum y)');
    yMax = input('   > ');
    
    [dunmmy,minIndex] = min(abs(y-yMin));
    [dunmmy,maxIndex] = min(abs(y-yMax));
    
    xi = x(minIndex:maxIndex,1);
    yi = y(minIndex:maxIndex,1);
    zi = z(minIndex:maxIndex,1);
    
    subplot(2,1,1); hold on; plot(yi(1),xi(1),'or'); hold on; plot(yi(end),xi(end),'or'); axis equal;
    subplot(2,1,2); hold on; plot(yi(1),zi(1),'or'); hold on; plot(yi(end),zi(end),'or'); axis equal;
    
    Min = min(yi); Max = max(yi);
    
elseif monoDimension == 3,
    subplot(2,1,1); plot(z,x); xlabel('z [m]'); ylabel('x [m]'); axis equal;
    subplot(2,1,2); plot(z,y); xlabel('z [m]'); ylabel('y [m]'); axis equal;
    disp(' ');
    disp('   Input lower boundary (minimum z)');
    zMin = input('   > ');
    disp('   Input upper boundary (maximum z)');
    zMax = input('   > ');
    
    [dunmmy,minIndex] = min(abs(z-zMin));
    [dunmmy,maxIndex] = min(abs(z-zMax));
    
    xi = x(minIndex:maxIndex,1);
    yi = y(minIndex:maxIndex,1);
    zi = z(minIndex:maxIndex,1);
    
    subplot(2,1,1); hold on; plot(zi(1),xi(1),'or'); hold on; plot(zi(end),xi(end),'or'); axis equal;
    subplot(2,1,2); hold on; plot(zi(1),yi(1),'or'); hold on; plot(zi(end),yi(end),'or'); axis equal;
    
    Min = min(zi); Max = max(zi);
    
end

disp(' ');
disp(['   Lower boundary = ' num2str(Min)]);
disp(['   Upper boundary = ' num2str(Max)]);

end
%%

%% 6- find interpolating line segment(s) with leat squares and...
function [t,dist,tLS,distLS,sample,traceVector,slip2trace,profileLenght] = LINE(xi,yi,zi,monoDimension,planeNormal,slipVector);

disp(' ');
disp('6- Calculate distance to plane');
close(1); figure(1);

% 5.1- find regression line with Principal Component Analysis

xiyizi = [xi yi zi];
%[coeff,~,roots] = princomp(xiyizi);
[coeff,~,roots] = pca(xiyizi);

firstEigenVector = coeff(:,1);
firstEigenValue = roots(1);
firstEigenValuePercent = 100*firstEigenValue/sum(roots);

disp(' ');
disp('    Interpolating line parameters:');
disp('    Eigenvector [plunge / trend]');
disp(['    ' plungetrendstr(firstEigenVector)]);
disp(' ');
disp('    Eigenvalue [value / percent]');
disp(['    ' num2str(firstEigenValue)  ' / ' num2str(firstEigenValuePercent) ]);

% 5.1.5- adjust best fit projection plane
disp(' ');
disp('Would you like to automatically adjust the best fit plane?');
disp('It is a good idea! 0 => no; 1 => yes');

ADJflag = -1;
while (ADJflag < 0) | (ADJflag > 1);
    ADJflag = round(input(' > '));
end

if ADJflag == 1,

    first = [xi(1); yi(1); zi(1)];
    last = [xi(end); yi(end); zi(end)];
    first2last = last - first;          % vector from first to last point in selected profile
    AuxPoint = first + first2last./2 + cross(first2last,planeNormal); % auxiliary point defined as:
    
    %       first ---------  last
    %           \           /
    %            \         /
    %             \       /
    %              \     /
    %               \   /
    %             AuxPoint

    xayaza = [xi yi zi ; repmat(AuxPoint',length(xi),1) + rand(length(xi),3).*sqrt(first2last'*first2last)*.05.*-0.005];   % add AuxPoint to PCA matrix
    %[coeff,~,~] = princomp(xayaza);   % solve... this gives the best fit plane through AuxPoint
    [coeff,~,~] = pca(xayaza);   % solve... this gives the best fit plane through AuxPoint

    % figure(1); plot3(xayaza(:,1),xayaza(:,2),xayaza(:,3),'.'); % JUST TO CHECK
    
    thirdEigenVector = coeff(:,3);   % vector perpendicular to plane

    planeNormal = thirdEigenVector;  % update planeNormal
    
    disp(' ');
    disp('Updated best fit plane is:');
    disp([' -> [dip / direction]: ' dipdirstr(planeNormal)]);

end

% 5.2- projection on interpolating plane
% (P = projection of first eigenvector on plane)
% (s = distances (position vector) from curve origin)
% (sPlane = distance parameter projected along P)

P = cross(planeNormal,cross(firstEigenVector,planeNormal));  % P = N x (R x N)

traceVector = P;
slip2trace = acos(dot(P,slipVector))*180/pi; if slip2trace > 90, slip2trace = 180 - slip2trace; end % this isn't 100% accurate because slipVector should be recalculated when adjusting the best fit plane

disp(' ');
disp('    Projection of first eigenvector on plane [plunge / trend / modulus]');
disp(['    ' plungetrendstr(P/(P'*P)) ' / ' num2str((P'*P))]);
disp(' ');
disp('    Slip vector');
disp(['    ' plungetrendstr(slipVector)]);
disp(' ');
disp('    Angle to slip vector');
disp(['    ' num2str(slip2trace)]);

s = [xi - xi(1) yi - yi(1) zi - zi(1)];

sPlane = s*P;

% 5.3- calculate dist from plane
baricenter = [mean(xi); mean(yi); mean(zi)];

disti = ((xi - baricenter(1))*planeNormal(1)...
    + (yi - baricenter(2))*planeNormal(2)...
    + (zi - baricenter(3))*planeNormal(3));

% 5.4- resample to uniformly spaced parameter t
% (t = uniformely sampled projected parameter)

meanSpacing = sPlane(end) / length(sPlane);

disp(' ');
disp(['    Input sampling rate exponent [m * 10^exp]. Mean spacing between projected data points is ' num2str(meanSpacing) ' m']);
disp('    0 -> 1 m');
disp('   -1 -> 1 dm');
disp('   -2 -> 1 cm');
disp('   -3 -> 1 mm');
disp('   -4 -> 0.1 mm');
disp('   -5 -> 0.01 mm');
disp('   -6 -> 1 ?m');
disp('   -7 -> 0.1 ?m');
disp('   -8 -> 0.01 ?m');

% NOTE: these predefined sampling rates are necessary to be able to compare
% and average multiple plots in the following

sample = 100;
while not(-8 <= sample && sample <= 0)
    sample = round(input(' > '));
end

sample = 10^sample;

disp(' ');
disp(['   Sampling rate [m] is ' num2str(sample)]);

disp(' ');
disp(['    Sampling frequency = ' num2str(1/sample) ' m-1']);

t = 0:sample:sPlane(end); t=t';  % this is the regularly spaced abscissa for the regularly resampled distance to plane "dist"

tLS = sPlane;   % this is the irregularly spaced abscissa for the irregular non-resampled distance to plane "distLS" - "LS" refers to the Lomb-Scargle power spectrum algorithm

distLS = disti; % this is the non-resampled distance to plane with abscissa "tLS" - "LS" refers to the Lomb-Scargle power spectrum algorithm

% ====================
% SIMPLE INTERPOLATION
% ====================
%
% regularInterpDist = interp1(sPlane,disti,t);
% dist = regularInterpDist;

% ====================
% MOVING AVERAGE FILER
% ====================
%
movAvDist = zeros(size(t));

disp(' ');
disp('    Moving average filter...');

for i = 2:length(t)-1;
    indexes = find(sPlane >= t(i-1) & sPlane <= t(i+1));
    movAvDist(i) = mean(disti(indexes));
end

disp('    done.');

disp(' ');
disp('    Filling NaNs by interpolation...');

NaNindexes = find(isnan(movAvDist));

for i = 1:length(NaNindexes);
    
    j = NaNindexes(i);
    left = j-1;
    right = j+1;
    
    if j == 1, movAvDist(1) = 0; break, end
    if j == length(movAvDist), movAvDist(end) = 0; break, end
    
    while isnan(movAvDist(left));
        left = left-1;
    end
    
    while isnan(movAvDist(right));
        right = right+1;
    end
    
    movAvDist(j) = movAvDist(left)*(right-j)/(right-left) + movAvDist(right)*(j-left)/(right-left);
end

disp('    done.');

dist = movAvDist;

% TAPER WAS HERE - MOVED TO PERIODOGRAM

profileLenght = t(end);

% RMS roughness, center-line average roughness (Power et al., 1988), and histogram

meanDist = mean(dist);  % this should be zero
varDist = var(dist);
sigmaRMSdist = std(dist,1);
distX = linspace(-3*sigmaRMSdist,3*sigmaRMSdist,31);
normPDFdist = normpdf(distX,meanDist,sigmaRMSdist)/sigmaRMSdist; % PDF must be scaled by 1/sigma so that the integral is still 1
max(normPDFdist)
RaDist = sum(abs(dist))/length(dist);

disp([' ']);
disp(['Mean of distance from mean axis (should be zero)  = ' num2str(meanDist) ' m']);
disp(['              Variance of distance from mean axis = ' num2str(varDist) ' m^2']);
disp(['   RMS of distance (RMS roughness) from mean axis = ' num2str(sigmaRMSdist) ' m']);
disp(['                    Center-line average roughness = ' num2str(RaDist) ' m']);

% PLOT

title(['Blue = regularly resampled data - Red = irregular non-resampled data - profile is ' num2str(profileLenght) ' m long']);

plot1 = subplot(2,2,1); plot(t,dist,'b'); xlim([0,profileLenght]); xlabel('t [m]'); ylabel('dist [m]'); hold on; plot(tLS,distLS,'r'); hold off
plot2 = subplot(2,2,3); plot(t,dist,'b'); xlim([0,profileLenght]); xlabel('t [m]'); ylabel('dist [m]'); hold on; plot(tLS,distLS,'r'); axis equal; hold off
plot3 = subplot(2,2,[2,4]); hist(dist,distX); xlim([-3*sigmaRMSdist,3*sigmaRMSdist]); xlabel('dist [m]'); ylabel('frequency'); hold on; plot(distX,normPDFdist,'r','LineWidth',2); hold off

disp(' ');
disp([' -> done.']);
disp([' slip2trace = ' num2str(slip2trace)]);


end
%%

%% 7- spectral analysis -> fourier_analysis_1D.m
function [Amplitude,Freq,method,lowerF,upperF,H,P1,sigmaRMS,Rsquare,taperPercent] = SPECTRUM(t,dist,tLS,distLS);

disp(' ');
disp('7- Power spectrum analysis');
close(1); figure(1);

% already calculated at previous step, but needed for plotting and scaling
varDist = var(dist);
sigmaRMSdist = std(dist,1);

% Spectral analysis of rough profile represented by t,dist curve - t and dist
% must be column vectors in metres. Input can be measured data with regular
% sampling interval or can be generated by script sythetic_selfaffine_profile.m

% Andrea Bistacchi 8/9/2009 - last review 14/12/2010
% Lomb-Scargle algorithm added on 28/7/2010 - DEACTIVATED

% % needed for Lomb-Scargle option - DEACTIVATED
% disp(' ');
% disp('Would you like to calculate the Lomb-Scargle periodogram?');
% disp('It takes forever! 0 = no, 1 = yes...');
% 
% LSflag = -1;
% while (LSflag < 0) | (LSflag > 1);
%     LSflag = round(input(' > '));
% end

N = size(dist,1);          % this is the number of data points - this is different from N in sythetic_selfaffine_profile.m scipt
profLength = t(end);        % length in metres of the profile
sampleInt = profLength/(N-1);	% sampling interval
sampleFreq = 1/sampleInt;	% sampling frequency
Nyquist = sampleFreq/2;	% Nyquist frequency

disp([' ']);
disp(['Length of profile  = ' num2str(profLength) ' m']);
disp(['Number of samples  = ' num2str(N)]);
disp(['Sampling interval  = ' num2str(sampleInt) ' m']);
disp(['Sampling frequency = ' num2str(sampleFreq) ' m^-^1']);
disp(['Nyquist frequency  = ' num2str(Nyquist) ' m^-^1']);

% (1) FFT - THIS SHOULD BE THE SAME AS THE PERIODGRAM - FOR DETAILS ON THIS
% AND ON THE FREQUENCY VECTOR SEE
% http://www.mathworks.com/support/tech-notes/1700/1702.html
% Scaling is pretty tricky and not well explained in Matlab documentation.
% Solved with Parseval's Theorem (see below).

% taper

disp(' ');
disp('   Input taper length [%]');
taperPercent = input('   > ')/100;

taperLength = t(end)*taperPercent;

disp(' ');
disp(['    Taper length = ' num2str(taperLength) ' m']);

% taper for regularly resampled data
leftTaperIndex = find(t <= taperLength);
rightTaperIndex = find(t >= t(end)-taperLength);
noTaperIndex = find(t > taperLength & t < t(end)-taperLength);

leftTaper = (1-cos(pi*t(leftTaperIndex)/taperLength))/2;
rightTaper = (1-cos(pi*(t(end)-t(rightTaperIndex))/taperLength))/2;
noTaper = ones(size(noTaperIndex));

taperFunction = [leftTaper ; noTaper ; rightTaper];

distTaper = dist.*taperFunction; % this is the resampled distance to plane with abscissa "t"

% periodogram

nfft = 2^nextpow2(N);	% use next highest power of 2 greater or equal to N in order to get an even FFT and a faster calculation
FFT = fft(distTaper,nfft)/N;	% lenfth(FFT) will be nfft+2; FFT(0) corrensponds to fPsd(0) = 0 and FFT(nfft+2) corrensponds to fPsd(nfft+2) = sampleFreq

NumUniquePts = 1+nfft/2;	% FFT is even (see above), hence only 1+nfft/2 points are unique, the rest are simmetrically redundant
FFT = FFT(1:NumUniquePts);	% discard the redundant simmetric points

modulusFFT = abs(FFT);
squaredModulusFFT = modulusFFT.^2;
nonscaledPSD = [squaredModulusFFT(1) ; squaredModulusFFT(2:end-1)*2 ; squaredModulusFFT(end)];  % multiply by 2 since we dropped half FFT; 0-frequency and Nyquist components are unique and must not be multiplied

fPsd = linspace(0,Nyquist,NumUniquePts); fPsd = fPsd';  % evenly spaced frequency column vector with NumUniquePts points

% Parseval's Theorem says that the area under the spectrum must equal
% the variance. For a spectrum S(f) at intervals of freq df the area is
% trapz(S)*df. Therefore, if the original signal was y(t), the scaling is
% by var(y)/(trapz(S)*df).

integrNonscaledPSD = trapz(fPsd(2:end),nonscaledPSD(2:end));   % this integral is calculated between 1/profile lenght frequency (approximately index = 2) and Nyquist frequency (index = end)
disp(['integrNonscaledPSD = ' num2str(integrNonscaledPSD)])
scalingPSD = varDist/integrNonscaledPSD;
PSD = nonscaledPSD * scalingPSD;
integrPSD = trapz(fPsd(2:end),PSD(2:end));
disp(['scalingPSD = ' num2str(scalingPSD)])
disp(['integrPSD scaled = ' num2str(integrPSD)])
disp(['varDist = ' num2str(varDist)])
disp(' ')

%PSD = nonscaledPSD * 2 * profLength; % this should be like in Candela's papers and it is almost equivalent to te above - it is explained in Mathematica notebook

% (3) MTM METHOD (Signal Processing Toolbox required)

[nonscaledPMTM,nonscaledPMTMconfidence,fPmtm] = pmtm(dist,[],nfft,sampleFreq);  % here the scaling is almost OK , but not perfect and the section below is needed

integrNonscaledPMTM = trapz(fPmtm(2:end),nonscaledPMTM(2:end));   % this integral is calculated between 1/profile lenght frequency (approximately index = 2) and Nyquist frequency (index = end)
disp(['integrNonscaledPMTM = ' num2str(integrNonscaledPMTM)])
scalingPMTM = varDist/integrNonscaledPMTM;
PMTM = nonscaledPMTM * scalingPMTM;
PMTMconfidence = nonscaledPMTMconfidence * scalingPMTM;
integrPMTM = trapz(fPmtm(2:end),PMTM(2:end));
disp(['scalingPMTM = ' num2str(scalingPMTM)])
disp(['integrPMTM scaled = ' num2str(integrPMTM)])
disp(['varDist = ' num2str(varDist)])
disp(' ')

% (4) WELCH METHOD (Signal Processing Toolbox required)

[nonscaledPWELCH,fWelch] = pwelch(dist,[],[],[],sampleFreq);  % here the scaling is almost OK , but not perfect and the section below is needed

integrNonscaledPWELCH = trapz(fWelch(1:end),nonscaledPWELCH(1:end));
disp(['integrNonscaledPWELCH = ' num2str(integrNonscaledPWELCH)])
scalingPWELCH = scalingPMTM;   % use scalingPMTM because Welch works on a smaller frequency interval
PWELCH = nonscaledPWELCH * scalingPWELCH;
integrPWELCH = trapz(fWelch(1:end),PWELCH(1:end));
disp(['scalingPWELCH = ' num2str(scalingPWELCH)])
disp(['integrPWELCH scaled = ' num2str(integrPWELCH)])
disp(['varDist = ' num2str(varDist)])
disp(' ')

% % (5) LOMB-SCARGLE METHOD (Press et al., 2001, translated in Matlab by C.
% % Saragiotis, Nov 2008) - optional - DEACTIVATED
% 
% if LSflag == 1,
%     [nonscaledPLS,fLS,alpha] = lomb(distLS,tLS);
%     PLS = nonscaledPLS * var(distLS) / trapz(fLS,nonscaledPLS); % scaling as above in method (1)
% end

% PLOTTING RESULTS

figure(1);
loglog(fPsd,PSD,'b'); hold on;
%loglog(fPsd,PSDlengthScaled,'.c'); hold on;
%loglog(fPsd,nonscaledPSD,'.k'); hold on;
loglog(fPmtm,PMTM,'r'); hold on;
loglog(fPmtm,PMTMconfidence,'--r'); hold on;
loglog(fWelch,PWELCH,'g'); hold on;
%if LSflag == 1, loglog(fLS,PLS,'m'); hold on; end - DEACTIVATED
title('Blue: variance-scaled squared modulus (periodogram) - Red: MTM - Green: Welch');
xlabel('Log spatial frequency [m^-^1]');
ylabel('Log power amplitude density [m^3]');


% LOG-LOG LEAST SQUARES FIT TO GET HURST EXPONENT

% select cut-off frequencies for least squares fit
disp(['Input spectral window (1/profLength = ' num2str(1/profLength) ' [m^-^1] - Nyquist frequency = ' num2str(Nyquist) ' [m^-^1])']);
lowerF = input('lower cut-off frequency (left limit): ');
upperF = input('upper cut-off frequency (right limit): ');
disp(' ');

% log-log linear fit of Periodogram

disp('Periodogram:');

[slopeP,HP,P1P,sigmaRMSP,RsquareP,PSD,fPsd] = loglogfit(PSD,fPsd,upperF,lowerF);

% log-log linear fit of MTM

disp('MTM:');

[slopeM,HM,P1M,sigmaRMSM,RsquareM,PMTM,fPmtm] = loglogfit(PMTM,fPmtm,upperF,lowerF);

% log-log linear fit of Welch

disp('Welch:');

[slopeW,HW,P1W,sigmaRMSW,RsquareW,PWELCH,fWelch] = loglogfit(PWELCH,fWelch,upperF,lowerF);

% summary

loglog([lowerF upperF],P1P*[lowerF upperF].^-(1+2*HP),'b','LineWidth',2); hold all;
loglog([lowerF upperF],P1M*[lowerF upperF].^-(1+2*HM),'r','LineWidth',2); hold all;
loglog([lowerF upperF],P1W*[lowerF upperF].^-(1+2*HW),'color',[0 .7 0],'LineWidth',2); hold all;
%loglog([lowerF upperF],P1W*[lowerF upperF].^-(1+2*HW),'g','LineWidth',2); hold all;

disp(' ');
disp( '                 Periodogram    MTM            Welch');
disp(['         Slope = ' num2str(slopeP,'%012g') '   ' num2str(slopeM,'%012g') '   ' num2str(slopeW,'%012g')]);
disp(['  Hurst exp. H = ' num2str(HP,'%012g') '   ' num2str(HM,'%012g') '   ' num2str(HW,'%012g')]);
disp(['  Intercept P1 = ' num2str(P1P,'%012g') '   ' num2str(P1M,'%012g') '   ' num2str(P1W,'%012g')]);
disp(['      sigmaRMS = ' num2str(sigmaRMSP,'%012g') '   ' num2str(sigmaRMSM,'%012g') '   ' num2str(sigmaRMSW,'%012g') '  calculated over lower-upper cut-off frequency interval']);
disp(['       Rsquare = ' num2str(RsquareP,'%012g') '   ' num2str(RsquareM,'%012g') '   ' num2str(RsquareW,'%012g')]);
disp(' ');
disp(['       varDist = ' num2str(varDist,'%012g') '       varDist/lenght = ' num2str(varDist/profLength,'%012g') ]);
disp(['  sigmaRMSdist = ' num2str(sigmaRMSdist,'%012g') '  sigmaRMSdist/lenght = ' num2str(sigmaRMSdist/profLength,'%012g') ]);
disp(' ');
disp(['  profLength = ' num2str(profLength,'%012g') ]);
disp(' ');

% % plot a line with specified Hurst exponent and intercept - DEACTIVATED
% 
% disp('Plot a line with specified Hurst exponent and intercept');
% specifiedH = input('Hurst exponent H: ');
% specifiedP1 = input('Intercept P1: ');
% 
% specifiedCoefs = [-1-2*specifiedH log10(specifiedP1)];
% specifiedLine = polyval(specifiedCoefs,LOGf);
% 
% figure(1);
% loglog(10.^LOGf,10.^specifiedLine,'c'); hold on;
% title(['Blue, Red, Green: power spectrum data - Black: least squares interpolation with H = ' num2str(H) ' - P1 = '  num2str(P1) ' - Cyan: specified line with H = ' num2str(specifiedH) ' - P1 = '  num2str(specifiedP1)]);
% 
% RMSDelta = sqrt(var(specifiedLine-interpLine)); % don't know if this is correct, but gives an idea
% normalizedRMSDelta = RMSDelta/(max(LogAmplitude-interpLine)-min(LogAmplitude-interpLine)); % don't know if this is correct, but gives an idea
% disp(['RMSDelta (difference between interpolated and specified line) = ' num2str(RMSDelta)]);
% disp(['normalized RMSDelta = ' num2str(normalizedRMSDelta*100) ' %']);
% disp(' ');

% select an algorithm
disp('Which reults do you like to save?');
disp('1- Power Spectral Density as variance-scaled square of modulus of FFT');
disp('2- Power Spectral Density estimate via the Thomson multitaper method (MTM)');
disp('3- Power Spectral Density estimate via the Welch method');
%if LSflag == 1, disp('4- Power Spectral Density estimate via Lomb-Scargle algorithm'); end  % needed for Lomb-Scargle option - DEACTIVATED

method = 0;
% % needed for Lomb-Scargle option - DEACTIVATED
% if LSflag == 1,
%     while (method < 1) | (method > 4);
%         method = round(input(' > '));
%     end
% elseif LSflag == 0,
while (method < 1) | (method > 3);
    method = round(input(' > '));
end
% end - DEACTIVATED

if method == 1, Amplitude = PSD; Freq = fPsd; H = HP; P1 = P1P; sigmaRMS = sigmaRMSP; Rsquare = RsquareP;
elseif method == 2, Amplitude = PMTM; Freq = fPmtm; H = HM; P1 = P1M; sigmaRMS = sigmaRMSM; Rsquare = RsquareM;
elseif method == 3, Amplitude = PWELCH; Freq = fWelch; H = HW; P1 = P1W; sigmaRMS = sigmaRMSW; Rsquare = RsquareW;
%elseif method == 4, Amplitude = PLS; Freq = fLS; - DEACTIVATED
end

end
%%

%% 8- Save processed data
function SAVE(Amplitude,Freq,H,P1,sigmaRMS,slip2trace,monoDimension,planeNormal,slipVector,Min,Max,traceVector,sample,taperPercent,method,lowerF,upperF,filename,displ,profileLenght,Rsquare);

disp(' ');
disp('8- Save processed data to .mat file');

[file, path] = uiputfile('*.mat');
save([path file]);

disp(' ');
disp([' -> file ' file ' successfully saved.']);

end
%%

%% 10- Power spectrum summary plot (load last results to plot).
function SUMMARY(Amplitude,Freq,H,P1,sigmaRMS,slip2trace,plotNotes,displ,lowerF,upperF);

disp(' ');
disp('10- Power spectrum summary plot (load last results to plot).');

figure(2);
loglog(Freq,Amplitude,'Color',[.5 .5 .5]); hold all;
loglog([lowerF upperF],P1*[lowerF upperF].^-(1+2*H),'LineWidth',2); hold all;
title('Summary of power spectrum data');
xlabel('Log spatial frequency [m^-^1]');
ylabel('Log power amplitude density [m^3]');
legend show

figure(3);
plot(H,P1,'s','MarkerSize',12,'LineWidth',2); hold all;
title('Summary of power spectrum linear fit parameters');
xlabel('H: Hurst exponent [ ]');
ylabel('P1: Power amplitude density @ 1m^-^1 frequency [m^3]');
legend show

figure(11);
plot(H,sigmaRMS,'s','MarkerSize',12,'LineWidth',2); hold all;
title('Summary of power spectrum linear fit parameters');
xlabel('H: Hurst exponent [ ]');
ylabel('sigmaRMS roughness [m]');
legend show

figure(4);
plot(H,slip2trace,'s','MarkerSize',12,'LineWidth',2); hold all;
title('Summary of H vs. angle between measured fault trace and slip vector');
xlabel('H: Hurst exponent [ ]');
ylabel('Angle between measured fault trace and slip vector [°]');
legend show

figure(5);
plot(P1,slip2trace,'s','MarkerSize',12,'LineWidth',2); hold all;
title('Summary of P1 vs. angle between measured fault trace and slip vector');
xlabel('P1: Power amplitude density @ 1m^-^1 frequency [m^3]');
ylabel('Angle between measured fault trace and slip vector [°]');
legend show

figure(12);
plot(sigmaRMS,slip2trace,'s','MarkerSize',12,'LineWidth',2); hold all;
title('Summary of sigmaRMS vs. angle between measured fault trace and slip vector');
xlabel('sigmaRMS roughness [m]');
ylabel('Angle between measured fault trace and slip vector [°]');
legend show

figure(6);
plot(H,displ,'s','MarkerSize',12,'LineWidth',2); hold all;
title('Summary of H vs. displacement');
xlabel('H: Hurst exponent [ ]');
ylabel('Displacement [m]');
legend show

figure(7);
plot(P1,displ,'s','MarkerSize',12,'LineWidth',2); hold all;
title('Summary of P1 vs. displacement');
xlabel('P1: Power amplitude density @ 1m^-^1 frequency [m^3]');
ylabel('Displacement [m]');
legend show

figure(13);
plot(sigmaRMS,displ,'s','MarkerSize',12,'LineWidth',2); hold all;
title('Summary of sigmaRMS vs. displacement');
xlabel('sigmaRMS roughness [m]');
ylabel('Displacement [m]');
legend show

disp(' ');
disp(['Plot Notes:']);
disp(plotNotes);

end
%%

%% 11- Power spectrum summary plot of multiple previously saved files.
function MULTIPLOT;

disp(' ');
disp('11- Power spectrum summary plot of multiple previously saved files.');
disp(' ');
disp('Load multiple previously saved files');

plotNotes =[{'plot No'}, {'monoDimension'},...
    {'plane normal D/D'},{'slip vector P/T'},...
    {'profileLenght'},...
    {'displ'},...
    {'profile Min'}, {'profile Max'},...
    {'trace vector'}, {'slip2trace'}, {'sampling'}, {'taper'},...
    {'spectrum method'}, {'lower freq'}, {'upper freq'}, {'H'}, {'P1'}, {'sigmaRMS'}, {'Rsquare'},...
    {'file name'}];

[filename, pathname] = uigetfile('*.mat','MultiSelect','on');
nfiles = max(max(max(size(filename))))
filelist = cell(nfiles);

Hsummary = [];
P1summary = [];
sigmaRMSsummary = [];
slip2traceSummary = [];

for i = 1:nfiles;
    filelist(i) = cellstr(strcat(pathname,char(filename(i))));
end

disp(' ');
disp([' -> ' num2str(nfiles) ' files selected.']);
disp(' ');
disp(['List of plotted files']);

for i = 1:nfiles;
    
    clear monoDimension planeNormal slipVector profileLenght displ Min Max traceVector slip2trace sample taperPercent method lowerF upperF H P1 sigmaRMS filename Rsquare;
    
    sample = 0; % needed because of conflict with sample.m function
    displ = 0;  % needed for joints where displ is not defined
    profileLenght = -1;  % needed for backward compatibility
    
    load(char(filelist(i)));
    
    if exist('P1','var')==0, % needed for backward compatibility with files where P1 wasn't saved
        P1 = A10/(10^-(1+2*H));
    end
    
    fMax = 1e2; fMin = 1e-2; %  -> note fMax and fMin
    sigmaRMS = sqrt(P1 * (fMax^-(2*H) - fMin^-(2*H)) / -(2*H)); % needed for backward compatibility with files where sigmaRMS was wrong
    
    thisPlotNotes = [{num2str(i)}, {num2str(monoDimension)},...
        {dipdirstr(planeNormal)},{plungetrendstr(slipVector)},...
        {num2str(profileLenght)},...
        {num2str(displ)},...
        {num2str(Min)}, {num2str(Max)},...
        {plungetrendstr(traceVector)}, {num2str(slip2trace)}, {num2str(sample)}, {num2str(taperPercent)},...
        {num2str(method)}, {num2str(lowerF)}, {num2str(upperF)}, {num2str(H)}, {num2str(P1)}, {num2str(sigmaRMS)}, {num2str(Rsquare)},...
        {filename}];
    
    plotNotes = [plotNotes; thisPlotNotes];
    
    Hsummary = [Hsummary H];
    P1summary = [P1summary P1];
    sigmaRMSsummary = [sigmaRMSsummary sigmaRMS];
    slip2traceSummary = [slip2traceSummary slip2trace];
    
    figure(2);
    loglog(Freq,Amplitude); hold all;
    title('Summary of power spectrum data');
    xlabel('Log spatial frequency [m^-^1]');
    ylabel('Log power amplitude density [m^3]');
    legend show
    
    figure(19);
    loglog([lowerF upperF],P1*[lowerF upperF].^-(1+2*H),'LineWidth',2); hold all;
    title('Summary of power spectrum data power law best fit');
    xlabel('Log spatial frequency [m^-^1]');
    ylabel('Log power amplitude density [m^3]');
    legend show
    
    figure(20);
    loglog(Freq,Amplitude,'Color',[.5 .5 .5]); hold all;
    loglog([lowerF upperF],P1*[lowerF upperF].^-(1+2*H),'LineWidth',2); hold all;
    title('Summary of power spectrum data power law best fit');
    xlabel('Log spatial frequency [m^-^1]');
    ylabel('Log power amplitude density [m^3]');
    legend show
    
    figure(3);
    %     plot(H,P1,'s','MarkerSize',12,'LineWidth',2); hold all;
    scatter(H,P1,100,displ,'filled','MarkerEdgeColor','k'); hold all;
    axis([0.5 1 1e-6 1e-3]);
    title('Summary of power spectrum linear fit parameters');
    xlabel('H: Hurst exponent [ ]');
    ylabel('Log P_1: Power amplitude density prefactor [m^3]');
    legend show
    
    figure(11);
    scatter(H,sigmaRMS,100,displ,'filled','MarkerEdgeColor','k'); hold all;
    title('Summary of power spectrum linear fit parameters');
    axis([0.5 1 0 0.7]);
    xlabel('H: Hurst exponent [ ]');
    ylabel('\sigma_R_M_S roughness [m]');
    legend show
    
    figure(4);
    %    plot(H,slip2trace,'s','MarkerSize',12,'LineWidth',2); hold all;
    scatter(H,slip2trace,100,displ,'filled','MarkerEdgeColor','k'); hold all;
    axis([0.5 1 0 90]);
    title('Summary of H vs. angle between measured fault trace and slip vector');
    xlabel('H: Hurst exponent [ ]');
    ylabel('Angle between fault trace and slip vector [°]');
    legend show
    
    figure(5);
    %     plot(P1,slip2trace,'s','MarkerSize',12,'LineWidth',2); hold all;
    scatter(P1,slip2trace,100,displ,'filled','MarkerEdgeColor','k'); hold all;
    axis([1e-6 1e-3 0 90]);
    title('Summary of P1 vs. angle between measured fault trace and slip vector');
    xlabel('Log P_1: Power amplitude density prefactor [m^3]');
    ylabel('Angle between fault trace and slip vector [°]');
    legend show
    
    figure(12);
    %     plot(P1,slip2trace,'s','MarkerSize',12,'LineWidth',2); hold all;
    scatter(sigmaRMS,slip2trace,100,displ,'filled','MarkerEdgeColor','k'); hold all;
    axis([0 0.7 0 90]);
    title('Summary of \sigma_RMS vs. angle between measured fault trace and slip vector');
    xlabel('\sigma_R_M_S roughness [m]');
    ylabel('Angle between fault trace and slip vector [°]');
    legend show
    
    figure(6);
    %     plot(H,displ,'s','MarkerSize',12,'LineWidth',2); hold all;
    scatter(H,displ,100,slip2trace,'filled','MarkerEdgeColor','k'); hold all;
    title('Summary of H vs. displacement');
    xlabel('H: Hurst exponent [ ]');
    ylabel('Displacement [m]');
    legend show
    
    figure(7);
    %     plot(P1,displ,'s','MarkerSize',12,'LineWidth',2); hold all;
    scatter(P1,displ,100,slip2trace,'filled','MarkerEdgeColor','k'); hold all;
    title('Summary of P1 vs. displacement');
    xlabel('P1: Power amplitude density @ 1m^-^1 frequency [m^3]');
    ylabel('Displacement [m]');
    legend show
    
    figure(13);
    %     plot(P1,displ,'s','MarkerSize',12,'LineWidth',2); hold all;
    scatter(sigmaRMS,displ,100,slip2trace,'filled','MarkerEdgeColor','k'); hold all;
    title('Summary of sigmaRMS vs. displacement');
    xlabel('sigmaRMS roughness [m]');
    ylabel('Displacement [m]');
    legend show
    
    disp([num2str(i) ' -> ' char(filelist(i))]);
    
end

disp(' ');
disp(['Plot Notes:']);
disp(plotNotes);

% plot showing which profile belongs to a single fault surface
disp(' ');
disp('Would you like a plot showing which profile belongs');
disp('to a single fault surface? 0 = no, 1 = yes.');

multiplot = -1;
while (multiplot < 0) | (multiplot > 1);
    multiplot = round(input(' > '));
end

if multiplot ==1,
    % assign profile to fault surface
    disp(' ');
    disp('Enter id numbers of profiles belonging to a single');
    disp('fault surface in a row vector (with square brackets),');
    disp('<enter> to finish row, <0> to exit.');
    
    profileIDs = -1;
    plotsArray = {};
    faultNumber = 0;
    while 1,
        profileIDs = input(' > ');
        if profileIDs == 0, break, end
        faultNumber = faultNumber+1;
        plotsArray{faultNumber,2} = profileIDs;
    end
    
    for i = 1:faultNumber,
        % assign correnspondent fault name,H, P1, sigmaRMS, and slip3trace data
        disp(' ');
        disp('Enter name of fault composed by profiles:');
        disp(plotsArray{i,2});
        
        disp(' ');
        faultName = input('Fault name (input as text string) > ');
        plotsArray{i,1} = faultName;
        
        profileIDs = plotsArray{i,2};
        
        plotsArray{i,2} = [plotsArray{i,2}; Hsummary(profileIDs)];
        plotsArray{i,2} = [plotsArray{i,2}; P1summary(profileIDs)];
        plotsArray{i,2} = [plotsArray{i,2}; slip2traceSummary(profileIDs)];
        plotsArray{i,2} = [plotsArray{i,2}; sigmaRMSsummary(profileIDs)];
        
        % rearrange data by slip2trace angle
        
        plotsArray{i,2} = (sortrows(plotsArray{i,2}',4))';
    end
    
    maxPlot16 = 0;
    maxPlot17 = 0;
    maxPlot18 = 0;
    
    for i = 1:faultNumber,
        % find max for polar plots
        maxPlot16 = max([maxPlot16 max(plotsArray{i,2}(2,:))]);
        maxPlot17 = max([maxPlot17 max(plotsArray{i,2}(3,:))]);
        maxPlot18 = max([maxPlot18 max(plotsArray{i,2}(5,:))]);
    end
    
    figure(16); polar(0,maxPlot16,'.'); hold all;
    figure(17); polar(0,maxPlot17,'.'); hold all;
    figure(18); polar(0,maxPlot18,'.'); hold all;
        
    for i = 1:faultNumber,
        % re-plot figures 3, 4, 5, 11 and 12, as 8, 9, 10, 14 and 15
        figure(8);
        plot(plotsArray{i,2}(2,:),plotsArray{i,2}(3,:),'--s','MarkerSize',12,'LineWidth',2); hold all;
        title('Summary of power spectrum linear fit parameters');
        xlabel('H: Hurst exponent [ ]');
        ylabel('P1: Power amplitude density @ 1m^-^1 frequency [m^3]');
        legend show
        
        figure(9);
        plot(plotsArray{i,2}(2,:),plotsArray{i,2}(4,:),'--s','MarkerSize',12,'LineWidth',2); hold all;
        title('Summary of H vs. angle between measured fault trace and slip vector');
        xlabel('H: Hurst exponent [ ]');
        ylabel('Angle between measured fault trace and slip vector [°]');
        legend show
        
        figure(10);
        plot(plotsArray{i,2}(3,:),plotsArray{i,2}(4,:),'--s','MarkerSize',12,'LineWidth',2); hold all;
        title('Summary of P1 vs. angle between measured fault trace and slip vector');
        xlabel('P1: Power amplitude density @ 1m^-^1 frequency [m^3]');
        ylabel('Angle between measured fault trace and slip vector [°]');
        legend show
        
        figure(14);
        plot(plotsArray{i,2}(2,:),plotsArray{i,2}(5,:),'--s','MarkerSize',12,'LineWidth',2); hold all;
        title('Summary of power spectrum linear fit parameters');
        xlabel('H: Hurst exponent [ ]');
        ylabel('sigmaRMS roughness [m]');
        legend show
        
        figure(15);
        plot(plotsArray{i,2}(5,:),plotsArray{i,2}(4,:),'--s','MarkerSize',12,'LineWidth',2); hold all;
        title('Summary of sigmaRMS vs. angle between measured fault trace and slip vector');
        xlabel('sigmaRMS roughness [m]');
        ylabel('Angle between measured fault trace and slip vector [°]');
        legend show
        
        figure(16);
        %    polar(plotsArray{i,2}(4,:)./180*pi,plotsArray{i,2}(2,:),'--s','MarkerSize',12,'LineWidth',2); hold all;
        polar(plotsArray{i,2}(4,:)./180*pi,plotsArray{i,2}(2,:),'--s'); hold all;
        axis([0 maxPlot16 0 maxPlot16]);
        title('Summary of H (Hurst exponent [ ]) vs. angle between measured fault trace and slip vector');
        legend show
        
        figure(17);
        %    polar(plotsArray{i,2}(4,:)./180*pi,plotsArray{i,2}(3,:),'--s','MarkerSize',12,'LineWidth',2); hold all;
        polar(plotsArray{i,2}(4,:)./180*pi,plotsArray{i,2}(3,:),'--s'); hold all;
        axis([0 maxPlot17 0 maxPlot17]);
        title('Summary of P1 (Power amplitude density @ 1m^-^1 frequency [m^3]) vs. angle between measured fault trace and slip vector');
        legend show
        
        figure(18);
        %    polar(plotsArray{i,2}(4,:)./180*pi,plotsArray{i,2}(5,:),'--s','MarkerSize',12,'LineWidth',2); hold all;
        polar(plotsArray{i,2}(4,:)./180*pi,plotsArray{i,2}(5,:),'--s'); hold all;
        axis([0 maxPlot18 0 maxPlot18]);
        title('Summary of sigmaRMS roughness [m] vs. angle between measured fault trace and slip vector');
        legend show
        
    end
    
    figure(8); legend(plotsArray{:,1});
    figure(9); legend(plotsArray{:,1});
    figure(10); legend(plotsArray{:,1});
    figure(14); legend(plotsArray{:,1});
    figure(15); legend(plotsArray{:,1});
    figure(16); legend(plotsArray{:,1});
    figure(17); legend(plotsArray{:,1});
    figure(18); legend(plotsArray{:,1});
    
end

end
%%

%% Direction cosine vs. geologic notation
function [plunge,trend] = plungetrend(columnvector);
trend = atan2(columnvector(1,1),columnvector(2,1))*180/pi;  % IMPORTANT. DO NOT USE ATAN! IT IS LIMITED TO -PI/2 PI/2.
plunge = -asin(columnvector(3,1))*180/pi;

if isnan(trend), trend = 0; elseif trend < 0, trend = trend + 360; end

end

function text = plungetrendstr(columnvector);
trend = atan2(columnvector(1,1),columnvector(2,1))*180/pi;  % IMPORTANT. DO NOT USE ATAN! IT IS LIMITED TO -PI/2 PI/2.
plunge = -asin(columnvector(3,1))*180/pi;

if isnan(trend), trend = 0; elseif trend < 0, trend = trend + 360; end

text = [num2str(plunge) ' / ' num2str(trend)];

end

function [dip,dir] = dipdir(columnvector);
dir = 180 + atan2(columnvector(1,1),columnvector(2,1))*180/pi;  % IMPORTANT. DO NOT USE ATAN! IT IS LIMITED TO -PI/2 PI/2.
dip = 90 - (-asin(columnvector(3,1))*180/pi);

if dip > 90, dir = dir + 180; dip = 180 - dip; end

if isnan(dir), dir = 0; elseif dir < 0, dir = dir + 360; end

end

function text = dipdirstr(columnvector);
dir = 180 + atan2(columnvector(1,1),columnvector(2,1))*180/pi;  % IMPORTANT. DO NOT USE ATAN! IT IS LIMITED TO -PI/2 PI/2.
dip = 90 - (-asin(columnvector(3,1))*180/pi);

if dip > 90, dir = dir + 180; dip = 180 - dip; end

if isnan(dir), dir = 0; elseif dir < 0, dir = dir + 360; end

text = [num2str(dip) ' / ' num2str(dir)];

end

function [slipPlunge,slipTrend] = rake2slip(dip,dir,rake); % rake according to Aki & Richards, 1980, p. 106
slipPlunge   = asin(sin(-rake*pi/180)*sin(dip*pi/180))*180/pi;  % along dip component of rake
relSlipTrend = atan2(sin(-rake*pi/180)*cos(dip*pi/180),cos(-rake*pi/180))*180/pi;  % along strike component of rake
slipTrend    = dir -90 + relSlipTrend;  % "-90" since rake refers to right-hand-rule strike

if isnan(slipTrend), slipTrend = 0; elseif slipTrend < 0, slipTrend = slipTrend + 360; end

end

function text = rake2slipstr(dip,dir,rake); % rake according to Aki & Richards, 1980, p. 106
slipPlunge   = asin(sin(-rake*pi/180)*sin(dip*pi/180))*180/pi;  % along dip component of rake
relSlipTrend = atan2(sin(-rake*pi/180)*cos(dip*pi/180),cos(-rake*pi/180))*180/pi;  % along strike component of rake
slipTrend    = dir -90 + relSlipTrend;  % "-90" since rake refers to right-hand-rule strike

if isnan(slipTrend), slipTrend = 0; elseif slipTrend < 0, slipTrend = slipTrend + 360; end

text = [num2str(slipPlunge) ' / ' num2str(slipTrend)];

end

function columnvector = lineation(plunge,trend);

end

function columnvector = pole(dip,dir);

end
%%

%%
% all code below this line comes from function
% lomb.m by C. Saragiotis, Nov 2008, available
% on Mathworks web site
%
% just copied from original m-file and pasted here
%
% it is called in section 7- spectral analysis
%%

function [P,f,alpha] = lomb(x,t,varargin)
% LOMB caculates the Lomb normalized periodogram (aka Lomb-Scargle, Gauss-Vanicek or Least-Squares spectrum) of a vector x with coordinates in t.
%   The coordinates need not be equally spaced. In fact if they are, it is
%   probably preferable to use PWELCH or SPECTRUM. For more details on the
%   Lomb normalized periodogram, see the excellent section 13.8 in [1], pp.
%   569-577.
%
%   This code is a transcription of the Fortran subroutine period in [1]
%   (pp.572-573). The calculation of the Lomb normalized periodogram is in
%   general a slow procedure. Matlab's characteristics have been taken into
%   account in order to make it fast for Matlab but still it is quite slow.
%   For a faster (but not exact) version see FASTLOMB.
%
%   If the 'fig' flag is on, the periodogram is created, along with some
%   default statistical significance levels: .001, .005, .01, .05, .1, .5.
%   If the user wants more significance levels, they can give them as input
%   to the function. Those will be red.
%
% SYNTAX
%   [P,f,alpha] = lomb(x,t,fig,hifac,ofac,a);
%   [P,f,alpha] = lomb(x,t,fig,hifac,ofac);
%   [P,f,alpha] = lomb(x,t,fig,hifac);
%   [P,f,alpha] = lomb(x,t,fig);
%   [P,f,alpha] = lomb(x,t);
%
% INPUTS
%   x:     the vector whose spectrum is wanted.
%   t:     coordinates of x (should have the same length).
%   fig:   if 0 (default), no figure is created.
%   hifac: the maximum frequency returned is
%               (hifac) x (average Nyquist frequency)
%          See "THE hifac PARAMETER" in the NOTES section below for more
%          details on the hifac parameter. Default is 1 (i.e. max frequency
%          is the Nyquist frequency)
%   ofac:  oversampling factor. Typically it should be 4 or larger. Default
%          is 4. See "INTERPRETATION AND SELECTION OF THE ofac PARAMETER"
%          in the NOTES section below for more details on the choice of
%          ofac parameter.
%   a:     additional significance levels to be drawn on the figure.
%
% OUTPUTS
%   P:     the Lomb normalized periodogram
%   f:     respective frequencies
%   alpha: statistical significance for each value of P
%
% NOTES
% A. INTERPRETATION AND SELECTION OF THE ofac PARAMETER [1]
%    The lowest independent frequency f to be examined is the inverse of
%    the span of the input data,
%               1/(tmax-tmin)=1/T.
%    This is the frequency such that the data can include one complete
%    cycle. In an FFT method, higher independent frequencies would be
%    integer multiples of 1/T . Because we are interested in the
%    statistical significance of ANY peak that may occur, however, we
%    should over-sample more finely than at interval 1/T, so that sample
%    points lie close to the top of ANY peak. This oversampling parameter
%    is the ofac. A value ofac >~4 might be typical in use.
%
% B. THE hifac PARAMETER [1]
%    Let fhi be the highest frequency of interest. One way to choose fhi is
%    to compare it with the Nyquist frequency, fc, which we would obtain, if
%    the N data points were evenly spaced over the same span T, that is
%               fc = N/(2T).
%    The input parameter hifac, is defined as fhi/fc. In other words, hifac
%    shows how higher (or lower) that the fc we want to go.
%
% REFERENCES
%   [1] W.H. Press, S.A. Teukolsky, W.T. Vetterling and B.P. Flannery,
%       "Numerical recipes in Fortran 77: the art of scientific computing",
%       2nd ed., vol. 1, Cambridge University Press, NY, USA, 2001.
%
% C. Saragiotis, Nov 2008


%% Inputs check and initialization
if nargin < 2, error('%s: there must be at least 2 inputs.',mfilename); end

[x,t,hifac,ofac,a_usr,f,fig] = init(x,t,varargin{:});
nf = length(f);

mx = mean(x);
x  = x-mx;
vx = var(x);
if vx==0, error('%s: x has zero variance',upper(mfilename)); end


%% Main

P = zeros(nf,1);
for i=1:nf
    wt  = 2*pi*f(i)*t;  % \omega t
    swt = sin(wt);
    cwt = cos(wt);
    
    %% Calculate \omega\tau and related quantities
    % I use some trigonometric identities to reduce the computations
    Ss2wt = 2*cwt.'*swt;            % \sum_t \sin(2\omega\t)
    Sc2wt = (cwt-swt).'*(cwt+swt);  % \sum_t \cos(2\omega\t)
    wtau  = 0.5*atan2(Ss2wt,Sc2wt);  %\omega\tau
    
    swtau = sin(wtau);
    cwtau = cos(wtau);
    
    % I use some trigonometric identities to reduce the computations
    swttau = swt*cwtau - cwt*swtau;  % \sin\omega(t-\tau))
    cwttau = cwt*cwtau + swt*swtau;  % \cos\omega(t-\tau))
    
    P(i) = ((x.'*cwttau)^2)/(cwttau.'*cwttau) + ((x.'*swttau)^2)/(swttau.'*swttau);
end
P = P/(2*vx);

%% Significance
M = 2*nf/ofac;
alpha = 1 - (1-exp(-P)).^M;              % statistical significance
alpha(alpha<0.1) = M*exp(-P(alpha<0.1)); % (to avoid round-off errors)

%% Figure
if fig
    figure
    styles = {':','-.','--'};
    
    a = [0.001 0.005 0.01 0.05 0.1 0.5];
    La = length(a);
    z = -log(1-(1-a).^(1/M));
    hold on;
    for i=1:La
        line([f(1),0.87*f(end)],[z(i),z(i)],'Color','k','LineStyle',styles{ceil(i*3/La)});
        text(0.9*f(end),z(i),strcat('\alpha = ',num2str(a(i))),'fontsize',8);
        %         lgd{i}=strcat('\alpha=',num2str(a(i)));
    end
    if ~isempty(a_usr)
        [tmp,ind] = intersect(a_usr,a);
        a_usr(ind)=[];
        La_usr = length(a_usr);
        z_usr  = -log(1-(1-a_usr).^(1/M));
        for i = 1:La_usr
            line([f(1),0.87*f(end)],[z_usr(i),z_usr(i)],'Color','r','LineStyle',styles{ceil(i*3/La_usr)});
            text(0.9*f(end),z_usr(i),strcat('\alpha = ',num2str(a_usr(i))),'fontsize',8);
            %         lgd{La+i}=strcat('\alpha=',num2str(a_usr(i)));
        end
        z = [z z_usr];
    end
    %     legend(lgd);
    plot(f,P,'k');
    title('Lomb-Scargle normalized periodogram')
    xlabel('f (Hz)'); ylabel('P(f)')
    xlim([0 f(end)]); ylim([0,1.2*max([z'; P])]);
    hold off;
end

end
%%

%% init (initialize)
function [x,t,hifac,ofac,a,f,fig] = init(x,t,varargin)
if nargin < 6, a = [];    % set default value for a
else           a = sort(varargin{4});
    a = a(:)';
end
if nargin < 5, ofac = 4;  % set default value for ofac
else           ofac = varargin{3};
end
if nargin < 4, hifac = 1; % set default value for hifac
else           hifac = varargin{2};
end
if nargin < 3, fig = 0;   % set default value for hifac
else           fig = varargin{1};
end

if isempty(ofac),  ofac  = 4; end
if isempty(hifac), hifac = 1; end
if isempty(fig),   fig   = 0; end

if ~isvector(x) ||~isvector(t),
    error('%s: inputs x and t must be vectors',mfilename);
else
    x = x(:); t = t(:);
    nt = length(t);
    if length(x)~=nt
        error('%s: Inputs x and t must have the same length',mfilename);
    end
end

[t,ind] = unique(t);    % discard double entries and sort t
x = x(ind);
if length(x)~=nt, disp(sprintf('WARNING %s: Double entries have been eliminated',mfilename)); end

T = t(end) - t(1);
nf = round(0.5*ofac*hifac*nt);
f = (1:nf)'/(T*ofac);
end

%%
% log-log linear fit of Periodogram

function [slope,H,P1,sigmaRMS,Rsquare,Amplitude,Freq] = loglogfit(Amplitude,Freq,upperF,lowerF);

%select frequency interval
[dunmmy,lowerFindex] = min(abs(Freq-lowerF));
lowerF = Freq(lowerFindex);
[dunmmy,upperFindex] = min(abs(Freq-upperF));
upperF = Freq(upperFindex);

disp(['lower cut-off frequency = ' num2str(lowerF)]);
disp(['upper cut-off frequency = ' num2str(upperF)]);
disp('');

% trim frequency and amplitude vectors

Freq = Freq(lowerFindex:upperFindex);
Amplitude = Amplitude(lowerFindex:upperFindex);

% logaritm of frequency and power amplitude vectors
LOGf = log10(Freq);
LogAmplitude = log10(Amplitude);

% least squares fit (backslash works with least squares for overdetermined systems)
A = [LOGf ones(size(LOGf))];
B = LogAmplitude;

interpCoefs = A\B;

interpLine = polyval(interpCoefs,LOGf);

ESS = sum((interpLine-mean(LogAmplitude)).^2); % Explained Sum of Squares
TSS = sum((LogAmplitude-mean(LogAmplitude)).^2); % Total Sum of Squares
RSS = sum((LogAmplitude-interpLine).^2); % Residual Sum of Squares

Rsquare = ESS/TSS; % Coefficient of determination

slope = interpCoefs(1);     % slope of power spectrum profile
H = -(interpCoefs(1)+1)/2;   % Hurst exponent -> interpCoefs(1) = -1 -2*H
P1 = 10^interpCoefs(2);   % P1 = intercept at frequency = 1
%fMax = 1e0; fMin = 1e-2; sigmaRMS = sqrt(P1 * (fMax^-(2*H) - fMin^-(2*H)) / -(2*H)); % sigmaRMS over the fMax-fMin interval
sigmaRMS = sqrt(P1 * (upperF^-(2*H) - lowerF^-(2*H)) / -(2*H)); % sigmaRMS over the fMax-fMin interval
end
