%% Discrete Elastic Rods (Bergou et al, SIGGRAPH, 2010) implementation
% Khalid Jawed, khalidjm@seas.ucla.edu
% May 2016
% A very fast C++ implementation can be found at
% http://www.cs.columbia.edu/cg/elastic_coiling/
% This MATLAB code is for demonstration only, and is very slow

%%
clear all;
clc;

FONT = 'Arial';
FONTSIZE = 12;
pWidth = 4; % inches
pHeight = 3;
colpos = [247 148 30;0 166 81; 237 28 36; 0 174 239; 120 254 50]/255; % colors

%%
fprintf('Discrete elastic rods\n');

%%
global m EI EA GJ dt x0 u nv ne d1 d2 m1 m2 x refLen voronoiRefLen tangent
global xCons nCons uUncons garr ctime RodLength r0
global consInd unconsInd
global ScaleSolver theta0 refTwist kappaBar
global tol maximum_iter
global Y G nu

%% Inputs
% number of vertices
nv = 50;

% Time step
dt = 1e-2;

% Rod Length
RodLength = 0.15;

% Density
rho = 1180;

% Cross-sectional radius of rod
r0 = 0.0032;

% Young's modulus
Y = 1.3e6;

% Poisson ratio
nu = 0.5;

% Shear modulus
G = Y/(2.0*(1.0+nu));

% gravity
% FOR DEMONSTRATION: let's set the gravity along both y and z directions so
% that the rod takes a 3D shape
g = [0, 0, -9.81];

% Tolerance on force function. This is multiplied by ScaleSolver so that we
% do not have to update it based on edge length and time step size
tol = 1e-7;

% Maximum number of iterations in Newton Solver
maximum_iter = 100;

% Total simulation time (it exits after t=totalTime)
totalTime = 10;

% Indicate whether images should be saved
saveImage = 0;

% How often the plot should be saved? (Set plotStep to 1 to show each plot)
plotStep = 1;

%% Utility quantities
ne = nv - 1;
EI = Y * pi * r0^4 / 4;
EA = Y * pi * r0^2;
GJ = G * pi * r0^4/2;
dm = pi * r0^2 * RodLength * rho / ne;

%% Geometry of the rod
% This section will have to updated based on the rod geometry
nodes = zeros(nv, 3);
for c=1:nv
    % nodes(c, 1) = - RodLength / 2 + (c-1) * RodLength / ne ;
    nodes(c, 1) = RodLength / 2 *sin( pi/nv * c - pi /2 );
    nodes(c, 3) = RodLength / 2 *cos( pi/nv * c - pi /2 );
end


%% Multiplier for Force & Jacobian
ScaleSolver = dm *(RodLength/ne) /dt^2;

%% Compute Mass
m = zeros(3*nv+ne, 1);
for c = 1:nv
    if c == 1
        m( 4 * (c-1) + 1 : 4 * (c-1) + 3) = dm/2;
    elseif c == nv
        m( 4 * (c-1) + 1 : 4 * (c-1) + 3) = dm/2;
    else
        m( 4 * (c-1) + 1 : 4 * (c-1) + 3) = dm;
    end
end
for c = 1:ne
    m( 4 * c ) = dm/2 * r0^2; % I = 1/2 m r ^ 2
end

%% gravity
garr = zeros(3*nv+ne, 1);
for c=1:nv
    garr( 4 * (c-1) + 1 : 4 * (c-1) + 3) = g;
end

%% Reference length and Voronoi length
refLen = zeros(ne, 1);
for c=1:ne
    dx = nodes(c+1, :) - nodes(c, :);
    refLen(c) = norm(dx);
end
voronoiRefLen = zeros(nv, 1);
for c=1:nv
    if c==1
        voronoiRefLen(c) = 0.5 * refLen(c);
    elseif c==nv
        voronoiRefLen(c) = 0.5 * refLen(c-1);
    else
        voronoiRefLen(c) = 0.5 * (refLen(c-1) + refLen(c));
    end
end

%% Reference director & material director
d1 = zeros(ne, 3); % reference director, u (or d1)
d2 = zeros(ne, 3); % reference director, v (or d2)
tangent = zeros(ne, 3); % tangent
for c=1:ne
    dx = nodes(c+1,:) - nodes(c,:);
    tangent(c,:) = dx / norm(dx);
end

% Figure out a good choice for d1(1)
t0 = tangent(1,:);
t1 = [0 0 -1];
d1Tmp = cross(t0, t1);
if (abs(d1Tmp) < 1.0e-6)
    t1 = [0 1 0];
    d1Tmp = cross(t0, t1);
end
d1(1,:) = d1Tmp;

d1_l = d1(1,:);
d2(1,:) = cross(t0, d1_l);
for c=2:ne
    t1 = tangent(c,:);
    d1_l = parallel_transport(d1_l, t0, t1);
    d1_l = (d1_l - dot(d1_l, t1) * t1);
    d1_l = d1_l / norm(d1_l);
    d1(c,:) = d1_l;
    d2_l = cross(t1, d1_l);
    d2(c,:) = d2_l;
    t0 = t1;
end

%% Initial
x0 = zeros(3*nv + ne, 1);
for c=1:nv
    x0( 4 * (c-1) + 1) = nodes(c,1);
    x0( 4 * (c-1) + 2) = nodes(c,2);
    x0( 4 * (c-1) + 3) = nodes(c,3);
end
x0(4:4:end) = 0; % theta
x = x0;

%% Constrained dofs
dummyInd = 1:length(x);
consInd = [1:7 (4*(nv-2)+1):4*nv-1]; % constrained dof
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

xCons = x(consInd); % first two nodes and one edge angle
nCons = length(xCons);
u = (x - x0) / dt;
uUncons = u(unconsInd); % unconstrained dof
leftLength = abs(x0(1) - x0(5));
rightLength = abs(x0(4*(nv-2)+1) - x0(4*(nv-1)+1));

%% Compute material director
theta0 = zeros(ne,1);
[m1, m2] = computeMaterialDirectors(d1, d2, theta0);
refTwist = zeros(ne,1);

%% Natural curvature computation
kappaBar = getkappa( x, m1, m2 );

%% Create director to save image
imageDirectory = date;
if (saveImage~=0)
    mkdir(imageDirectory);
end

%% Time marching
Nsteps = round(totalTime/dt); % number of time steps

ctime = 0;

endZ = zeros(Nsteps,1);

for timeStep=1: Nsteps
    fprintf('t=%.2f\n', ctime);

    middleZ = 0.5*(x(floor(nv/2)*4+3) + x(floor(nv/2)*4-1));
    
    xUncons = x(unconsInd) + uUncons * dt; % Guess
    [xUncons, error] = objfun(xUncons);
    
    zref = -RodLength + 0.1*sin(ctime);
    xCons = xCons;
    
    omega1 = -sin(10*ctime);
    kp = 10;
    % omega1 = (middleZ- zref) * kp * 0;
    xCons(5) = xCons(1) + leftLength * cos(omega1 * ctime);
    xCons(7) = xCons(3) + leftLength * sin(omega1 * ctime);
    
    omega2 = -omega1;
    xCons(8) = xCons(12) + rightLength * cos(omega1 * ctime + pi);
    xCons(10) = xCons(14) + rightLength * sin(omega2 * ctime + pi);


    x(consInd) = xCons;
    
    x(unconsInd) = xUncons;
    u = (x - x0) / dt; % velocity
    ctime = ctime + dt; % current time
    uUncons = u(unconsInd);
    
    % Theta
    theta0 = x(4:4:end);
    
    % Update material and reference directors
    tangent = computeTangent(x);
    [d1, d2] = computeTimeParallel(d1, x0, x);
    refTwist = getRefTwist(d1, tangent, refTwist); % Compute reference twist
    [m1, m2] = computeMaterialDirectors(d1, d2, theta0);
    
    % Update x0
    x0 = x;
    
    if (mod(timeStep, plotStep) ==0)
        plotrod(x, d1, d2, m1, m2, ctime, saveImage, imageDirectory);
    end
    
    
    obs(timeStep) = middleZ;
    tar(timeStep) = zref;
end

h1 = figure(2);
plot( (1: Nsteps) * dt, [obs; tar], 'o-', 'Color', colpos(1,:));
xlabel('Time, t [s]','Fontname', FONT,'FontSize',FONTSIZE);
ylabel('z-coord of last node, \delta_z [m]','Fontname', FONT,'FontSize',FONTSIZE);
box on
set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
    'PaperSize', [pWidth pHeight]);
saveas(h1, 'Height_Time.pdf');
