function [Robot, M, CoM, I] = Robot_dV()
%   Returns D-H table of parameters, mass, center of mass and
%   inertia matrices for robotic arm
%
%   Robot=[d v a alpha offset;
%          d v a alpha offset;
%          . . .   .   offset;
%          d v a alpha offset];
%
%   M = [m1, m2, m3, ...];
%
%   CoM=[xcm1 ycm1 zcm1;
%        xcm2 ycm2 zcm2;
%        ...];
%   CoM --> r ^i _i,Ci
%
%   I=[Il1 (3x3);
%      Il2 (3x3);
%      Il3 (3x3);
%      ...];

% D-H table of parameters 

syms q1 q2 q3 a2 real
Robot = [ 0     q1     0      pi/2     0;
          0     q2     a2     pi/2     0;
          q3    0      0      0        0];


% Construction of links 
% Link 1 
% L1 = 0.4;
% r1 = 0.045;
% 
% % Link 2
% L2 = 0.25;
% b2 = 0.1;
% h2 = 0.04;
% 
% % Link 3
% L5 = 0.1;
% b5 = 0.06;
% h5 = 0.03;
% 
% rho = 2750; 

% Mass matrix
% M = [pi*r1^2*L1*rho  b2*h2*L2*rho  b3*h3*L3*rho  b4*h4*L4*rho  b5*h5*L5*rho];

syms m1 m2 m3 real

M = [m1 m2 m3];

% Center of mass matrix
% lc1 = L1/2;
% lc2 = L2/2;
% lc3 = L3/2;
% lcx4 = a4 + h4/2; 
% lcz4 = 0.05;
% lc5 = d5 + 0.2;

% CoM = [0        -lc1    0;
       % -lc2     0       0;
       % -lc3     0       0;
       % -lcx4    0       -lcz4;
       % 0        0       -lc5];


syms lc1 lc2 lc3 real

CoM = [0  -lc1     0;
       -lc2  0     0;
       0     0  -lc3];


% Moments of inertia

% I = {[M(1)/12 * (3*r1^2+L1^2) 0 0; 0 M(1)/12 * (3*r1^2+L1^2) 0; 0 0 M(1)/2 * r1^2]
%      [M(1)/12 * (b2^2+L2^2) 0 0; 0 M(1)/12 * (L2^2+h2^2) 0; 0 0 M(1)/12 * (L2^2+b2^2)]
%      [M(1)/12 * (b3^2+h3^2) 0 0; 0 M(1)/12 * (L3^2+h3^2) 0; 0 0 M(1)/12 * (L3^2+b3^2)]
%      [M(1)/12 * (b4^2+h4^2) 0 0; 0 M(1)/12 * (L4^2+b4^2) 0; 0 0 M(1)/12 * (L4^2+h4^2)]
%      [M(1)/12 * (L5^2+h5^2) 0 0; 0 M(1)/12 * (L5^2+b5^2) 0; 0 0 M(1)/12 * (h5^2+b5^2)]};


syms Ix [3, 1] real
syms Ixy [3, 1] real
syms Ixz [3, 1] real
syms Iy [3, 1] real
syms Iyz [3, 1] real
syms Iz [3, 1] real


% I = {[Ix1 Ixy1 Ixz1; Ixy1 Iy1 Iyz1; Ixz1 Iyz1 Iz1]
%      [Ix2 Ixy2 Ixz2; Ixy2 Iy2 Iyz2; Ixz2 Iyz2 Iz2]
%      [Ix3 Ixy3 Ixz3; Ixy3 Iy3 Iyz3; Ixz3 Iyz3 Iz3]};

I = {[Ix1 0 0; 0 Iy1 0; 0 0 Iz1]
     [Ix2 0 0; 0 Iy2 0; 0 0 Iz2]
     [Ix3 0 0; 0 Iy3 0; 0 0 Iz3]};

end

