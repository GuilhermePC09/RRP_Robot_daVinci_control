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


syms m1 m2 m3 real

M = [m1 m2 m3];

% Center of mass matrix

syms lc1 lc2 lc3 real

CoM = [0  -lc1     0;
       -lc2  0     0;
       0     0  -lc3];


% Moments of inertia

syms Ix [3, 1] real
syms Ixy [3, 1] real
syms Ixz [3, 1] real
syms Iy [3, 1] real
syms Iyz [3, 1] real
syms Iz [3, 1] real


I = {[Ix1 Ixy1 Ixz1; Ixy1 Iy1 Iyz1; Ixz1 Iyz1 Iz1]
     [Ix2 Ixy2 Ixz2; Ixy2 Iy2 Iyz2; Ixz2 Iyz2 Iz2]
     [Ix3 Ixy3 Ixz3; Ixy3 Iy3 Iyz3; Ixz3 Iyz3 Iz3]};

% I = {[Ix1 0 0; 0 Iy1 0; 0 0 Iz1]
%      [Ix2 0 0; 0 Iy2 0; 0 0 Iz2]
%      [Ix3 0 0; 0 Iy3 0; 0 0 Iz3]};

end

