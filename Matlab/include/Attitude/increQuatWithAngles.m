function [ qkp1 ] = increQuatWithAngles( qk, theta )
% Increment Quaternion vector with angular increments
% 
% @arg
% 
% @return
% 
% @author: Matt Marti
% @date: 2019-03-04

% Obtain norm of input
thetax = theta(1);
thetay = theta(2);
thetaz = theta(3);
theta = sqrt(thetax^2 + thetay^2 + thetaz^2);

% Sine and Cosine things
halftheta = 0.5*theta;
s = insSinc(halftheta);
c = 2*(cos(halftheta) - 1);

% Quaternion update
A = [ c,         s*thetaz, -s*thetay, s*thetax ; ...
     -s*thetaz,  c,         s*thetax, s*thetay ; ...
      s*thetay, -s*thetax,  c,        s*thetaz ; ...
     -s*thetax, -s*thetay, -s*thetaz, c        ];
qkp1 = qk +  0.5*A*qk;

end

