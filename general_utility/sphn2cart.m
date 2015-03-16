function [x,y,z] = sphn2cart(r,p,t)
%SPHN2CART Transform spherical to Cartesian coordinates.
%   [X,Y,Z] = SPHN2CART(r,p,t) transforms corresponding elements of
%   data stored in spherical coordinates (azimuth R, theta T,
%   radius R) to Cartesian coordinates X,Y,Z.  The arrays T, P, and
%   R must be the same size (or any of them can be scalar).  T and
%   P must be in radians.
%
%   P is the counterclockwise angle in the xy plane measured from the
%   positive x axis. T is the angle from the z axis.
%
%   Class support for inputs T,P,R:
%      float: double, single
%
%   See also CART2SPHN

z = r .* cos(t);
rsint = r .* sin(t);
x = rsint .* cos(p);
y = rsint .* sin(p);
