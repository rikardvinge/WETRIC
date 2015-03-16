function [r,p,t] = cart2sphn(x,y,z)
%CART2SPHN Transform Cartesian to spherical coordinates.
%   [R,P,T] = CART2SPHN(X,Y,Z) transforms corresponding elements of
%   data stored in Cartesian coordinates X,Y,Z to spherical
%   coordinates (azimuth R, theta T, and radius R).  The arrays
%   X,Y, and Z must be the same size (or any of them can be scalar).
%   T and P are returned in radians.
%
%   P is the counterclockwise angle in the xy plane measured from the
%   positive x axis.  T is the angle from the z-axis.
%
%   Class support for inputs X,Y,Z:
%      float: double, single
%
%   See also SPHN2CART

hypotxy = hypot(x,y);
r = hypot(hypotxy,z);
t = acos( z./ r);
p = atan2(y,x);
