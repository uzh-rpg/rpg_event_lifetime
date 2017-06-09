function r = iscolinear(p1,p2,p3)
%
% r = iscolinear(p1,p2,p3)
%
% [Input]:
% p1,p2,p3:     2D points (x,y)
%
% [Output]:
% r:            {0,1} boolean if points are colinear
%               r = 0 if not colinear
%               r = 1 if colinear
%
% This function decides if the three 2D points p1, p2, p3 = (xi,yi) are 
% colinear. If the 3 points are colinear, then the norm of the cross
% product of the two vectors p12 = p2-p1 and p13 = p3-p1 is equal to zero.

if(length(p1) ~= 2 || length(p2) ~= 2 || length(p3) ~= 2)
    error('Points must have dimension 2')
end

% extend to homogenous coordinates for cross product 
p1(3) = 1; p2(3) = 1; p3(3) = 1;

% return 1 if colinear (norm is bigger that eps)
r =  norm(cross(p2-p1, p3-p1)) < eps;

end
