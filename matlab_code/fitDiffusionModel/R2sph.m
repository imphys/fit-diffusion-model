function [alpha1,alpha2,alpha3] = R2sph(EV)
% angles are pitch, roll, yaw, see:
% http://planning.cs.uiuc.edu/node103.html

EV(:,3) = cross(EV(:,1),EV(:,2)); %right-handed coordinate-system

alpha2 = asin(EV(1,3));

if abs(abs(alpha2)-pi/2)>eps
    alpha1 = atan2( -EV(2,3), EV(3,3));
    alpha3 = atan2(-EV(1,2), EV(1,1));
else
    alpha3 = 0;
    alpha1 = atan2(EV(3,2), EV(2,2));
end

end