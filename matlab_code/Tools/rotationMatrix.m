function [out, doutdangl] = rotationMatrix(angles)
% [out, doutdangl] = rotationMatrix(angles)
% Computes the rotation matrix from the column vector angles.
% The second output is a matrix with derivatives (n x n x (n*(n-1)/2)),
% where the third dimension is specifies the angle to which the derivative is taken.
%
% Created by Dirk Poot, University of Antwerp
% 7-1-2009

n = ceil(sqrt(size(angles,1)*2));
if n*(n-1)/2~=size(angles,1)
    error('A correct number of angles should be provided ( n*(n-1)/2  for dimension n).\nClosest n = %d, with expected number of angles of %d', n, n*(n-1)/2);
end;

out = 1;
indx = 1;
if nargout>1
    doutdangl = zeros(0,1);
end
for k=2:n
    % increase size of rotation. 
    out = [out zeros(k-1,1);zeros(1,k-1) 1]; %#ok: increase = beneficial for speed.
    if nargout>1
        doutdangl = reshape([reshape([doutdangl zeros(size(doutdangl,1),1)],k-1,[]);zeros(1,(indx-1)*k)],[],k);
    end;
    for k2=1:k-1
        % apply rotation in elementary direction (could do sparse, but then
        % I think it is better when several elementary rotations are combined)
        base = eye(k);
        base([k2 k],[k2 k]) = [cos(angles(indx)) sin(angles(indx));-sin(angles(indx)) cos(angles(indx))];
        if nargout>1
            dbase = zeros(k);
            dbase([k2 k],[k2 k]) = [-sin(angles(indx)) cos(angles(indx));-cos(angles(indx)) -sin(angles(indx))];
            doutdangl = [doutdangl * base; out * dbase];
        end;
        indx = indx+1;
        out = out * base;
    end;
end;
if nargout>1
    doutdangl = permute(reshape(doutdangl,n,size(angles,1),n),[1 3 2]);
end;
