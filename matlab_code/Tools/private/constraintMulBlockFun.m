function [c, ceq, Gc, Gceq] = constraintMulBlockFun(x, constraints , npar )
% [x] = constraintMulBlockFun(x, constraints , npar )
% 'Nonlinear' constraints multiplication function. 
% Handles all constraints of constrained optimzations.
% Constraints = scalar structure with fields:
%     Aineq, bineq, Aeq, beq, lb, ub, nonlcon
% which specify the constraints of a single voxel (see fmincon for meaning of each).
% Created by Dirk Poot, Erasmus MC, 21-10-2011;
x = reshape(x, npar, [] );
c = bsxfun(@minus, constraints.Aineq* x, constraints.bineq);

if isempty( constraints.Aeq )
    ceq = [];
else
    ceq = bsxfun(@minus, constraints.Aeq* x, constraints.beq);
end;
    
if ~isempty(constraints.nonlcon)
    [cn, ceqn, Gcn, Gceqn] = constraints.nonlcon( x );
    c = [c(:);cn(:)];
    ceq = [ceq(:);ceqn(:)];
end;
if nargout>2
    rp = repmat( {constraints.Aineq'} ,1,size(x,2));
    Gc = blkdiag(rp{:});
    rpeq = repmat({constraints.Aeq'},1,size(x,2));
    Gceq = blkdiag(rpeq{:});
    if ~isempty(constraints.nonlcon)
        Gc = [Gc Gcn];
        Gceq = [Gceq Gceqn];
    end;
end;
