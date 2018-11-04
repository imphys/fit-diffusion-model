function [LL, gL] = lineSearchThetaCostfun( theta0, deltatheta, alpha, opts )
% [LL, gL] = lineSearchThetaCostfun( theta0, deltatheta, alpha, opts )
% Evaluates the fit_MRI cost function at theta0 + alpha * deltatheta. 
%
% Created by Dirk Poot, Erasmus MC, 2-2-2012
curtheta = theta0 + alpha * deltatheta;
if nargout==1
    data = predictAllImagesfun( curtheta, opts );
    LL = projectAllAndLL(data, opts);
    LL = sum(LL(:));
    if opts.doRegularize
        LL = LL + opts.spatialRegularizer{1}(curtheta);
    end;
else
    [data , datagrd]= predictAllImagesfun( curtheta, opts , deltatheta);
    [LL, projgradifo, LLgrads] = projectAllAndLL(data, opts);
    LL= sum(LL);
    gL = sum( reshape( datagrd.*LLgrads ,[],1) );
    if opts.doRegularize
        [LLR, GLLR ] = opts.spatialRegularizer{1}(curtheta);
        LL = LL + LLR;
        gL = gL + deltatheta(:)'*GLLR(:);
    end;
end;