function [ Hx ] = PSFMulSelectedBlock_default( x, hLLifoBlk )
% [ Hx ] = PSFMulSelectedBlock_default( x, hLLifoBlk )
% Multiply x with the point spread function.
% INPUTS:
%  x : numimg x numvoxelsinblock x numvect
%  hLLifoBlk : output of PSFMulprepareBlock_default
% OUTPUTS:
% Hx : numimg x numvoxelsinblock
%      Hx(i,j) = sum_klm d2 L /d im_i(k) d im_i(l) * (d im_i(k)/d f_i(j) ) * (d im_i(l)/d f_i(m) ) * x(i,m)
%      with L = total log likelihood, im_i = project{i,1}( f_i )
%          with f = fun( theta ), f_i = f(i,:,:,:)
% Created by Dirk Poot, Erasmus MC, 19-10-2011

% Multiply each row of x with its corresponding psf
% Hx = x;
% for k = 1: size(x,1)
%     Hx(k,:) = x(k,:) * hLLifoBlk.psfBlk{k};
% end;
% However, we changed psfBlk to be size(psfBlk) = numim x numtr x numtr
% so we can do:
szx = [size(x) 1];
Hx = reshape( sum( bsxfun( @times, reshape(x,[szx(1),szx(2),1,szx(3)]) , hLLifoBlk.psfWithinBlock),2),szx);
 
function [ Hx ] = PSFMulLargeBlock_default( x, hLLifoBlk)
szx = [size(x) 1];
Hx = reshape( sum( bsxfun( @times, reshape(x,[szx(1),szx(2),1,szx(3)]), hLLifoBlk.psfOutsideBlock),2),szx(1), size(hLLifoBlk.psfOutsideBlock,3), szx(3) );


%%%%%%%%%%%%%%%%%%%  VALIDATION HELPER FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%

function validateDerivHess(fun, hessmulfun, x)
[f1] = fun(x);
[f2,g2] = fun(x);
[f,g, H] = fun(x);
goodresults = true;
if ~isequal(f1,f2,f) || ~isequal(g2,g)
    warning('Function value or gradient not equal if the number of requested outputs changes. This should not happen.'); % (When the difference is just roundoff errors it might occaisionally be acceptable).
    goodresults = false;
end;
% Using gradest from the 'DERIVESTsuite' by John D'Errico
% http://www.mathworks.com/matlabcentral/fileexchange/13490-automatic-numerical-differentiation
[gest , err] = gradest(fun, x);
gerr = (g(:)-gest(:));
if any( abs(gerr)./err(:) >2  &  abs(gerr)./abs(g(:)+gest(:))>1e-5 )
    disp('[ gradient,   numerical_gradient,  difference_in_gradient    absolute_difference_in_gradient/numerical_error]: ')
    fprintf( '   %15f  ,  %15f  , %15f  ,    %15f\n',[g(:) gest(:)  gerr  abs(gerr)./err(:)]');
    warning('error in gradient computation seems to be too large: ');
    goodresults = false;
end;
if ~isempty(hessmulfun)
    H = hessmulfun(H, eye(numel(x)));
end;
% Using jacobianest from the 'DERIVESTsuite' by John D'Errico
% http://www.mathworks.com/matlabcentral/fileexchange/13490-automatic-numerical-differentiation
[Hest , errHest] = jacobianest( @(x) reshape(shuffleoutputs(fun, 2, 2, {x}),[],1) , x);
Herr = (H-Hest);
if any( abs(Herr(:))./(errHest(:)+1e-4*max(errHest(:))) >2 & abs(Herr(:))./abs(H(:)+Hest(:)+1e-6*max(H(:)))>1e-4)
    warning('error in hessian computation seems to be too large.');
    goodresults = false;
end;
if goodresults
    disp('The gradient and hessian appear to be computed correctly');
end;
