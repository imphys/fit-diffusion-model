function [laplaceFun, laplaceHessMulFun] = laplaceRegularizer_prepare( regularizerfuncs, explicitHessian, thetaselLrg, subindexsel, spacing , weights)
% [Hx] = laplaceRegularizerHessMul(hessinfo, x  , spacing , weights)
% Multiplies a vector with the hessian of the the log likelihood 
% regularization value for the Laplace regularizer in the spatial dimensions.
% INPUTS:
% hessinfo  : third output of laplaceRegularizer
% x         : vector or multicolumn matrix with parameter vectors that should be multiplied by 
%             the hessian
% spacing   : spacing(i): spatial distance between samples in dimension i
% weights   : weights(i): weight of each parameter i  ( = weigth of x(i, ..) )
%
% Created by Dirk Poot, Erasmus MC, 22-3-2011

if ~explicitHessian
    % use the default regularization function when hessian is not
    % explicitly requested.
    laplaceFun = [];
    laplaceHessMulFun = [];
    return;
end;

szx = size( thetaselLrg );
ismatweight = all(size(weights)==size(thetaselLrg,1)*[1 1]);

[dummy, L ] = laplaceMulND(zeros([szx(2:end) 1]),1, spacing(2:end));
% L = [Lo, Li]
% tht = [thtO'; thtI']
% f = .5*trace( (tht'*L')' * W *(L * tht)' )
%   =  trace( (thtO*Lo')' * W * (Lo * thtO')' )  % => C0
%    + trace( (thtO*Lo')' * W * (Li * thtI')' )  % \
%     +trace( (thtI*Li')' * W * (Lo * thtO')' )  % \+ => G0
%    + trace( (thtI*Li')' * W * (Li * thtI')' )  % => H

%  trace( (thtO*Lo')' * W * (Li * thtI')' )
%  = sum_ijklm   Lo_ij thtO_kj W_kl Li_im thtI_lm
%  = sum_lm    (sum_ijk Lo_ij thtO_kj W_kl Li_im ) thtI_lm
%  = sum_lm    (sum_i sum_k( sum_j( Lo_ij thtO_kj) W_kl ) Li_im ) thtI_lm

thetaselLrg(:,subindexsel) = 0;
C1o = thetaselLrg(:,:)*L';
if ismatweight
    C1w = weights * C1o;
    wm = weights;
else
    C1w = bsxfun( @times, C1o, weights(:));
    wm = diag(weights(:));
end;
C0 =.5*( C1w(:)'*C1o(:));
Li = L(:,subindexsel);
G0 = (C1w * Li);

H = kron( Li'*Li, wm);

laplaceFun = makelaplaceFun(C0, G0(:), H);
laplaceHessMulFun = [];
if 0 
    %% test if function is correct:
    xtst = randn(szx(1),numel(subindexsel));
    [f,g,h]= laplaceFun(xtst);
    xtstf = thetaselLrg;
    xtstf(:,subindexsel) = xtst;
    [ft,gt,ht] = regularizerfuncs{1}(xtstf);
    [f ,ft]
    [g, gt(:,subindexsel)]
    %% numerical derivative:
    step = .001;
    a = step*[eye(8) -eye(8)];
    fnum = zeros(size(a,1),2);
    gnum = cell(size(a,1),2);
    fnum2 = zeros(size(a,1),2);
    gnum2 = cell(size(a,1),2);
    fnum3 = zeros(size(a,1),2);
    for k=1:size(a,2);
        [fnum(k),gnum{k}]= laplaceFun(xtst+a(:,k));
        xtstf(:,subindexsel) = xtst+a(:,k);
        [fnum2(k),gnum2{k},ht] = regularizerfuncs{1}(xtstf);
        gnum2{k} = gnum2{k}(:,subindexsel);
        fnum3(k) = 100*xtstf(end,:)*L'*L*xtstf(end,:)'/2;
    end;
    [fnum*[1;-1]/(2*step) , ...
    fnum2*[1;-1]/(2*step) ,...
    fnum3*[1;-1]/(2*step)]
    Hnum = ([gnum{:,1}]-[gnum{:,2}])/(2*step)
    Hnum2 = ([gnum2{:,1}]-[gnum2{:,2}])/(2*step)
end;

function laplaceFun = makelaplaceFun(C0, G0, H)
laplaceFun = @(x)  laplaceEvalFun(x, C0, G0, H);

function [f,g,H] = laplaceEvalFun( x, C0, G0, H)
g = G0 + H*x(:);
f = C0 + G0'*x(:) + .5*(x(:)'*H*x(:));

% function laplaceHessMulFun