function [regulariser, border] = parse_spatialRegularizer( opts )
% [regulariser] = parse_spatialRegularizer( opts )
% Helper function that parses spatialRegulariser[LS] arguments in opts
% and creates a single spatial regulariser function.
%
% Created by Dirk Poot, Erasmus MC
% created: 17-1-2013, moved from fit_MRI and updated to allow multiple
%                     regularizing terms (which are summed up)
    
if ~iscell(opts.spatialRegularizer)
    opts.spatialRegularizer = {opts.spatialRegularizer , [] };
    if ~iscell(opts.spatialRegularizerWeights)
        opts.spatialRegularizerWeights = {opts.spatialRegularizerWeights};
    end;
end;
if size(opts.spatialRegularizer,2)<2
    opts.spatialRegularizer(1,2)={[]};
end;
if isa(opts.spatialRegularizerLS,'function_handle')
    opts.spatialRegularizerLS = {opts.spatialRegularizerLS ,[]};
end;
if ~isempty(opts.spatialRegularizerLS)
    for k=1:size(opts.spatialRegularizerLS,1)
        spatialRegularizerLS = opts.spatialRegularizerLS{1};
        spatialRegularizerLS_JacMul = opts.spatialRegularizerLS{2};
        spatialRegularizerFun = make1arg_anonfun(@LSfun2normalfun, spatialRegularizerLS, spatialRegularizerLS_JacMul);
        spatialRegularizerHessMul = LSfun2normalfun('hessmulfun',spatialRegularizerLS, spatialRegularizerLS_JacMul);
        opts.spatialRegularizer(end+1,:) = { spatialRegularizerFun , spatialRegularizerHessMul};
    end;
end;

k=1;
while k<=size(opts.spatialRegularizer,1)
    if ischar( opts.spatialRegularizer{k,1} ) 
%         if numel(opts.spatialRegularizerWeights)<k
        [opts.spatialRegularizer{k,1:3}, opts.spatialRegularizeBorder] = ...
            build_regularisation_function( opts.spatialRegularizer{k,1}, opts.spatialRegularizerWeights{k} , opts.voxelSpacing);
    elseif isempty(opts.spatialRegularizer{k,1})
        opts.spatialRegularizer(k,:) = [];  %remove row;
        k=k-1;
    elseif ~isa(opts.spatialRegularizer{1},'function_handle')
        error('unsuported spatialRegularizer value');
    end;
    k=k+1;
end;
if size(opts.spatialRegularizer,1)>1
    regulariser = makeaddfun( opts.spatialRegularizer );
else
    regulariser = opts.spatialRegularizer;
end;
border = opts.spatialRegularizeBorder;
   
function [spatialRegularizerfun, spatialRegularizer_hessmul, spatialRegularizer_prepare , border] = build_regularisation_function( spatialRegularizer, spatialRegularizerWeights , voxelSpacing)
spatialRegularizer_prepare = [];
switch spatialRegularizer
    case {'TV','total variation'}
        sp = [inf;voxelSpacing(:)];
        spatialRegularizerfun= make1arg_anonfun(@totalVariationVecRegularizer, sp , spatialRegularizerWeights,[],[],1) ;
        spatialRegularizer_hessmul = totalVariationVecRegularizer([] , [], [], [], [], 2) ;
        border = 1;
    case {'TV_cyclic'}
        sp = [inf;voxelSpacing(:)];
        spatialRegularizerfun= make1arg_anonfun(@totalVariationVecRegularizer, sp , spatialRegularizerWeights,[],[],1,1,true) ;
        spatialRegularizer_hessmul = totalVariationVecRegularizer([] , [], [], [], [], 2) ;
        border = 1;        
    case {'laplacian'}
        sp = [inf;voxelSpacing(:)];
        % TODO: make laplacian regularizer work with preconditioner and laplaceRegularizerHessMul
        spatialRegularizerfun =  make1arg_anonfun(@laplaceRegularizer, sp , spatialRegularizerWeights, 1 ) ;
        spatialRegularizer_hessmul =  make2arg_anonfun(@laplaceRegularizerHessMul, sp , spatialRegularizerWeights) ;
        spatialRegularizer_prepare =  @(regularizerfuncs, explicitHessian, thetaselLrg, subindexsel ) laplaceRegularizer_prepare( regularizerfuncs, explicitHessian, thetaselLrg, subindexsel, sp , spatialRegularizerWeights) ;
        border = 2;
    case {'laplacian_cyclic'} %TODO: handle differently; This only works correctly for global optimization. 
        sp = [inf;voxelSpacing(:)];
        % TODO: make laplacian regularizer work with preconditioner and laplaceRegularizerHessMul
        spatialRegularizerfun =  make1arg_anonfun(@laplaceRegularizer, sp , spatialRegularizerWeights, 1 , true) ;
        spatialRegularizer_hessmul =  make2arg_anonfun(@laplaceRegularizerHessMul, sp , spatialRegularizerWeights, true) ;
        spatialRegularizer_prepare =  @(regularizerfuncs, explicitHessian, thetaselLrg, subindexsel ) error('purposedly invalid; should not be called');
        border = 2;    
    case {'laplacianlog'}
        sp = [inf;voxelSpacing(:)];
        selRI = repmat({':'},[1 numel(sp)]);
        selRI = {selRI,selRI};
        selRI{1}{1} = 1:size(spatialRegularizerWeights,1);
        spatialRegularizerfun =  make1arg_anonfun(@laplaceLogRegularizer, sp , spatialRegularizerWeights, selRI , false) ;
        spatialRegularizer_hessmul =  spatialRegularizerfun('hessmulfun');
        spatialRegularizer_prepare =  @(regularizerfuncs, explicitHessian, thetaselLrg, subindexsel ) error('purposedly invalid; should not be called');
        border = 2;    
    case {'laplacian_Block'}
        sp = [inf;voxelSpacing(:)];
        % TODO: make laplacian regularizer work with preconditioner and laplaceRegularizerHessMul
        spatialRegularizerfun =  @(theta, opti, Binfo) Laplace_Regul_blocks(theta, sp , spatialRegularizerWeights, 1, opti ,Binfo ) ;
        spatialRegularizer_hessmul =  make2arg_anonfun(@laplaceRegularizerHessMul, sp , spatialRegularizerWeights) ;
        spatialRegularizer_prepare =  @(regularizerfuncs, explicitHessian, thetaselLrg, subindexsel ) laplaceRegularizer_prepare( regularizerfuncs, explicitHessian, thetaselLrg, subindexsel, sp , spatialRegularizerWeights) ;
        border = 2;
    otherwise
        error('Unrecognised spatialRegularizer method');
end;

function [regul] = makeaddfun( spatialRegularizer )
funs = spatialRegularizer(:,1);
funs_hessmul = spatialRegularizer(:,2);
regul = {@(val) addfun( val, funs), ...
         @(hinfo, x) addfun_hessmul( hinfo, x, funs_hessmul)};
    
function [f,g,h] = addfun( val, funs )
f=zeros(numel(funs),1);
g=0;
h = cell(numel(funs),1);
for k = 1 : numel(funs)
    if nargout<=1
        [f(k)] = funs{k}(val);
    else
        if nargout>2
            [f(k), gk, h{k}] = funs{k}(val);
        else
            [f(k), gk] = funs{k}(val);
        end;
        g = g + gk;
    end;
end;
f = sum(f);

function [hx]= addfun_hessmul( hinfo, x, funs_hessmul)
hx = 0;
for k=1:numel(funs_hessmul)
    if isempty(funs_hessmul{k})
        hx = hx + hinfo{k}*x;
    else
        hx = hx + funs_hessmul{k}( hinfo{k}, x);
    end;
end;