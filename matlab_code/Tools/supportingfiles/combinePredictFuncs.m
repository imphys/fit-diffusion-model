function [A, dA] = combinePredictFuncs( funcs , nparfuncs,  theta , mixing, offset, fields)
%     [A, dA] = combinePredictFuncs( funcs , nparfuncs, theta , mixing , offset , fields)
% or: [A, dA] = combinePredictFuncs( funcs , nparfuncs, theta , thetamodiffun ,[], fields)
%
% function that can be used as measurement prediction function in 
% fit_MRI, CramerRaoLowerBound_MRI, etc. Combines several prediction
% functions. Convenient to test combinations of measurements that have some
% shared variables. If you want to repeat the same function multiple times
% 'repeatMeasurementFun' might be slightly more effecient.
%
% INPUTS:
%  funcs   : cell array with functions that predicts outputs of a single experiment.
%            (that is: each function can be used with fit_MRI)
%  nparfuncs: column vector of which element (k,1) specifies size( theta_fun{k} , 1) 
%             which is the theta input of function k. 
%             If fields is provided a second column is required that
%             specifies the number of fields passed to funcs. 
%             use nparfuncs(k,2)==-1 to specify that no field
%             should be passed to funcs{k}. 
%  theta   : (reduced) parameter vector, (which is optimized by fit_MRI)
%  mixing  : mixing matrix
%  offset  : vector with offset for thetafun (see below)
%  thetamodiffun : theta modification function that allows non linear
%                   warping. Same interface as an element of funcs. See
%                   below for useage.
%  fields  : provides fields arguments to funcs. The second column of nparfuncs
%            specifies the number of fields passed to funcs{k}. If the
%            fields argument is not present or nparfuncs(k,2)==-1 no field
%            argument is passed to funcs{k}. The fields should be provided
%            in order and can only be direcly passed through. (It is never
%            required to modify fields in this function since the fields
%            argument can be constructed arbitrarily before calling this
%            function)
%
%  thetafun_full = mixing * theta + offset;
%  or: [thetafun_full , grad_thetafun_full] = thetamodiffun( theta );
%  theta_fun{k} = mat2cell( thetafun_full, nparfuncs, size(thetafun_full,2) );
%  A_fun = [ funcs{  1}( theta_fun{ 1} ) );
%             ....
%            funcs{end}( theta_fun{end} ) )];
%
% E.g. using joint estimation of FA_scale: 
%
% Created by Dirk Poot, TUDelft, 4-2-2014
% Original version based on repeatMeasurementFun

hasfields = nargin>=6;

% compute thetafun_full:
if isa( mixing, 'function_handle')
    if nargout>1
        [thetafun_full , grad_thetafun_lrg ] = mixing( theta(:,:) );
    else
        thetafun_full = mixing( theta(:,:) );
    end;
else
    if nargin<5 || isempty(offset)
        thetafun_full = mixing * theta(:,:) ;
    else
        thetafun_full = bsxfun( @plus, mixing * theta(:,:) , offset);
    end;
end;
% split for each function:
theta_fun =  mat2cell( thetafun_full, nparfuncs(:,1), size(thetafun_full,2) ) ;

% compute functions (and derivative)
A_fun = cell(1,numel(funcs));
if nargout>1
    dA_fun =cell(size(A_fun));
    if hasfields
        fieldidx_start = 1;
        for k=1:numel(funcs)
            if nparfuncs(k,2)==-1
                [A_fun{k} , dA_fun{k} ]= funcs{k}( theta_fun{ k } );
            else
                fieldidx_end = fieldidx_start + nparfuncs(k,2);
                fields_sel = fields( fieldidx_start :fieldidx_end-1,:);
                [A_fun{k} , dA_fun{k} ]= funcs{k}( theta_fun{ k } , fields_sel);
                fieldidx_start = fieldidx_end;
            end;
        end;
    else
        for k=1:numel(funcs)
            [A_fun{k} , dA_fun{k} ] = funcs{k}( theta_fun{k} );
        end;
    end;
    A = vertcat(A_fun{:});
    szdAfun = [size(A,1), size(thetafun_full,2) , sum(nparfuncs(:,1)) ];
    dA_funl = zeros( szdAfun , class(dA_fun{1}) );
    startA = 1;
    startP = 1;
    for k=1:numel(funcs)
        endA = startA+size(A_fun{k},1);
        endP = startP+nparfuncs(k);
        dA_funl(startA:endA-1,:,startP:endP-1) = dA_fun{k};
        startA = endA;
        startP = endP;
    end;
    if isa( mixing, 'function_handle')
        error('not correctly implemented yet'); % next lines are not yet correct but just copied from repeatMeasurementFun
        dA = zeros( [szdAfun(1)*szdAfun(2), size(theta,1) ] , class(dA_fun) );
        for k = 1 : size(thetafun_full,2)
            dA( :, k , : ) = reshape( dA_funl( :,: , k ,: ,: ), szdAfun(1)*szdAfun(2), size(thetafun_full,1) ) * reshape( grad_thetafun_lrg(:, k , :) , size(thetafun_full,1), size(theta,1));
        end;
    else
        dA = reshape( reshape( dA_funl , szdAfun(1)*szdAfun(2), size(mixing,1) ) * mixing, [szdAfun(1), szdAfun(2), size(mixing,2)]);
    end;    
else
    if hasfields
        fieldidx_start = 1;
        for k=1:numel(funcs)
            if nparfuncs(k,2)==-1
                A_fun{k} = funcs{k}( theta_fun{ k } );
            else
                fieldidx_end = fieldidx_start + nparfuncs(k,2);
                fields_sel = fields( fieldidx_start :fieldidx_end-1,:);
                A_fun{k} = funcs{k}( theta_fun{ k } , fields_sel);
                fieldidx_start = fieldidx_end;
            end;
        end;
    else
        for k=1:numel(funcs)
            A_fun{k} = funcs{k}( theta_fun{ k } );
        end;
    end;
    A = vertcat(A_fun{:});
end;

function test_fun;
%%
disp(' ');
disp(repmat('*',[1 100]));
disp(' ');
disp('Testing combinePredictFuncs');
disp(' ');

Tr = [10;10;10;20]; % all repetition times in ms
nominal_alphas=[3;20; 40;80 ]; % all nominal flip angles.
sigma = 1;
thttst = [500;1100;1];
numreps = 4;

predfun = @(par) pred_S_T1_alpha( par , Tr, nominal_alphas );
[val, grad ]= predfun(thttst);
[gradnum]= jacobianest( predfun, thttst );
[C,I,J] = CramerRaoLowerBound_MRI( thttst, predfun , sigma);
disp('Difference between numerical and analytical gradient (original function):')
disp( squeeze(grad )-gradnum)

predfuncs = repmat({predfun},[1 numreps]);
mixing = eye( numreps * size(thttst,1) ); 
nonrepparidx = size(thttst,1)*(1:numreps-1); % dont' repeat last parameter
mixing(:, nonrepparidx )=[]; mixing( nonrepparidx, end )=1;
predfun_lrg = @(par) combinePredictFuncs( predfuncs , 3*ones(1,numreps), par , mixing);
thttst_lrg = thttst([repmat( 1:size(thttst,1)-1 ,1,numreps) size(thttst,1)]);
[vall, gradl ]= predfun_lrg(thttst_lrg);
[gradnuml]= jacobianest( predfun_lrg, thttst_lrg );
disp('Difference between numerical and analytical gradient (expanded function/mixing):')
disp( squeeze(gradl )-gradnuml)
[Cl,Il,Jl] = CramerRaoLowerBound_MRI( thttst_lrg, predfun_lrg , sigma);

disp('sqrt diag CRLB of original problem:');
disp(sqrt( C(I==J,:) ))

disp('sqrt diag CRLB of expanded problem:');
disp(sqrt( Cl(Il==Jl,:) ))
%%

