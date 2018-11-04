function [varargout] = shuffleoutputs(fun, nout, select, funargs, sel, perm, resh )
%[varargout] = shuffleoutputargs(fun, nout, select [, args, sel, perm, resh] )
% Alows the caller to shuffle, select and reshape the output arguments of a function,
% called with arguments args. 
%
% INPUTS:
% fun    : function to be called
% nout   : number of output arguments requested from fun.
% select : a selection vector (logical or integer)
% args   : optional cell array with optional extra arguments passed to fun.
% sel    : optional cell array to select from outputs matrices. (indexing)
%          sel{k} is a cell array appied to varargout{k} (so after selection by select)
% perm   : optional cell array with permutation vectors, (permute function)
%          perm{k} is applied to varargout{k} (so after selection by select)
% resh   : optional cell array with reshape arguments, (reshape function)
%          resh{k} is applied to varargout{k} (so after selection and permutation)
%
% Any string arguments in sel and resh are evaluated before beiing applied. The size
% of the inputs and outputs is stored in variables, so can be used by the evaluation.
% Sizes of inputs  : ['I' #input 'd' #dimension] 
% Sizes of outputs : ['O' #output 'd' #dimension] 
%    (outputs of function, so before selection by select)
% E.g. size of input 1 in dimension 2 is called 'I1d2'.
% 
% Example:
%  fun    = @MyFun;
%  nout   = 2;
%  select = 2;
%  args   = {In1, In2, In3};
%  sel    = {{':','1:I1d2',1:4}}
%  perm   = {[2 1 3]};
%  resh   = {{'I1d2','O2d1*2',4}}
% Now 
%  outp = shuffleoutputs(fun, nout, select, funargs, sel, perm, resh )
% does the same as:
%  [Out1, Out2] = MyFun( In1, In2, In3);
%  tmp1 = Out2( : , 1:size(In1,2) , 1:4 );
%  tmp2 = permute(tmp1, [2 1 3]);
%  outp = reshape(tmp2, size(In1,2) , size(Out2,1)*2 , 4);
%
% OUTPUTS: 
% the shuffled outputs of function.
%
% Created by Dirk Poot, University of Antwerp
% 11-3-2009

if nargin<3
    select = 1:nout;
end;
if nargin<4
    funargs = {};
end;
noutp = cell(nout,1);
[noutp{:}] = fun(funargs{:});
varargout = noutp(select);
if nargin>=5 && ~isempty(sel)
    sel = parseargs(sel, funargs, noutp);
    for k=1:min(numel(sel),numel(varargout))
        if ~isempty(sel{k})
            varargout{k} = varargout{k}(sel{k}{:});
        end;
    end;
end;
if nargin>=6 && ~isempty(perm)
    for k=1:min(numel(perm),numel(varargout))
        if ~isempty(perm{k})
            varargout{k} = permute(varargout{k},perm{k});
        end;
    end;
end;
if nargin>=7 && ~isempty(resh)
    resh = parseargs(resh, funargs, noutp);
    for k=1:min(numel(resh),numel(varargout))
        if ~isempty(resh{k})
            varargout{k} = reshape(varargout{k},resh{k});
        end;
    end;
end;

function inp = parseargs(inp, funargs, noutp)
anyischar = false;
for k=1:numel(inp)
    for k2=1:numel(inp{k})
        if ischar(inp{k}{k2}) && ~isequal(inp{k}{k2},':')
            anyischar = true;
            break;
        end;
    end;
end;
if anyischar
    % Assign sizes of inputs to variables:
    for k = 1:numel(funargs)
        s = size(funargs{k});
        for k2 = 1:numel(s)
            eval(['I' num2str(k) 'd' num2str(k2) '=s(k2);']);
        end;
    end;
    % Assign sizes of outputs to variables:
    for k = 1:numel(noutp)
        s = size(noutp{k});
        for k2 = 1:numel(s)
            eval(['O' num2str(k) 'd' num2str(k2) '=s(k2);']);
        end;
    end;
    % evaluate string arguments:
    for k=1:numel(inp)
        for k2=1:numel(inp{k})
            if ischar(inp{k}{k2}) && ~isequal(inp{k}{k2},':')
                inp{k}{k2} = eval(inp{k}{k2});
            end;
        end;
    end;
end;