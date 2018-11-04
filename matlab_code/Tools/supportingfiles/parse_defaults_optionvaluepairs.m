function out = parse_defaults_optionvaluepairs( defaults , varargin )
% out = parse_defaults_optionvaluepairs( defaults , option1, value1, ...)
% out = parse_defaults_optionvaluepairs( defaults , option_value_struct)
% Helper function to parse option value pairs and update the fields in the defaults structure.
%
% INPUTS:
% defaults : structure with default values for all possible options.
% option_i : string with name of options, case sensitive matched to the fields in 
%            defaults. Errors if the option is not found in the defaults structure.
%            For scalar structures in defaults, the subfields can be 
%            assigned by specifying 'option.sub_field_name' as option_i.
% value_i  : anything (that can be assigned to a field of a structure).
% option_value_struct : structure with valid fields (each field should also be in defaults)
%            The fields override the default values.
% 
% OUTPUTS:
% out      : structure defaults, with all fields specified by the options 
%            replaced by the corresponding values.
%            No checking of the values is performed.
%
% Created by Dirk Poot (d.poot@erasmusmc.nl), Erasmus medical center, Rotterdam

% %          There is one special field 'addUnknownFields'.
% %   Reason not to include this option:
% %   (Performance and) We do not want to allow options that we do not know
% %   how to use/interpret as this will cause too many silent bugs (ie.
% %   mistyping an option will cause it to be ignored.)
% addUnknownFields = isfield(defaults,'addUnknownFields') && defaults.addUnknownFields;
% if addUnknownFields
%    defaults = rmfield(defaults,'addUnknownFields');
% end;

out = defaults;
if mod(nargin,2)~=1
    % if not odd number of input arguments => no well defined option-value pairs; maybe option_value_struct 
    if nargin==2 && isstruct(varargin{1}) && numel(varargin{1})==numel(defaults)
        % option, values are already in structure.
        fieldnm = fieldnames(varargin{1});
        for k=1:numel(fieldnm)
            if ~isfield(out,fieldnm{k})
                allfieldnamesout = fieldnames(out);
                k_invalid = k;
                for k2=k+1:numel(fieldnm)
                    if ~isfield(out,fieldnm{k2})
                        k_invalid(end+1) = k2; %#ok<AGROW> ; only for displaying an error. That is not time critical.
                    end;
                end;
                errorstring = disperror( fieldnm(k_invalid), allfieldnamesout, true);
                error( errorstring );
            end;
            if isstruct(out.(fieldnm{k})) && isstruct(varargin{1}.(fieldnm{k}))
                % if field of default is structure and field in options is structure, merge structures by calling self again:
                out.(fieldnm{k}) = parse_defaults_optionvaluepairs( out.(fieldnm{k}) , varargin{1}.(fieldnm{k}) );
            else
                out.(fieldnm{k}) = varargin{1}.(fieldnm{k});
            end;
        end;
        return;
    else
        error('arguments should come in option+value pairs.');
    end;
end;
npairs = (nargin-1)/2;
isopt = [true(1,npairs);false(1,npairs)];
for k=1:npairs
    if ~ischar(varargin{k*2-1}) 
        isopt(1,k) = false;
        break;
    end;
end;
[allfieldnm , nnormal] = expandedfieldnames(defaults);

[isopt(isopt), idx ] = ismember( varargin(isopt), allfieldnm);
haserror = ~all(isopt(1,:));
if haserror
    k2= find( ~isopt(1,:) );
    invalidoptions = varargin( k2*2-1 );
    errorstring = disperror( invalidoptions, allfieldnm, false);
    error( errorstring );
    return; % actually 
end;

for k=1:numel(idx)
    if idx(k) <= nnormal
        if isstruct( out.( allfieldnm{ idx(k) } ) )  && isstruct( varargin{k*2} )
            out.( allfieldnm{ idx(k) } ) = parse_defaults_optionvaluepairs( out.( allfieldnm{ idx(k) } ) , varargin{k*2} );
        else
            out.( allfieldnm{ idx(k) } ) = varargin{k*2};
        end;
    else
        eval( ['out.' allfieldnm{ idx(k) } '=varargin{k*2};']);
    end;
end;



function [list, norig] = expandedfieldnames( s )
list = fieldnames( s );
norig = numel(list);
adlist = {};
for k=1:norig
    if isstruct(s.(list{k})) && numel(s.(list{k}))==1
        adlist{end+1} = expandedfieldnames( s.(list{k}) );
        for k2 = 1:numel(adlist{end})
            adlist{end}{k2} = [list{k} '.' adlist{end}{k2}];
        end;
    end;
end;
if ~isempty(adlist)
    list = vertcat(list,adlist{:});
end;

function errorstring = disperror( invalidoptions, alloptions, isstruct)
allfieldsstr = sprintf('"%s", ', alloptions{:});
alloptsstr = cell(1,numel(invalidoptions));
reploptstr = cell(1,numel(invalidoptions));
for k2= 1:numel(invalidoptions)
    if ischar(invalidoptions{k2})
        alloptsstr{k2} = ['"' invalidoptions{k2} '", '];
        if exist('levenshtein_distance.m','file')==2
            [dist] =  levenshtein_distance( invalidoptions{k2},alloptions,'caseCost',.2);
            replopt = ( dist <= 1.1*min(dist) );
            reploptstr{k2} = sprintf('"%s", ',alloptions{replopt});
        end;
    else
        alloptsstr{k2} = sprintf('class of option %d: %s, ',k2, class(invalidoptions{k2}));
    end;
end;
alloptsstr= [alloptsstr{:}];
reploptstr = [reploptstr{:}];
%     error('Each option should be a string identifying a valid option.\nAll valid options: %s\nProvided options : %s', allfieldsstr(1:end-2), alloptsstr(1:end-2));
if isstruct
    startstr = 'Each field of your option structure should identify a valid option.' ;
else
    startstr = 'Each option should be a string identifying a valid option.' ;
end;
if isempty(reploptstr)
    errorstring = sprintf( '%s\n\nProvided invalid options : %s\nAll valid options        : %s', startstr, alloptsstr(1:end-2), allfieldsstr(1:end-2) );
else
    errorstring = sprintf( '%s\n\nProvided invalid options : %s\nLikely valid options     : %s\n\nAll valid options        : %s', startstr, alloptsstr(1:end-2),reploptstr(1:end-2), allfieldsstr(1:end-2) );
end;