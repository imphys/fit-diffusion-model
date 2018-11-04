function [ dist ] = levenshtein_distance( str1 , str2 , varargin )
% [ dist, idx ] = levenshtein_distance( str1 , str2 , options )
% Compute the minimum distance between str1 and str2
%
% INPUTS:
%   str1 : a single string or a cell array with strings of which the 
%          distance to str2 has to be computed.
%   str2 : a single string or a cell array with with strings.
%          The minimum distance of str1 elements to all str2 elements is
%          returned.
% options; as option-value pairs or as options structure
%   insertion : cost of inserting a character of str2 into str1
%   deletion  : cost of deleting a character of str1
%   replacement: cost of replacing a not equal character of str1 with a
%                character of str2
%   caseCost  : if specified the cost of replacing a lowercase value with
%                an uppercase value and vice-versa. Overrides 'replacement'
%                in this situation.
%   replace_values : an  n x 2 character matrix with pairs of characters
%                    for which you want to specify an alternative cost
%   replace_costs  : an  n element vector with the alternative cost of the
%                    replacements specified in replace_values
%
% OUTPUTS:
%   dist   : a numel(str1) x numel(str2) matrix with minimum distance of str1 to str2
%            for cell str1 or str2 inputs. size is 1 if a string is provided as input. 
%
% Examples:
% [ dist ] = levenshtein_distance( 'test' , 'testing' )
%  > dist =
%  >    3
% [ dist ] = levenshtein_distance( 'test' , {'testing' ,'TEST'},'caseCost',.3 )
%  > dist =
%  >   3.0   1.2
%
% Created by Dirk Poot, Erasmus MC, 21-11-2012
% Based on code from 
% http://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Octave_And_MATLAB

opts.insertion = 1;
opts.deletion = 1;
opts.replacement = 1;
opts.replace_values = [];
opts.replace_costs = [];
opts.caseCost = []; 
% TODO:
%   maxDist   : Stop comparing when cost reaches or exceeds this value.
%               Any returned distance larger or equal to this value is a 
%               lower bound to the actual distance.
%   - dont know how to efficiently reduce number of computations.

if nargin>2
    opts = parse_defaults_optionvaluepairs( opts, varargin{:});
end;
if ~isempty(opts.caseCost) && opts.caseCost~=opts.replacement
    % Case difference is special case.
    lcase = ('a':'z')';
    ucase = ('A':'Z')';
    opts.replace_values = [lcase ucase ; ucase lcase; opts.replace_values];
    opts.replace_costs = [ opts.caseCost * ones(numel(ucase)*2,1); opts.replace_costs];
end;

if iscell(str1) 
    n1 = numel(str1);
else
    n1 = 1;
end;
if iscell(str2) 
    n2 = numel(str2);
else
    n2 = 1;
end;
dist = inf(n1,n2);
for idx2 = 1 : n2
    % get string to currently compare with
    if iscell(str2)
        str2i = str2{idx2};
    else
        str2i = str2;
    end;
    for idx1 = 1 : n1
        % get string to currently compare with
        if iscell(str1)
            str1i = str1{idx1};
        else
            str1i = str1;
        end;

        % initialize costs
        L1 = length(str1i) + 1;
        L2 = length(str2i) + 1;
        cost = zeros(L1,L2);

        cost(:,1) = ( (0:L1-1) * opts.deletion)';
        cost(1,:) =   (0:L2-1) * opts.insertion;

        % compute character replacement costs:
        score = opts.replacement * bsxfun(@ne, str1i(:), reshape( str2i ,1,[]) );
        if ~isempty( opts.replace_values )
            [I,J] = find(score);
            [tf, loc] = ismember( [reshape(str1i(I),[],1) reshape(str2i(J),[],1)],opts.replace_values ,'rows');
            replcost = opts.replace_costs( loc(tf) );
            Irepl = sub2ind( size(score), I(tf), J(tf) );
            score(Irepl) = replcost;
        end;

        % Fill cost matrix with minimum cost:
        for i2 = 2:L2
            for i1 = 2:L1;
                m1 = cost(i1-1 , i2-1) + score(i1-1, i2-1); % replace character
                m2 = cost(i1-1 , i2)   + opts.deletion;     % delete character from str1
                m3 = cost(i1   , i2-1) + opts.insertion;    % insert character into str1
                cost(i1,i2) = min( [m1, m2, m3] );
            end
        end

        % store final minimum cost:
        dist(idx1, idx2) = cost(L1,L2);
    end;
end;