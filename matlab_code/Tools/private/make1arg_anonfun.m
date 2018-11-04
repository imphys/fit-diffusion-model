function [funout] = make1arg_anonfun( fun, varargin)
funout = @(p) fun(p, varargin{:});
