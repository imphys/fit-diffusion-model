function [funout] = make2arg_anonfun( fun, varargin)
funout = @(p1,p2) fun(p1,p2, varargin{:});
