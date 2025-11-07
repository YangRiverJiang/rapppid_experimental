function x = ins_vec_el(vec, pos, el)
% insert an element into vector at defined position(s)
%
% INPUT: 
%  	vec:  	vector where an element should be inserted, [| vector]
%  	pos:   	indices where element should be inserted, places where
%               elements are standing after inserting
% 	el:    	element to insert
% OUPUT:
%  	x:      (new) vector with inserted element(s), [| vector]
% 
% Revision:
%   2025/09/12, MFWG: replace loop with a vectorized version
% 
% This function belongs to raPPPid, Copyright (c) 2023, M.F. Glaner
% *************************************************************************

% number of elements
n = numel(vec);
m = numel(pos);

% initialize output vector
x = zeros(n+m, 1);

% logical mask for insert positions
mask = false(n+m, 1);
mask(pos) = true;

% Fill result
x(mask)  = el;
x(~mask) = vec;

end

