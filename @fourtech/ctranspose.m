function ctranspose(f) %#ok<*INUSD>
%CTRANSPOSE   FOURTECH objects are not transposable, so this method will throw an
%error.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:FOURTECH:ctranspose:notpossible', ...
    'FOURTECH objects are not transposable.')
    
end
