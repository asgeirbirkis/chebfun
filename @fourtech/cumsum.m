function f = cumsum(f, m, dim)
%CUMSUM   Indefinite integral of a FOURTECH.
%   CUMSUM(F) is the indefinite integral of the FOURTECH F, whose mean
%   is zero, with the constant of integration chosen so that F(-1) = 0.
%   If the mean of F is not zero then an error is thrown since the indefinite
%   integral would no longer be periodic.
%
%   CUMSUM(F, M) will compute the Mth definite integral with the constant of
%   integration chosen so that each intermediary integral evaluates to 0 at
%   -1.
%   Thus, CUMSUM(F, 2) is equivalent to CUMSUM(CUMSUM(F)).
%
%   CUMSUM(F, M, 2) will take the Mth cumulative sum over the columns F an
%   array-valued FOURTECH.
%
% See also DIFF, SUM.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the FOURTECH G of length n is represented as
%       \sum_{k=-(n-1)/2}^{(n-1)/2} c_k exp(i*pi*kx)
% its integral is represented with a FOURTECH of length n given by
%       \sum_{k=-(n-1)/2}^{(n-1)/2} b_k exp(i*pi*kx)
% where b_0 is determined from the constant of integration as
%       b_0 = \sum_{k=-(n-1)/2}^{(n-1)/2} (-1)^k/(i*pi*k) c_k;
% with c_0 := 0. The other coefficients are given by
%       b_k = c_k/(i*pi*k). 
%
% If the FOURTECH G of length n is represented as
%       \sum_{k=-n/2+1}^{n/2-1} c_k exp(i*pi*kx) + c(n/2)cos(n*pi/2x)
% then first set c(n) = 0.5*c(n) and define a = [0.5*c(n/2) c] so that we
% have the equivalent expansion:
%       \sum_{k=-n/2}^{n/2} a_k exp(i*pi*kx)
% The integral of this is represented with a FOURTECH of length n+1 given by
%       \sum_{k=-n/2}^{n/2} b_k exp(i*pi*kx)
% where b_0 is determined from the constant of integration as
%       b_0 = \sum_{k=-n/2}^{n/2} (-1)^k/(i*pi*k) a_k;
% with a_0 := 0. The other coefficients are given by
%       b_k = a_k/(ik).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trivial case of an empty FOURTECH:
if ( isempty(f) )
    return
end

if ( nargin < 2 || isempty(m) )
    % Order of intergration not passed in. Assume 1 by default:
    m = 1; 
elseif ( m == 0 )
    % Nothing to do here!
    return
end    

% Sum with respect to the continuous variable by default:
if ( nargin < 3 )
    dim = 1;
end

if ( dim == 1 )
    % Take difference across 1st dimension:
    f = cumsumContinuousDim(f, m);
else
    % Take difference across 2nd dimension:
    f = cumsumFiniteDim(f, m);
end

end

function f = cumsumContinuousDim(f, m)
% CUMSUM over the continuous dimension.

    % Initialize storage:
    c = f.coeffs; % Obtain Fourier coefficients {c_k}
    numCoeffs = size(c,1);
    
    fIsEven = mod(numCoeffs,2) == 0;

    % Check that the mean of the fourtech is zero.  If it is not, then
    % throw an error.
    if any(abs(c((numCoeffs+1-fIsEven)/2,:)) > f.vscale.*f.epslevel)
        error('CHEBFUN:FOURTECH:cumsum:meanNotZero', 'Indefinite integrals are only possible for FOURTECH objects with zero mean.');
    end
    
    % Force the mean to be exactly zero.
    if fIsEven
        c(numCoeffs/2,:) = 0;
        % Expand the coefficients to be symmetric (see above discussion).
        c(numCoeffs,:) = 0.5*c(numCoeffs,:);
        c = [c(numCoeffs,:);c];
        highestDegree = numCoeffs/2;
    else
        c((numCoeffs+1)/2,:) = 0;
        highestDegree = (numCoeffs-1)/2;
    end
    
    % Loop for integration factor for each coefficient:
    sumIndicies = (highestDegree:-1:-highestDegree).';
    integrationFactor = (-1i./sumIndicies/pi).^m;
    % Zero out the one corresponding to the constant term.
    integrationFactor(highestDegree+1) = 0;
    c = bsxfun(@times,c,integrationFactor);
    % If this is an odd order cumsum and there are an even number of
    % coefficients then zero out the cofficient corresponding to sin(N/2x)
    % term, since this will be zero on the Fourier grid.
    if mod(m,2) == 1 && fIsEven
        c(1,:) = 0;
        c(numCoeffs+1,:) = 0;
    end
    
    % Fix the constant term.    
    c(highestDegree+1,:) = -sum(bsxfun(@times,c,(-1).^sumIndicies));
    
    % If the original fourtech had an even number of coefficients then
    % shrink the coefficent vector corresponding to its indefinite integral
    % back to its original size since it was increased by one above to make
    % the integration code slicker.
    if fIsEven 
        c = c(2:end,:);
    end
            
    % Recover values and attach to output:
    f.values = f.coeffs2vals(c);
    f.values(:,f.isReal) = real(f.values(:,f.isReal));
    
    f.coeffs = c;

    % Update vscale
    f.vscale = max(abs(f.values), [], 1);

    % Update epslevel:
    f.epslevel = updateEpslevel(f);
    
    % Simplify (as suggested in Chebfun ticket #128)
    f = simplify(f);
    
    % Ensure f(-1) = 0:
    lval = get(f, 'lval');
    f.coeffs(end,:) = f.coeffs(end,:) - lval;
    f.values = bsxfun(@minus, f.values, lval);
    
end

function f = cumsumFiniteDim(f, m)
% CUMSUM over the finite dimension.

    for k = 1:m
        f.values = cumsum(f.values, 2);
        f.coeffs = cumsum(f.coeffs, 2);
        newVscale = max(abs(f.values), [], 1);
        epslevelApprox = sum(f.epslevel.*f.vscale, 2)/sum(newVscale, 2); % TODO: Is this right?        
        f.epslevel = updateEpslevel(f, epslevelApprox);
        f.vscale = newVscale;
    end

end
