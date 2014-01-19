function pass = test_airy(pref)

% Grab some preferences
if ( nargin == 0 )
    pref = chebpref();
end

% Choose a domain:
a = -1;
b = 5;

% Make alinear chebfun:s
x = chebfun(@(x) x, [a, b], pref);

% Make a test space:
xx = linspace(a, b, 100);

% Test using x and (1+1i)*x;
for im = [0 1]
    
    % Initialise k:
    k = 1;

    % Test each of the K options of AIRY().
    for K = 0:3

        % [TODO] Implement scale = 1 (requires SINGFUNS)
        for scale = 0

            F = @(x) airy(K, (1+im*1i)*x, scale);
            f = F(x);
            pass(im+1,k) = norm(feval(f,xx) - F(xx), inf) < 20*max(f.epslevel.*f.vscale);
            k = k + 1;

        end

    end

end

end

