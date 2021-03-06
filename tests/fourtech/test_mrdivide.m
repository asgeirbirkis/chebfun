% Test file for fourtech/mrdivide.m

function pass = test_mrdivide(pref)

% Get preferences.
if ( nargin < 1 )
    pref = fourtech.techPref();
end

testclass = fourtech();

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Random number to use as a scalar constant.
alpha = -0.194758928283640 + 0.075474485412665i;


%%
% Check division of a FOURTECH object by a numeric array.

f_op = @(x) [exp(sin(pi*x)) 3./(4-cos(pi*x))];
f = testclass.make(f_op, [], pref);
pass(1) = isnan(f / 0);

g = f / alpha;
g_exact = @(x) f_op(x) ./ alpha;
pass(2) = norm(feval(g, x) - g_exact(x), inf) < ...
    50*max(g.vscale.*g.epslevel);

% A "least-squares" case where the solution is obvious.
I = eye(2);
g = f / I;
err = g*I - f;
pass(3) = max(max(abs(feval(err, x)))) < 10*max(g.vscale.*g.epslevel);

% A less trivial least-squares case for which we still know the answer.
A = [1 1];
g = f / A;
g_exact = @(x) (exp(sin(pi*x)) + 3./(4-cos(pi*x)))/2;
pass(4) = norm(feval(g, x) - g_exact(x), inf) < ...
    10*max(g.vscale.*g.epslevel);

%%
% Check division of a numeric array by a FOURTECH object.

f = testclass.make(@(x) cos(sin(pi*x)));
g = alpha / f;
pass(5) = abs(innerProduct(f, g) - alpha) < 10*max(g.vscale.*g.epslevel);

f = testclass.make(@(x) [sin(2*pi*x) cos(2*pi*x)]);
g = [1 1]/f;
g_exact = @(x) (sin(2*pi*x) + cos(2*pi*x));
pass(6) = norm(feval(g, x) - g_exact(x), inf) < ...
    10*max(g.vscale.*g.epslevel);
    
%%
% Check error conditions.

% Catch dimension mismatch errors.
try
    g = f / [1 2 3];
    pass(7) = false;
catch ME
    pass(7) = strcmp(ME.identifier, 'CHEBFUN:FOURTECH:mrdivide:size');
end

% Can't do f/g if both f and g are fourtech objects.
try
    f = testclass.make(@(x) sin(pi*x));
    g = testclass.make(@(x) cos(pi*x));
    h = f / g;
    pass(8) = false;
catch ME
    pass(8) = strcmp(ME.identifier, ...
        'CHEBFUN:FOURTECH:mrdivide:fourtechDivFourtech');
end

% Can't call mrdivide on a FOURTECH and a non-FOURTECH or non-double
% object.
try
    f = testclass.make(@(x) sin(x));
    g = f / true;
    pass(9) = false;
catch ME
    pass(9) = strcmp(ME.identifier, 'CHEBFUN:FOURTECH:mrdivide:badArg');
end

end
