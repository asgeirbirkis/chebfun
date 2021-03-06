% Test file for fourtech/plus.m

function pass = test_plus(pref)

% Get preferences.
if ( nargin < 1 )
    pref = fourtech.techPref();
end

testclass = fourtech();

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% A random number to use as an arbitrary additive constant.
alpha = -0.194758928283640 + 0.075474485412665i;

%%
% Check operation in the face of empty arguments.

f = testclass.make();
g = testclass.make(@(x) x, [], pref);
pass(1) = (isempty(f + f) && isempty(f + g) && isempty(g + f));

%%
% Check addition with scalars.

f_op = @(x) exp(sin(pi*x));
f = testclass.make(f_op, [], pref);
pass(2:3) = test_add_function_to_scalar(f, f_op, alpha, x);

%%
% Check addition of two FOURTECH objects.

f_op = @(x) zeros(size(x));
f = testclass.make(f_op, [], pref);
pass(4:5) = test_add_function_to_function(f, f_op, f, f_op, x);

f_op = @(x) exp(cos(pi*x)) - 1;
f = testclass.make(f_op, [], pref);

g_op = @(x) sin(100*pi*x);
g = testclass.make(g_op, [], pref);
pass(6:7) = test_add_function_to_function(f, f_op, g, g_op, x);

g_op = @(x) sin(cos(10*pi*x));
g = testclass.make(g_op, [], pref);
pass(8:9) = test_add_function_to_function(f, f_op, g, g_op, x);

%%
% Check operation for array-valued FOURTECH objects.

f_op = @(x) [zeros(size(x)) zeros(size(x)) zeros(size(x))];
f = testclass.make(f_op, [], pref);
pass(10:11) = test_add_function_to_function(f, f_op, f, f_op, x);

f_op = @(x) [sin(10*pi*x) sin(cos(pi*x)) exp(cos(pi*x))];
f = testclass.make(f_op, [], pref);
pass(12:13) = test_add_function_to_scalar(f, f_op, alpha, x);

g_op = @(x) [sin(pi*x) exp(1i*pi*x).*exp(1i*pi*x) cos(pi*x)];
g = testclass.make(g_op, [], pref);
pass(14:15) = test_add_function_to_function(f, f_op, g, g_op, x);

% This should fail with a dimension mismatch error.
g_op = @(x) sin(10*x);
g = testclass.make(g_op, [], pref);
try
    h = f + g; %#ok<NASGU>
    pass(16) = false;
catch ME
    pass(16) = strcmp(ME.message, 'Matrix dimensions must agree.');
end

%%
% Check that direct construction and PLUS give comparable results.

tol = 10*eps;
f = testclass.make(@(x) sin(pi*cos(3*pi*x)), [], pref);
g = testclass.make(@(x) cos(pi*sin(10*pi*x)) - 1, [], pref);
h1 = f + g;
h2 = testclass.make(@(x) sin(pi*cos(3*pi*x)) + (cos(pi*sin(10*pi*x)) - 1), [], pref);

pass(17) = norm(h1.coeffs - h2.coeffs, inf) < tol;

%%
% Check that adding a FOURTECH and an unhappy FOURTECH gives an
% unhappy result.

f = testclass.make(@(x) cos(pi*x));    % Happy
g = testclass.make(@(x) cos(x));   % Unhappy
h = f + g;  % Add unhappy to happy.
pass(18) = (~g.ishappy) && (~h.ishappy);
h = g + f;  % Add happy to unhappy.
pass(19) = (~g.ishappy) && (~h.ishappy);

%%
% Test addition of array-valued scalar to array-valued FOURTECH.

f = testclass.make(@(x) exp([sin(pi*x) cos(pi*x) -sin(pi*x).^2]));
g = f + [1 2 3];
g_exact = @(x) [exp(sin(pi*x))+1 exp(cos(pi*x))+2 exp(-sin(pi*x).^2)+3];
err = feval(g, x) - g_exact(x);
pass(20) = norm(err(:), inf) < 10*max(g.vscale.*g.epslevel);

%%
% Test scalar expansion in FOURTECH argument.

f = testclass.make(@(x) sin(pi*x));
g = f + [1 2 3];
g_exact = @(x) [(1 + sin(pi*x)) (2 + sin(pi*x)) (3 + sin(pi*x))];
err = feval(g, x) - g_exact(x);
pass(21) = isequal(size(g.coeffs, 2), 3) && norm(err(:), inf) < ...
    10*max(g.vscale.*g.epslevel);

end

% Test the addition of a FOURTECH F, specified by F_OP, to a scalar ALPHA using
% a grid of points X in [-1  1] for testing samples.
function result = test_add_function_to_scalar(f, f_op, alpha, x)
    g1 = f + alpha;
    g2 = alpha + f;
    result(1) = isequal(g1, g2);
    g_exact = @(x) f_op(x) + alpha;
    result(2) = norm(feval(g1, x) - g_exact(x), inf) <= ...
        50*max(g1.vscale.*g1.epslevel);
end

% Test the addition of two FOURTECH objects F and G, specified by F_OP and
% G_OP, using a grid of points X in [-1  1] for testing samples.
function result = test_add_function_to_function(f, f_op, g, g_op, x)
    h1 = f + g;
    h2 = g + f;
    result(1) = isequal(h1, h2);
    h_exact = @(x) f_op(x) + g_op(x);
    norm(feval(h1, x) - h_exact(x), inf);
    result(2) = norm(feval(h1, x) - h_exact(x), inf) <= ...
        50*max(h1.vscale.*h1.epslevel);
end
