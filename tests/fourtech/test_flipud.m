% Test file for fourtech/flipud.m

function pass = test_flipud(pref)

% Get preferences.
if ( nargin < 1 )
    pref = fourtech.techPref();
end
    
testclass = fourtech();

%%
% Conduct a few very straightforward tests.
f = testclass.make(@(x) sin(pi*x), [], pref);
g = testclass.make(@(x) -sin(pi*x), [], pref);
h = flipud(f);
pass(1) = norm(g.coeffs - h.coeffs, inf) < 10*h.vscale.*h.epslevel;

f = testclass.make(@(x) [sin(sin(pi*x)), exp(1i*pi*x)], [], pref);
g = testclass.make(@(x) [-sin(sin(pi*x)), exp(-1i*pi*x)], [], pref);
h = flipud(f);
pass(2) = norm(g.coeffs - h.coeffs, inf) < 100*max(h.vscale.*h.epslevel);

end
