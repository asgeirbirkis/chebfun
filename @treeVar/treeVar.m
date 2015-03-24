classdef  (InferiorClasses = {?chebfun}) treeVar
    %TREEVAR   A class for analysing syntax trees of ODEs in Chebfun.
    %   The TREEVAR class allows Chebfun to analyse the syntax trees of ODEs in
    %   CHEBFUN. Its current use is to enable Chebfun to automatically convert
    %   (systems of) higher order ODEs to coupled first order systems. This is
    %   particularly useful for initial-value problems (IVPs), as that allows
    %   Chebfun to call one of the built-in MATLAB solvers for solving IVPs via
    %   time-stepping, rather than globally via spectral methods and Newton's
    %   method in function space.
    %
    %   This class is not intended to be called directly by the end user.
    %
    %   T = TREEVAR(ID, DOMAIN), where ID is a Boolean vector corresponding to 
    %   the order of variables in the problem, and DOMAIN is interval that the 
    %   problem is specified on, returns the TREEVAR object T, which stores the 
    %   ID and the DOMAIN. See example below for how the ID vector is specified.
    %
    %   T = TREEVAR() is the same as above, but with the default ID = 1, and
    %   DOMAIN = [-1, 1]. This is useful for quick testing purposes.
    %
    %   Example 1: Construct TREEVAR object for the scalar IVP
    %       u'' + sin(u) = 0
    %   on the interval [0, 10]:
    %       u = treeVar(1, [0 10]);
    %
    %   Example 2: Construct TREEVAR objects for the coupled IVP
    %       u'' + v = 1, u + v' = x
    %   on the interval [0, 5]:
    %       u = treeVar([1 0], [0 5]);
    %       v = treeVar([0 1], [0 5]);
    %
    % See also CHEBOP, CHEBOP/SOLVEIVP.
    
    % Copyright 2014 by The University of Oxford and The Chebfun Developers.
    % See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TREEVAR class description:
    %
    % The TREEVAR class is used by the CHEBOP class to convert higher order
    % ODEs to coupled systems of first order ODEs, which can then be solved
    % using one of the built-in MATLAB solvers, such as ODE113. This is done
    % by evaluating the (anonymous) functions in the .OP field of the CHEBOP
    % with TREEVAR arguments, which will construct a syntax tree of the
    % mathematical expression describing the operator. By then analysing the
    % syntax tree and restructuring it appropriately, conversion to a first-
    % order system is made possible.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        % The syntax tree of a TREEVAR variable, starting from an initial
        % variable. Each syntax tree is a MATLAB struct, which contains the
        % following fields:
        %    METHOD:    The method leading to the construction of the variable.
        %    NUMARGS:   Number of arguments to the method that constructed the
        %        variable.
        %    DIFFORDER: The differential order of the TREEVAR, which
        %        represents how many times the base variable(s) have been
        %        differentiated when before we arrive at the current TREEVAR.
        %        Note that DIFFORDER is vector valued; for example, the
        %        sequence
        %            u = treeVar([1 0 0], [0 1]);
        %            v = treeVar([0 1 0], [0 1]);
        %            w = treeVar([0 0 1], [0 1]);
        %            f = diff(u) + diff(w, 2);
        %        will lead to f.tree.DIFFORDER == [1 0 2].
        %    HEIGHT: The height of the syntax tree, i.e., the number of
        %        operations between the base variables(s) and the current
        %        variable.
        %    MULTCOEFF: The multiplication in front of the variable, which can
        %        either be a CHEBFUN or a scalar. For example, the sequence
        %            u = treeVar();
        %            v = sin(x)*u'l
        %        will have v.multcoeff == sin(x).
        %    ID: A Boolean vector, whose ith element is equal to 1 if the
        %        TREEVAR variable was constructed from the ith base variable,
        %        0 otherwise. For example, the sequence
        %            u = treeVar([1 0 0], [0 1]);
        %            v = treeVar([0 1 0], [0 1]);
        %            w = treeVar([0 0 1], [0 1]);
        %            f = u + 2*w;
        %        will lead to f.tree.ID == [1 0 1].
        %    HASTERMS: Indiciates whether a TREEVAR is constructed from a
        %        sequence of computations that include multiple terms. For
        %        example, the sequence
        %            u = treeVar(1, [0 1]);
        %            v = cos(x).*u;
        %            w = cos(u) + u;
        %        leads to v.tree.hasTerms = 0, w.tree.hasTerms = 1.
        tree
        method = 'constr';
        numArgs = 0;
        diffOrder
        height = 0;
        ID
        hasTerms = 0;
        left
        center
        right
        % The domain of the problem we're solving when constructing the
        % TREEVAR objects.
        domain
    end
    
    methods
        
        function obj = treeVar(IDvec, domain)
            % The TREEVAR constructor. See documentation above for calling
            % sequences to the constructor.
            
            if ( nargin > 0 )
                % Store the domain.
                obj.domain = domain;
            else
                % Default ID and domain.
                IDvec = 1;
                obj.domain = [-1 1];
            end
            
            obj.diffOrder = 0*IDvec;
            obj.ID = logical(IDvec);
            % Initialise a syntax tree for a base variable:
%             obj.tree  = struct('method', 'constr', 'numArgs', 0, ...
%                 'diffOrder', 0*IDvec, 'height', 0, 'multCoeff', 1, ...
%                 'ID', logical(IDvec),'hasTerms', 0);
        end
        
        function f = abs(f)
            f = univariate(f, 'abs');
        end
        
        function f = acos(f)
            f = univariate(f, 'acos');
        end
        
        function f = acosd(f)
            f = univariate(f, 'acosd');
        end
        
        function f = acot(f)
            f = univariate(f, 'acot');
        end
        
        function f = acoth(f)
            f = univariate(f, 'acoth');
        end
        
        function f = acsc(f)
            f = univariate(f, 'acsc');
        end
        
        function f = acscd(f)
            f = univariate(f, 'acscd');
        end
        
        function f = acsch(f)
            f = univariate(f, 'acsch');
        end
        
        function f = airy(f)
            f = univariate(f, 'airy');
        end
        
        function f = asec(f)
            f = univariate(f, 'asec');
        end
        
        function f = asecd(f)
            f = univariate(f, 'asecd');
        end
        
        function f = asech(f)
            f = univariate(f, 'asech');
        end
        
        function f = asin(f)
            f = univariate(f, 'asin');
        end
        
        function f = asind(f)
            f = univariate(f, 'asind');
        end
        
        function f = asinh(f)
            f = univariate(f, 'asinh');
        end
        
        function f = atan(f)
            f = univariate(f, 'atan');
        end
        
        function f = atand(f)
            f = univariate(f, 'atand');
        end
        
        function f = atanh(f)
            f = univariate(f, 'atanh');
        end
        
        function f = cos(f)
            f = univariate(f, 'cos');
        end
        
        function f = cosd(f)
            f = univariate(f, 'cosd');
        end
        
        function f = cosh(f)
            f = univariate(f, 'cosh');
        end
        
        function f = cot(f)
            f = univariate(f, 'cot');
        end
        
        function f = cotd(f)
            f = univariate(f, 'cotd');
        end
        
        function f = coth(f)
            f = univariate(f, 'coth');
        end
        
        function f = csc(f)
            f = univariate(f, 'csc');
        end
        
        function f = cscd(f)
            f = univariate(f, 'cscd');
        end
        
        function f = csch(f)
            f = univariate(f, 'csch');
        end
        
        function f = cumsum(f)
            %CUMSUM   Not supported.
            % We don't support integral equations with our first order
            % reformulation. However, we could accidentally end up here in case
            % of first order integral equation, where the conditions are
            % specified via N.LBC/RBC. Throw a meaningful error message in this
            % case.
            error('CHEBFUN:TREEVAR:cumsum:notSupported', ['First order ' ...
                'reformulation does not support integral equations.\nPlease ' ...
                'specify conditions via N.BC rather than N.LBC/RBC.'])
        end
        
        function f = diff(f, k)
            %DIFF   Derivative of a TREEVAR.
            
            % By default, compute first derivative:
            if ( nargin < 2 )
                k = 1;
            end
            
            % The derivative syntax tree.
            f.left = f;
            f.right = k;
            f.method = 'diff';
            f.diffOrder = f.diffOrder + k*f.ID;
            f.numArgs = 2;
%             f.tree = struct('method', 'diff', 'numArgs', 2, ...
%                 'left', f.tree, 'right', k, ...
%                 'diffOrder', f.tree.diffOrder + k*f.tree.ID, ...
%                 'height', f.tree.height + 1, ...
%                 'ID', f.tree.ID, ...
%                 'multCoeff', f.tree.multCoeff, ...
%                 'hasTerms', f.tree.hasTerms);
        end
        
        function disp(u)
            %DISP   Display a TREEVAR.
            
            if ( length(u) == 1 )
                % Scalar case.
%                 disp('treeVar with tree:')
%                 disp(u.tree);
%                 disp('and the domain:')
%                 disp(u.domain);
                    details(u)
            else
                % Systems case.
                disp('Array-valued treeVar, with trees:');
                for treeCounter = 1:length(u)
                    fprintf('tree %i\n', treeCounter)
                    disp(u(treeCounter).tree);
                end
                disp('and the domain:')
                disp(u(1).domain);
            end
        end
        
        function f = exp(f)
            f = univariate(f, 'exp');
        end
        
        function f = expm1(f)
            f = univariate(f, 'expm1');
        end
        
        function f = log(f)
            f = univariate(f, 'log');
        end
        
        function f = log10(f)
            f = univariate(f, 'log10');
        end
        
        function f = log2(f)
            f = univariate(f, 'log2');
        end
        
        function f = log1p(f)
            f = univariate(f, 'log1p');
        end
        
        function h = minus(f, g)
            %-   Subtraction of TREEVAR objects.
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                % (CHEBFUN/SCALAR) - TREEVAR
                h.tree = treeVar.bivariate(f, g.tree, 'minus', 1);
            elseif ( ~isa(g, 'treeVar') )
                % TREEVAR - (CHEBFUN/SCALAR)
                h.tree = treeVar.bivariate(f.tree, g, 'minus', 0);
            else
                % TREEVAR - TREEVAR
                h.tree = treeVar.bivariate(f.tree, g.tree, 'minus', 2);
            end
            h.domain = updateDomain(f, g);
        end
        
        function h = mrdivide(f, g)
            %/   Matrix division of TREEVAR objects
            %
            % This method only supports (SCALAR/TREEVAR)/(SCALAR/TREEVAR), i.e.
            % not (TREEVAR/CHEBFUN)/(TREEVAR/CHEBFUN).
            if ( isnumeric(f) || isnumeric(g) )
                h = rdivide(f, g);
                h.domain = updateDomain(f, g);
            else
                error('Dimension mismatch');
            end
        end
        
        
        function h = mtimes(f, g)
            %*   Matrix multiplication of TREEVAR objects.

            % This method only supports SCALAR/TREEVAR*SCALAR/TREEVAR, i.e. not
            % CHEBFUN/TREEVAR*CHEBFUN/TREEVAR.
            if ( isnumeric(f) || isnumeric(g) )
                h = times(f, g);
                h.domain = updateDomain(f, g);
            else
                error('Dimension mismatch');
            end
        end
        
        function f = pow2(f)
            f = univariate(f, 'pow2');
        end
        
        function h = power(f, g)
            %.^   Power of a TREEVAR.
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                % (CHEBFUN/SCALAR).^TREEVAR
                h.tree = treeVar.bivariate(f, g.tree, 'power', 1);
            elseif ( ~isa(g, 'treeVar') )
                % TREEVAR.^(CHEBFUN/SCALAR)
                h.tree = treeVar.bivariate(f.tree, g, 'power', 0);
            else
                % TREEVAR.^TREEVAR
                h.tree = treeVar.bivariate(f.tree, g.tree, 'power', 2);
            end
            h.domain = updateDomain(f, g);
        end
        
        function plot(treeVar)
            %PLOT   Plot of a TREEVAR syntax tree.
            %
            % See also TREEVAR.PLOTTREE.
            treeVar.plotTree(treeVar);
        end
        
        function h = plus(f, g)
            %+   Addition of TREEVAR objects.
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                % (CHEBFUN/SCALAR)+TREEVAR
                h.tree = treeVar.bivariate(f, g.tree, 'plus', 1);
            elseif ( ~isa(g, 'treeVar') )
                % TREEVAR + (CHEBFUN/SCALAR)
                h.tree = treeVar.bivariate(f.tree, g, 'plus', 0);
            else
                % TREEVAR + TREEVAR
                h.left = f;
                h.right = g;
                h.diffOrder = max(f.diffOrder, g.diffOrder);
                h.height = max(f.height, g.height) + 1;
                h.hasTerms = f.hasTerms || g.hasTerms;
                h.method = 'plus';
                h.numArgs = 2;
%         treeOut = struct('method', method, 'numArgs', 2, ...
%             'left', leftTree, 'right', rightTree, ...
%             'diffOrder', max(leftTree.diffOrder, rightTree.diffOrder), ...
%             'ID', leftTree.ID | rightTree.ID, ...
%             'height', max(leftTree.height, rightTree.height) + 1, ...
%             'hasTerms', isPM || leftTree.hasTerms || rightTree.hasTerms);

%                 h.tree = treeVar.bivariate(f.tree, g.tree, 'plus', 2);
            end
            h.domain = updateDomain(f, g);
        end
        
        function s = print(treeVar)
            %PRINT   Text rendering of a TREEVAR syntax tree.
            %
            % See also TREEVAR.PRINTTREE.
            s = treeVar.printTree(treeVar);
        end
        
        function h = rdivide(f, g)
            %./   Division of TREEVAR objects.
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                % (CHEBFUN/SCALAR)./TREEVAR
                h.tree = treeVar.bivariate(f, g.tree, 'rdivide', 1);
            elseif ( ~isa(g, 'treeVar') )
                % TREEVAR./(CHEBFUN/SCALAR)
                h.tree = treeVar.bivariate(f.tree, g, 'rdivide', 0);
            else
                % TREEVAR./TREEVAR
                h.tree = treeVar.bivariate(f.tree, g.tree, 'rdivide', 2);
            end
            h.domain = updateDomain(f, g);
        end
        
        function f = sec(f)
            f = univariate(f, 'sec');
        end
        
        function f = secd(f)
            f = univariate(f, 'secd');
        end
        
        function f = sech(f)
            f = univariate(f, 'sech');
        end

        function f = sin(f)
            f = univariate(f, 'sin');
%             f.center = f;
%             f.method = 'sin';
%             f.numArgs = 1;
%             f.height = f.height + 1;
        end
        
        function f = sind(f)
            f = univariate(f, 'sind');
        end
        
        function f = sinh(f)
            f = univariate(f, 'sinh');
        end

        function f = sqrt(f)
            f = univariate(f, 'sqrt');
        end
        
        function f = tan(f)
            f = univariate(f, 'tan');
        end
        
        function f = tand(f)
            f = univariate(f, 'tand');
        end
        
        function f = tanh(f)
            f = univariate(f, 'tanh');
        end
        
        function h = times(f, g)
            %.*   Multiplication of treeVar objects.
            
            % Initialize an empty TREEVAR
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                % (CHEBFUN/SCALAR).*TREEVAR
                h = g;
                h.method = 'times';
                h.numArgs = 2;
                h.left = f;
                h.right = g;
                h.height = h.height + 1;
                
%                 treeOut = struct('method', method, 'numArgs', 2, ...
%                     'left', leftTree, 'right', rightTree, ...
%                     'diffOrder', rightTree.diffOrder, ...
%                     'height', rightTree.height + 1, ...
%                     'ID', rightTree.ID, ...
%                     'hasTerms', isPM || rightTree.hasTerms);
%                 h.tree = treeVar.bivariate(f, g.tree, 'times', 1);
            elseif ( ~isa(g, 'treeVar') )
                % TREEVAR.*(CHEBFUN/SCALAR)
                h.tree = treeVar.bivariate(f.tree, g, 'times', 0);
            else
                % TREEVAR.*TREEVAR
                h.tree = treeVar.bivariate(f.tree, g.tree, 'times', 2);
            end
            h.domain = updateDomain(f, g);
        end
        
        function f = uminus(f)
            f = univariate(f, 'uminus');
        end
        
        function f = uplus(f)
            f = univariate(f, 'uplus');
        end
        
        % Construct syntax trees for univariate methods
        treeOut = univariate(treeIn, method)
    end
    
    methods ( Access = private )
        
        function dom = updateDomain(f, g)
            %UPDATEDOMAIN   Update domain in case we encounter new breakpoints.
            if ( isnumeric(f) )
                dom = g.domain;
            elseif ( isnumeric(g) )
                dom = f.domain;
            else
                dom = union(f.domain, g.domain);
            end
        end
        
    end
    
    methods ( Static = true )
        
        % Plot a syntax tree
        plotTree(tree, varargin)
        
        % Print a syntax tree
        s = printTree(tree, ind, indStr)
        
        % Returns how the results of evaluating BCs should be sorted
        idx = sortConditions(funIn, domain)
        
        % Convert higher order anonymous functions to first order systems
        [funOut, indexStart, problemDom, coeffs, totalDiffOrders] = ...
            toFirstOrder(funIn, rhs, domain)
        
    end
    
    methods ( Static = true, Access = private )
        
        % Construct syntax trees for bivariate methods
        treeOut = bivariate(leftTree, rightTree, method, type)
        
        % Convert expressions like 5*(diff(u) + u) to 5*diff(u) + 5*u
        newTree = expandTree(tree, maxOrder)
        
        % Split syntax trees into derivative part and non-derivative part
        [newTree, derTree] = splitTree(tree, maxOrder)
        
        % Convert the infix form of an expression to an anonymous function
        anonFun = toAnon(infix, varArray)
        
        % Convert infix expressions to anonymous function suited for ODE solvers
        funOut = toRHS(infix, varArray, coeff, indexStart, totalDiffOrders);
        
        % Convert a syntax tree to infix form
        [infix, varArray] = tree2infix(tree, diffOrders, varCounter, varArray)
        
    end

end
