function [Jx Jt] = JacobianCalc
%%% Symbolically calculate the jacobian of the poisson equations
%%% for implementation into C++

n = 3; % 3 Points to simulate, BC's + 1 intermediate point
nc = 7; % Number of PDEs

% Initialise arrays
x = sym('x', [n nc]); 
t = sym('t', [n nc]); 

f = x;
for i=1:n*nc
    f(i) = 0;
end

% Constants
C.z = sym('z', [3 1]);
C.beta = sym('beta', [3 1]);
C.alpha = sym('alpha', [3 1]);
C.b = sym('Cb', [3 1]);
C.gammax = sym('gammax', [1 1]);

% Index positions
C.V = 1;
C.C = [2 3 4];
C.A = [5 6 7];
C.nC = size(C.C, 2);
C.nA = size(C.A, 2);

syms gammax dx Kw Ka F_CONST R_CONST T_CONST
syms V1 V2

% Electrode reaction terms
R.i0 = sym('i0');
R.E0 = sym('E0');
R.beta = sym('betaBV');
R.n = sym('n_e');
R.n_const = sym('n_const');
R.bv_const = sym('bv_const');

% Make 3rd component diffuse only
C.z(3) = 0;
C.beta(3) = 0;

%% Poisson
temp = 0;
for i = 1:C.nC
	temp = temp + C.z(i) * x(1, C.C(i));
end
f(1, C.V) = V1 - 2*x(1, C.V) + x(2, C.V) - C.gammax * temp;

temp = 0;
for i = 1:C.nC
	temp = temp + C.z(i) * x(2, C.C(i));
end
f(2, C.V) = x(1, C.V) - 2*x(2, C.V) + x(3, C.V) - C.gammax * temp;

temp = 0;
for i = 1:C.nC
	temp = temp + C.z(i) * x(3, C.C(i));
end
f(3, C.V) = x(2, C.V) - 2*x(3, C.V) + V2 - C.gammax * temp;

%% NP - diffusion, migration
i = 1;
alpha = C.alpha(i) * (-x(1, C.C(i)) + x(2, C.C(i)));
beta = C.beta(i) * ( ...
    x(1, C.C(i)) * (V1 - 2*x(1, C.V) + x(2, C.V))  ...
    -x(2, C.C(i)) * (x(1, C.V) - x(2, C.V)));
f(1, C.C(i)) = t(1, C.C(i)) - alpha - beta;

alpha = C.alpha(i) * (x(1, C.C(i)) - 2*x(2, C.C(i)) + x(3, C.C(i)));
beta = C.beta(i) * ( ...
    x(1, C.C(i)) * (x(1, C.V) - x(2, C.V))  ...
    +x(2, C.C(i)) * (x(1, C.V) - 2*x(2, C.V) + x(3, C.V))  ...
    -x(3, C.C(i)) * (x(2, C.V) - x(3, C.V)));
f(2, C.C(i)) = t(2, C.C(i)) - alpha - beta;

alpha = C.alpha(i) * (x(2, C.C(i)) - 2*x(3, C.C(i)) + C.b(i));
beta = C.beta(i) * ( ...
    x(2, C.C(i)) * (x(2, C.V) - x(3, C.V))  ...
    +x(3, C.C(i)) * (x(2, C.V) - 2*x(3, C.V) + V2)  ...
    -C.b(i) * (x(3, C.V) - V2));
f(3, C.C(i)) = t(3, C.C(i)) - alpha - beta;

%% NP - diffusion, migration, butler volmer
i = 2;

alpha = C.alpha(i) * (-x(1, C.C(i)) + x(2, C.C(i)));
beta = C.beta(i) * ( ...
    x(1, C.C(i)) * (V1 - 2*x(1, C.V) + x(2, C.V))  ...
    -x(2, C.C(i)) * (x(1, C.V) - x(2, C.V)));
E = R.E0 - R.n_const * log(x(1,C.C(i))/x(1,C.C(i+1)));
BV = R.i0 * (exp(-R.beta*R.bv_const*(E-V1)) ...
	- exp((1-R.beta)*R.bv_const) * (E-V1));
f(1, C.C(i)) = t(1, C.C(i)) - alpha - beta - BV / R.n / F_CONST / dx;

alpha = C.alpha(i) * (x(1, C.C(i)) - 2*x(2, C.C(i)) + x(3, C.C(i)));
beta = C.beta(i) * ( ...
    x(1, C.C(i)) * (x(1, C.V) - x(2, C.V))  ...
    +x(2, C.C(i)) * (x(1, C.V) - 2*x(2, C.V) + x(3, C.V))  ...
    -x(3, C.C(i)) * (x(2, C.V) - x(3, C.V)));
f(2, C.C(i)) = t(2, C.C(i)) - alpha - beta;

alpha = C.alpha(i) * (x(2, C.C(i)) - 2*x(3, C.C(i)) + C.b(i));
beta = C.beta(i) * ( ...
    x(2, C.C(i)) * (x(2, C.V) - x(3, C.V))  ...
    +x(3, C.C(i)) * (x(2, C.V) - 2*x(3, C.V) + V2)  ...
    -C.b(i) * (x(3, C.V) - V2));
f(3, C.C(i)) = t(3, C.C(i)) - alpha - beta;

%% NP - diffusion, butler volmer
i = 3;
alpha = C.alpha(i) * (-x(1, C.C(i)) + x(2, C.C(i)));
E = R.E0 - R.n_const * log(x(1,C.C(i-1))/x(1,C.C(i)));
BV = R.i0 * (exp(-R.beta*R.bv_const*(E-V1)) ...
	- exp((1-R.beta)*R.bv_const) * (E-V1));
f(1, C.C(i)) = t(1, C.C(i)) - alpha + BV / R.n / F_CONST / dx;
   
alpha = C.alpha(i) * (x(1, C.C(i)) - 2*x(2, C.C(i)) + x(3, C.C(i)));
f(2, C.C(i)) = t(2, C.C(i)) - alpha;

alpha = C.alpha(i) * (x(2, C.C(i)) - 2*x(3, C.C(i)) + C.b(i));
f(3, C.C(i)) = t(3, C.C(i)) - alpha;

%% Algebraics
for i = 1:n
	f(i, C.A(1)) = x(i, C.A(1)) - x(i, C.A(2)) - x(i, C.A(3));
	f(i, C.A(2)) = x(i, C.A(2)) - Kw / x(i, C.A(1));
	f(i, C.A(3)) = x(i, C.A(3)) - (Ka * x(i, C.C(3))) / (Ka + x(i, C.A(1)));
end

%% Reorganise for jacobian
ft = transpose(f); xt = transpose(x); tt = transpose(t);
ft = [ft(:,1); ft(:,2); ft(:,3)];
xt = [xt(:,1); xt(:,2); xt(:,3)];
tt = [tt(:,1); tt(:,2); tt(:,3)];

Jx = jacobian(ft, xt);
Jt = jacobian(ft, tt);
end
