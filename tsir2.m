function [x,cged,ferr,nbe,cbe,sirit,gmres_midits,gmresits,switch_iter_mid,switch_iter] = tsir2(A,b,precf,precw,precr,iter_max,rho_thresh,x,xact,u_old)
%TSIR2   Three-stage iterative refinement in three precisions used within the MSIR function.
%     [x,cged,ferr,nbe,cbe,sirit,gmres_midits,gmresits,switch_iter_mid,switch_iter] = tsir2(A,b,precf,precw,precr,iter_max,rho_thresh,x,xact,u_old) 
%     solves Ax = b using iterative refinement (with at most 3*iter_max ref. steps), with
%     LU factors computed in precision precf:
%       * half if precf = 0,
%       * single if precf = 1,
%       * double if precf = 2,
%     working precision precw:
%       * half if precw = 0,
%       * single if precw = 1,
%       * double if precw = 2,
%     and residuals computed at precision precr:
%       * single if precr = 1,
%       * double if precr = 2,
%       * quad if precr = 4
%
% Note: requires Cleve Laboratory, Advanpix multiprecision toolbox, and
% chop library (https://github.com/higham/chop)

n = length(A);

if precf == 1
   fprintf('**** Factorization precision is single.\n')
   ufs = 'single';
elseif precf == 2
   fprintf('**** Factorization precision is double.\n')
   ufs = 'double';
else
   fprintf('**** Factorization precision is half.\n')
   ufs = 'half';
end

if precw == 0
   fprintf('**** Working precision is half.\n')
   fp.format = 'h';
   chop([],fp);
   gtol=1e-2;
   u = float_params('h');
elseif precw == 2
   fprintf('**** Working precision is double.\n')
   A = double(A);
   b = double(b);
   gtol = 1e-10;
   u = eps('double');
else
    fprintf('**** Working precision is single.\n')
    A = single(A);
    b = single(b);
    gtol=1e-6;
    u = eps('single');
end

if precr == 1
   fprintf('**** Residual precision is single.\n')
elseif precr == 2
   fprintf('**** Residual precision is double.\n')
else
   fprintf('**** Residual precision is quad.\n')
   mp.Digits(34);
end

[~,~,~,xmax] = float_params(ufs);

%Compute LU factorization
if precf == 1
    [L,U,P] = lu(single(A));
    LL = single(double(P')*double(L));
elseif precf == 2
    [L,U,P] = lu(double(A));
    LL = double(double(P')*double(L));
else
    Ah = chop(A);
    [L,U,p] = lutx_chop(Ah);
    
    I = chop(eye(n)); P = I(p,:);
    LL = P'*L;

    % If L or U factors contain NAN or INF, try again with scaling
    if ( sum(sum(isinf(single(LL))))>0 || sum(sum(isnan(single(LL))))>0 || sum(sum(isinf(single(U))))>0 || sum(sum(isnan(single(U))))>0 )
        
        [Ah,R,C] = scale_diag_2side(A);
        mu = (0.1)*xmax;
        Ah = mu*Ah;
        
        Ah = chop(Ah);
        [L,U,p] = lutx_chop(Ah);
        I = chop(eye(n)); P = I(p,:);
        LL = P'*L; 
        LL = (1/mu)*diag(1./diag(R))*(LL);
        U = (U)*diag(1./diag(C));
    else
        t1 = lp_matvec(P,chop(b));
        t1 = trisol(L,t1);
    end
    
end

%Store initial solution in working precision
if precw == 0
    x = chop(x);
elseif precw == 2
    x = double(x);
else
    x = single(x);
end

% Initialization
cged = 0;

gmres_midits = []; gmresits = [];
gmres_miderr = []; gmreserr = [];

sir = 1; gmresir_mid = 0; gmresir_tsir = 0; gmresir_mid_dummy_tsir = 0;

gmresmid_maxiter = round(0.1*n);

dex = u^(-1);

switch_iter_mid = 0; switch_iter = 0;

ferrs = []; nbes = []; cbes = [];
ferrsg = []; nbesg = []; cbesg = [];
ferrg = []; nbeg = []; cbeg = [];

%Run SIR
[x,cged,switch_iter_mid, ferrs, nbes, cbes] = tsir_sir(A,b,precf,precw,precr,iter_max, L,U,P,x, xact, rho_thresh,u,u_old);

%If SIR didn't converge, run SGMRESIR
if ~cged
    gmresir_mid_dummy_tsir = 1;
    [x,cged,switch_iter1,gmres_midits,gmres_miderr, ferrsg, nbesg, cbesg] = tsir_sgmresir(A,b,precf,precw,precr,iter_max,LL,U,x, xact, rho_thresh,gtol,u,gmresmid_maxiter,u_old);
    switch_iter = switch_iter_mid + switch_iter1;
end

%If SGMRESIR didn't converge, run GMRESIR
if ~cged
    gmresir_tsir = 1;
    [x,cged,i,gmresits,gmreserr, ferrg, nbeg, cbeg] = tsir_gmresir(A,b,precf,precw,precr,iter_max,LL,U,x, xact, rho_thresh,gtol,u);
end

%Concatenate error vectors for plotting
ferr = [ferrs, ferrsg, ferrg];
nbe = [nbes, nbesg, nbeg];
cbe = [cbes, cbesg, cbeg];

sirit = num2str(numel(ferrs));

end