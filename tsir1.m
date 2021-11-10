function [x,cged,ferr,nbe,cbe,sirit,gmres_midits,gmresits,switch_iter_mid,switch_iter,switch_iter_last] = tsir1(A,b,precf,precw,precr,iter_max,rho_thresh,x,xact)

n = length(A);

if precf == 1
  %  fprintf('**** Factorization precision is single.\n')
    ufs = 'single';
elseif precf == 2
  %  fprintf('**** Factorization precision is double.\n')
    ufs = 'double';
else
  %  fprintf('**** Factorization precision is half.\n')
    ufs = 'half';
end

if precw == 0
  %  fprintf('**** Working precision is half.\n')
    uws = 'half';
    uws1 = 'h';
    fp.format = 'h';
    chop([],fp);
    gtol=1e-2;
    u = float_params('h');
elseif precw == 2
  %  fprintf('**** Working precision is double.\n')
    uws = 'double';
    A = double(A);
    b = double(b);
    gtol = 1e-10;
    u = eps('double');
else
  %  fprintf('**** Working precision is single.\n')
    uws = 'single';
    A = single(A);
    b = single(b);
    gtol=1e-6;
    u = eps('single');
end

if precr == 1
  %  fprintf('**** Residual precision is single.\n')
    urs = 'single';
elseif precr == 2
  %  fprintf('**** Residual precision is double.\n')
    urs = 'double';
else
  %  fprintf('**** Residual precision is quad.\n')
    urs = 'quad';
    mp.Digits(34);
end

[uh,xmins,xmin,xmax] = float_params(ufs);



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
        LL = P'*L; %(double(P')*double(L));
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

%Compute size of initial errors
% ferri = double(norm(mp(double(x),34)-mp(xact,34),'inf')/norm(mp(xact,34),'inf'));
% res = double(b) - double(A)*double(x);
% nbei = double(norm(mp(res,34),'inf')/(norm(mp(double(A),34),'inf')*norm(mp(double(x),34),'inf')+ norm(mp(double(b),34),'inf')));
% temp = double( abs(mp(res,34)) ./ (abs(mp(double(A),34))*abs(mp(double(x),34)) + abs(mp(double(b),34))) );
% temp(isnan(temp)) = 0; % Set 0/0 to 0.
% cbei = max(temp);


% Initialization
cged = 0;

gmres_midits = []; gmresits = [];
gmres_miderr = []; gmreserr = [];

sir = 1; gmresir_mid = 0; gmresir_tsir = 0; gmresir_mid_dummy_tsir = 0;

gmresmid_maxiter = round(0.1*n);

dex = u^(-1);

switch_iter_mid = 0; switch_iter = 0; switch_iter_last = 0;

rho_threshmax = 0;
ferrs = []; nbes = []; cbes = [];
ferrsg = []; nbesg = []; cbesg = [];
ferrg = []; nbeg = []; cbeg = [];

%Run SIR
[x,cged,switch_iter_mid, ferrs, nbes, cbes] = tsir_sir(A,b,precf,precw,precr,iter_max, L,U,P,x, xact, rho_thresh,u);

%If SIR didn't converge, run SGMRESIR
if ~cged
    gmresir_mid_dummy_tsir = 1;
    [x,cged,switch_iter1,gmres_midits,gmres_miderr, ferrsg, nbesg, cbesg] = tsir_sgmresir(A,b,precf,precw,precr,iter_max,LL,U,x, xact, rho_thresh,gtol,u,gmresmid_maxiter);
    switch_iter = switch_iter_mid + switch_iter1;
end

%If SGMRESIR didn't converge, run GMRESIR
if ~cged
    gmresir_tsir = 1;
    [x,cged,switch_iter2,gmresits,gmreserr, ferrg, nbeg, cbeg] = tsir_gmresir1(A,b,precf,precw,precr,iter_max,LL,U,x, xact, rho_thresh,gtol,u,gmresmid_maxiter);
    switch_iter_last = switch_iter + switch_iter2;
end

%Concatenate error vectors for plotting
% ferr = [ferri, ferrs, ferrsg, ferrg];
% nbe = [nbei, nbes, nbesg, nbeg];
% cbe = [cbei, cbes, cbesg, cbeg];
ferr = [ferrs, ferrsg, ferrg];
nbe = [nbes, nbesg, nbeg];
cbe = [cbes, cbesg, cbeg];

sirit = num2str(numel(ferrs));

end