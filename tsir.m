function x = tsir(A,b,precf,precw,precr,iter_max,rho_thresh,savename)
%TSIR   Three-stage iterative refinement in three precisions.
%     x = tsir3(A,b,precf,precw,precr,iter_max,rho_thresh) solves Ax = b using 
%     iterative refinement (with at most 3*iter_max ref. steps), with
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

if precf ~=0 && precf ~=1 && precf ~= 2, error('precf should be 0, 1 or 2'), end
if precw ~=0 && precw ~=1 && precw ~= 2, error('precw should be 0, 1 or 2'), end
if precr ~=1 && precr ~= 2 && precr ~= 4, error('precr should be 1, 2, or 4'), end

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

xact = double(mp(double(A),34)\mp(double(b),34));

[uh,xmins,xmin,xmax] = float_params(ufs);



%Compute LU factorization
if precf == 1
    [L,U,P] = lu(single(A));
    LL = single(double(P')*double(L));
    x =  U\(L\(P*single(b)) );
elseif precf == 2
    [L,U,P] = lu(double(A));
    LL = double(double(P')*double(L));
    x =  U\(L\(P*double(b)) );
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
        x = U\(LL\b);
    else
        t1 = lp_matvec(P,chop(b));
        t1 = trisol(L,t1);
        x = trisol(U,t1);
    end
    
end

%Note: when kinf(A) is large, the initial solution x can have 'Inf's in it
%If so, default to using 0 as initial solution
if ( sum(isinf(single(x)))>0 || sum(isnan(single(x)))>0)
    x =  zeros(size(b,1),1);
    fprintf('**** Warning: x0 contains Inf or NaN. Using 0 vector as initial solution.\n')
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
ferri = double(norm(mp(double(x),34)-mp(xact,34),'inf')/norm(mp(xact,34),'inf'));
res = double(b) - double(A)*double(x);
nbei = double(norm(mp(res,34),'inf')/(norm(mp(double(A),34),'inf')*norm(mp(double(x),34),'inf')+ norm(mp(double(b),34),'inf')));
temp = double( abs(mp(res,34)) ./ (abs(mp(double(A),34))*abs(mp(double(x),34)) + abs(mp(double(b),34))) );
temp(isnan(temp)) = 0; % Set 0/0 to 0.
cbei = max(temp);


% Initialization
cged = 0;

gmres_midits = []; gmresits = [];
gmres_miderr = []; gmreserr = [];

sir = 1; gmresir_mid = 0; gmresir = 0; gmresir_mid_dummy = 0;

gmresmid_maxiter = 0.1*n;

dex = u^(-1);

switch_iter_mid = 0;

rho_threshmax = 0;
ferrs = []; nbes = []; cbes = [];
ferrsg = []; nbesg = []; cbesg = [];
ferrg = []; nbeg = []; cbeg = [];

%Run SIR
[x,cged,switch_iter_mid, ferrs, nbes, cbes] = tsir_sir(A,b,precf,precw,precr,iter_max, L,U,P,x, xact, rho_thresh,u);

%If SIR didn't converge, run SGMRESIR
if ~cged
    gmresir_mid_dummy = 1;
    [x,cged,switch_iter1,gmres_midits,gmres_miderr, ferrsg, nbesg, cbesg] = tsir_sgmresir(A,b,precf,precw,precr,iter_max,LL,U,x, xact, rho_thresh,gtol,u,gmresmid_maxiter);
    switch_iter = switch_iter_mid + switch_iter1;
    
end

%If SGMRESIR didn't converge, run GMRESIR
if ~cged
    gmresir = 1;
    [x,cged,i,gmresits,gmreserr, ferrg, nbeg, cbeg] = tsir_gmresir(A,b,precf,precw,precr,iter_max,LL,U,x, xact, rho_thresh,gtol,u);
end


%Concatenate error vectors for plotting
ferr = [ferri, ferrs, ferrsg, ferrg];
nbe = [nbei, nbes, nbesg, nbeg];
cbe = [cbei, cbes, cbesg, cbeg];

%Construct output with number of steps and iterations for printing
sirit = num2str(numel(ferrs));
output = sirit;

if ~isempty(gmres_midits)
    
    output = strcat(output, ', (');
    for k = 1:numel(gmres_midits)
        if k == 1
            output = strcat(output, num2str(gmres_midits(k)));
        else
            output = strcat(output, ', ', num2str(gmres_midits(k)));
        end
    end
    output = strcat(output, ')');
end

if ~isempty(gmresits)
    output = strcat(output, ', (');
    for k = 1:numel(gmresits)
        if k == 1
            output = strcat(output, num2str(gmresits(k)));
        else
            output = strcat(output, ', ', num2str(gmresits(k)));
        end
    end
    output = strcat(output, ')');
end

if cged == 0
    output = '-';
end

fprintf('TSIR Iteration/Step Count: %s\n',output)


%Generate plot

%plot ferr, nbe, cbe
fig1 = figure();
semilogy(0:numel(ferr)-1, ferr, '-rx');
hold on
semilogy(0:numel(nbe)-1, nbe, '-bo');
hold on
semilogy(0:numel(cbe)-1, cbe, '-gv');
hold on
semilogy(0:numel(nbe)-1, double(u)*ones(numel(nbe),1), '--k');


if gmresir_mid_dummy
    if isempty(gmres_midits)==0
        hold on
        semilogy(switch_iter_mid,ferr(switch_iter_mid+1),'p','MarkerSize',14,'MarkerFaceColor','magenta','MarkerEdgeColor',[0,0,0]);
        hold on
        semilogy(switch_iter_mid,nbe(switch_iter_mid+1),'p','MarkerSize',14,'MarkerFaceColor','magenta','MarkerEdgeColor',[0,0,0]);
        hold on
        semilogy(switch_iter_mid,cbe(switch_iter_mid+1),'p','MarkerSize',14,'MarkerFaceColor','magenta','MarkerEdgeColor',[0,0,0]);
    end
end
if gmresir
    if isempty(gmresits)==0
        hold on
        semilogy(switch_iter,ferr(switch_iter+1),'p','MarkerSize',14,'MarkerFaceColor','yellow','MarkerEdgeColor',[0,0,0]);
        hold on
        semilogy(switch_iter,nbe(switch_iter+1),'p','MarkerSize',14,'MarkerFaceColor','yellow','MarkerEdgeColor',[0,0,0]);
        hold on
        semilogy(switch_iter,cbe(switch_iter+1),'p','MarkerSize',14,'MarkerFaceColor','yellow','MarkerEdgeColor',[0,0,0]);
    end
end


%Ensure only integers labeled on x axis
atm = get(gca,'xticklabels');
m = str2double(atm);
xlab = [];
num = 1;
for i = 1:numel(m)
    if ceil(m(i)) == m(i)
        xlab(num) = m(i);
        num = num + 1;
    end
end
set(gca,'xticklabels',xlab);
set(gca,'xtick',xlab);
xlabel({'refinement step'},'Interpreter','latex');
set(gca,'FontSize',14)
a = get(gca,'Children');
set(a,'LineWidth',1);
set(a,'MarkerSize',10);

tt = strcat('TSIR');
title(tt,'Interpreter','latex');

h = legend('ferr','nbe','cbe');
set(h,'Interpreter','latex');

%Save figure as pdf
if ~isempty(savename)
    saveas(gcf, strcat(savename,'.pdf'));
end

end