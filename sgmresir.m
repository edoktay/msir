function x = sgmresir(A,b,precf,precw,precr,iter_max,savename)
%SGMRESIR   GMRES-BASED iterative refinement in three precisions, uses uniform precision within GMRES.
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
    % fprintf('**** Factorization precision is single.\n')
    ufs = 'single';
elseif precf == 2
    % fprintf('**** Factorization precision is double.\n')
    ufs = 'double';
else
    % fprintf('**** Factorization precision is half.\n')
    ufs = 'half';
end

if precw == 0
    % fprintf('**** Working precision is half.\n')
    uws = 'half';
    uws1 = 'h';
    fp.format = 'h';
    chop([],fp);
    gtol=1e-2;
    u = float_params('h');
elseif precw == 2
    % fprintf('**** Working precision is double.\n')
    uws = 'double';
    A = double(A);
    b = double(b);
    gtol = 1e-10;
    u = eps('double');
else
    % fprintf('**** Working precision is single.\n')
    uws = 'single';
    A = single(A);
    b = single(b);
    gtol=1e-6;
    u = eps('single');
end

if precr == 1
    % fprintf('**** Residual precision is single.\n')
    urs = 'single';
elseif precr == 2
    % fprintf('**** Residual precision is double.\n')
    urs = 'double';
else
    % fprintf('**** Residual precision is quad.\n')
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
        LL = P'*L;
        LL = (1/mu)*diag(1./diag(R))*(LL);
        U = (U)*diag(1./diag(C));
        x = U\(LL\b);
    else
        t1 = lp_matvec(P,chop(b));
        t1 = trisol(L,t1);
        x = trisol(U,t1);
    end
    
end

%Compute condition number of A, cond(A), and
%cond(A,x) for the exact solution
kinfA = cond(mp(double(A),34),'inf');
condAx = norm(abs(inv(mp(double(A),34)))*abs(mp(double(A),34))*abs(xact),inf)/norm(xact,inf);
condA = norm(abs(inv(mp(double(A),34)))*abs(mp(double(A),34)),'inf');


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


%Initialization
cged = 0;

iter = 0; dx = 0; rd = 0;

gmres_midits = [];
gmres_miderr = [];


while ~cged
    
    %Compute size of errors, quantities in bound
    ferr(iter+1) = double(norm(mp(double(x),34)-mp(xact,34),'inf')/norm(mp(xact,34),'inf'));
    res = double(b) - double(A)*double(x);
    nbe(iter+1) = double(norm(mp(res,34),'inf')/(norm(mp(double(A),34),'inf')*norm(mp(double(x),34),'inf')+ norm(mp(double(b),34),'inf')));
    temp = double( abs(mp(res,34)) ./ (abs(mp(double(A),34))*abs(mp(double(x),34)) + abs(mp(double(b),34))) );
    temp(isnan(temp)) = 0; % Set 0/0 to 0.
    cbe(iter+1) = max(temp);
    
    iter = iter + 1;
    if iter > iter_max, break, end
    
    %Check convergence
    if isnan(ferr(iter)) || isnan(nbe(iter)) || isnan(cbe(iter))
        cged = 0;
        break;
    end
    if max([ferr(iter) nbe(iter) cbe(iter)]) <= u
        cged = 1;
        break;
    end
    
    
    %Compute residual vector
    if precr == 1
        rd = single(b) - single(A)*single(x);
    elseif precr == 2
        rd = double(b) - double(A)*double(x);
    else
        rd = mp(double(b),34) - mp(double(A),34)*mp(double(x),34);
    end
    
    %Scale residual vector
    norm_rd = norm(rd,'inf');
    rd1 = rd/norm_rd;
    
    
    %Call GMRES to solve for correction term
    if precw == 0
        [d, err, its, ~] = gmres_hh( A, chop(zeros(n,1)), chop(rd1), LL, U, n, 1, gtol);
    elseif precw == 2
        [d, err, its, ~] = gmres_dd( A, zeros(n,1), double(rd1), LL, U, n, 1, gtol);
    else
        [d, err, its, ~] = gmres_ss( A, single(zeros(n,1)), single(rd1), LL, U, n, 1, gtol);
    end
    
    if precw == 1
        d = single(d);
    elseif precw == 2
        d = double(d);
    end
    
    
    %Record number of iterations gmres took
    gmresmid_lastiter = its;
    gmres_midits = [gmres_midits,its];
    
    %Record final relative (preconditioned) residual norm in GMRES
    gmres_miderr = [gmres_miderr,err(end)];
    
    %Record relative (preconditioned) residual norm in each iteration of
    %GMRES (so we can look at convergence trajectories if need be)
    gmres_miderrvec{iter} = err;
    
    
    
    %Update solution
    if precw == 0
        x = chop(x + chop(chop(norm_rd)*chop(d)));
    elseif precw == 2
        x = x + norm_rd*double(d);
    else
        x = x + single(norm_rd)*single(d);
    end

    
    
end

%Construct output with number of steps and iterations for printing
output = '(';
for k = 1:numel(gmres_midits)
    if k == 1
        output = strcat(output, num2str(gmres_midits(k)));
    else
        output = strcat(output, ', ', num2str(gmres_midits(k)));
    end
end
output = strcat(output, ')');

if cged == 0
    output = '-';
end

fprintf('SGMRES-IR Iteration/Step Count: %s\n',output)


%Generate plot

%plot ferr, nbe, cbe
fig1 = figure();
semilogy(0:iter-1, ferr, '-rx');
hold on
semilogy(0:iter-1, nbe, '-bo');
hold on
semilogy(0:iter-1, cbe, '-gv');
hold on
semilogy(0:iter-1, double(u)*ones(iter,1), '--k');



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


tt = strcat('SGMRES-IR');
title(tt,'Interpreter','latex');

h = legend('ferr','nbe','cbe');
set(h,'Interpreter','latex');

%Save figure as pdf
if ~isempty(savename)
    saveas(gcf, strcat(savename,'.pdf'));
end

end
