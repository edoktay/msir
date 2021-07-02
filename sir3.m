function x = sir3(A,b,precf,precw,precr,iter_max)
% SIR3   LU-based iterative refinement in three precisions.
%     x = sir3(A,b,precf,precw,precr, iter_max) solves Ax = b using LU-based
%     iterative refinement (with at most iter_max ref. steps), with
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
% Note: requires Cleve Laboratory and Advanpix multiprecision toolbox, and
% chop library (https://github.com/higham/chop)

if precf ~=0 && precf ~=1 && precf ~= 2, error('precf should be 0, 1 or 2'), end
if precw ~=0 && precw ~=1 && precw ~= 2, error('precw should be 0, 1 or 2'), end
if precr ~=1 && precr ~= 2 && precr ~= 4, error('precr should be 1, 2, or 4'), end

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
    fp.format = 'h';
    chop([],fp);
end

if precw == 0
    fprintf('**** Working precision is half.\n')
    fp.format = 'h';
    chop([],fp);
    u = float_params('h');
elseif precw == 2
    fprintf('**** Working precision is double.\n')
    A = double(A);
    b = double(b);
    u = eps('double');
else
    fprintf('**** Working precision is single.\n')
    A = single(A);
    b = single(b);
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

xact = double(mp(double(A),34)\mp(double(b),34));  % Exact solution
[~,~,~,xmax] = float_params(ufs);

% Compute LU factorization
if precf == 1
    [L,U,P] = lu(single(A));
    x =  U\(L\(P*single(b)) );
elseif precf == 2
    [L,U,P] = lu(double(A));
    x =  U\(L\(P*double(b)) );
else
    Ah = chop(A);
    [L,U,p] = lutx_chop(Ah);
    I = chop(eye(n)); P = I(p,:);
    t1 = lp_matvec(P,chop(b));
    t1 = trisol(L,t1);
    x = trisol(U,t1); 
end

% Note: when condition number of A in infinity norm is large, the initial 
% solution x can have 'Inf's in it. If so, default to using 0 as initial solution
if ( sum(isinf(single(x)))>0 || sum(isnan(single(x)))>0)
    
    % retry with scaling
    if (precf == 0)
        fprintf('**** Warning: x0 contains Inf or NaN. Retrying LU factorization with scaling.\n')
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
    end
    
    if ( sum(isinf(single(x)))>0 || sum(isnan(single(x)))>0)
        x =  zeros(size(b,1),1);
        fprintf('**** Warning: x0 contains Inf or NaN. Using 0 vector as initial solution.\n')
    end
    
end

% Store initial solution in working precision
if precw == 0
    x = chop(x);
elseif precw == 2
    x = double(x);
else
    x = single(x);
end

cged = false; iter = 0;

while ~cged
    
    % Compute size of errors, quantities in bound   
    ferr(iter+1) = double(norm(mp(double(x),34)-mp(xact,34),'inf')/norm(mp(xact,34),'inf'));
    mu(iter+1) = norm(double(A)*(mp(double(x),34)-mp(xact,34)),'inf')/(norm(mp(double(A),34),'inf')*norm(mp(double(x),34)-mp(xact,34),'inf')); 
    res = double(b) - double(A)*double(x);
    nbe(iter+1) = double(norm(mp(res,34),'inf')/(norm(mp(double(A),34),'inf')*norm(mp(double(x),34),'inf')+ norm(mp(double(b),34),'inf')));
    temp = double( abs(mp(res,34)) ./ (abs(mp(double(A),34))*abs(mp(double(x),34)) + abs(mp(double(b),34))) );
    temp(isnan(temp)) = 0; % Set 0/0 to 0.
    cbe(iter+1) = max(temp);
    
    iter = iter + 1;
    if iter > iter_max, break, end
    
    % Check convergence
    if max([ferr(iter) nbe(iter) cbe(iter)]) <= u, break, end
    
    % Compute residual vector
    if precr == 1
        rd = single(b) - single(A)*single(x);
    elseif precr == 2
        rd = double(b) - double(A)*double(x);
    else
        rd = mp(double(b),34) - mp(double(A),34)*mp(double(x),34);
    end
    
    % Scale residual vector
    norm_rd = norm(rd,'inf');
    rd1 = rd/norm_rd;
    
    % Solve for correction term
    if precf == 1
        d =  U\(L\(P*single(rd1)) );
    elseif precf == 2
        d =  U\(L\(P*double(rd1)) );
    else
        t1 = lp_matvec(P,rd1);
        t1 = trisol(L,t1);
        d = trisol(U,t1);
    end
    
    if precw == 1
        d = single(d);
    elseif precw == 2
        d = double(d);
    end
    
    xold = x;
    
    % Update solution
    if precw == 0
        x = chop(x + chop(chop(norm_rd)*chop(d)));
    elseif precw == 2
        x = x + norm_rd*double(d);
    else
        x = x + single(norm_rd)*single(d);
    end
    
    dx = norm(x-xold,'inf')/norm(x,'inf');
    
    % Check if dx contains infs, nans, or is 0
    if dx == Inf || isnan(double(dx))
        break;
    end
    
end

fprintf('\nNumber of iterations: %d\n',iter-1)

% Generate error plot
figure();
semilogy(0:iter-1, ferr, '-rx','LineWidth',1);
hold on
semilogy(0:iter-1, nbe, '-bo','LineWidth',1);
hold on
semilogy(0:iter-1, cbe, '-gv','LineWidth',1);
hold on
semilogy(0:iter-1, double(u)*ones(iter,1), '--k','LineWidth',1);

% Ensure only integers labeled on x axis
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
xlabel({'refinement step'},'Interpreter','latex','FontSize', 14);
tt = strcat('SIR'); 
title(tt,'Interpreter','latex','FontSize', 14);
h = legend('ferr','nbe','cbe');
set(h,'Interpreter','latex','FontSize', 14);

end
