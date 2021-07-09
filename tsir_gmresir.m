function [x,cged,i,gmresits,gmreserr,ferr,nbe,cbe] = tsir_gmresir(A,b,precf,precw,precr,iter_max,LL,U,x, xact, rho_thresh,gtol,u)
%TSIR_GMRESIR   GMRES-based iterative refinement in three precisions used within the TSIR function.
%     Solves Ax = b using GMRES-based
%     iterative refinement where the preconditioner is applied in double the working precision
%     (with at most iter_max ref. steps), with
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

%Initialization
x0 = x;
cged = 0;

gmresits = [];
gmreserr = [];

ferr = [];
nbe = [];
cbe = [];
phi = [];

dex = u^(-1);
rho_threshmax = 0;

for i = 1:iter_max
    
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
    
    
    %Solve for correction term
    %Call GMRES to solve for correction term
    if precw == 0
        [d, err, its, ~] = gmres_hs( A, chop(zeros(n,1)), chop(rd1), LL, U, n, 1, gtol);
    elseif precw == 2
        [d, err, its, ~] = gmres_dq( A, zeros(n,1), double(rd1), LL, U, n, 1, gtol);
    else
        [d, err, its, ~] = gmres_sd( A, single(zeros(n,1)), single(rd1), LL, U, n, 1, gtol);
    end
    
    %Record number of iterations gmres took
    gmresits = [gmresits,its];
    
    %Record final relative (preconditioned) residual norm in GMRES
    gmreserr = [gmreserr,err(end)];
    
    %Record relative (preconditioned) residual norm in each iteration of
    %GMRES (so we can look at convergence trajectories if need be)
    gmreserrvec{i} = err;
    
    if precw == 1
        d = single(d);
    elseif precw == 2
        d = double(d);
    end
    
    
    if i < iter_max
        %Update solution
        if precw == 0
            x = chop(x + chop(chop(norm_rd)*chop(d)));
        elseif precw == 2
            x = x + norm_rd*double(d);
        else
            x = x + single(norm_rd)*single(d);
        end
        
        %Compute size of errors
        ferr(i) = double(norm(mp(double(x),34)-mp(xact,34),'inf')/norm(mp(xact,34),'inf'));
        res = double(b) - double(A)*double(x);
        nbe(i) = double(norm(mp(res,34),'inf')/(norm(mp(double(A),34),'inf')*norm(mp(double(x),34),'inf')+ norm(mp(double(b),34),'inf')));
        temp = double( abs(mp(res,34)) ./ (abs(mp(double(A),34))*abs(mp(double(x),34)) + abs(mp(double(b),34))) );
        temp(isnan(temp)) = 0; % Set 0/0 to 0.
        cbe(i) = max(temp);
        
        %Check convergence
        if isnan(ferr(i)) || isnan(nbe(i)) || isnan(cbe(i))
            cged = 0;
            break;
        end
        if max([ferr(i) nbe(i) cbe(i)]) <= u
            cged = 1;
            return;
        end
    end
    
    
    
end



end