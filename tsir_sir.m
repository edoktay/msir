function [x,cged,switchit, ferr, nbe, cbe] = tsir_sir(A,b,precf,precw,precr,iter_max, L,U,P,x, xact, rho_thresh,u)
%TSIR_SIR   LU-based iterative refinement in three precisions used within the TSIR function.
%     Solves Ax = b using LU-based
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
% Note: requires Cleve Laboratory, Advanpix multiprecision toolbox, and
% chop library (https://github.com/higham/chop)

n = length(A);

%Initialization
x0 = x;
cged = 0;

ferr = [];
nbe = [];
cbe = [];
phi = [];
switchit = 0;

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
    
    %Check if correction has nans or infs; if so, return without converging
    if ( sum(isinf(single(d)))>0 || sum(isnan(single(d)))>0)
        switchit = 0;
        cged = 0;
        return;
    else      
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
                break;
            end
            
        end
        
        %Compute quantities needed to decide whether to switch to SGMRESIR
        norm_ddex = norm(d,'inf')/norm(dex,'inf');
        rho_threshmax = max(rho_threshmax, norm_ddex);
        norm_dx = norm_rd*norm(d,'inf')/norm(x,'inf');
        norm_dx = double(norm_dx);
        phi(i) = norm_dx/(1-rho_threshmax);
        dex = d;
        
        %Check whether we should switch to SGMRESIR
        if ( (norm_dx <= u) || (norm_ddex >= rho_thresh) ||  phi(i) < u )
            
            %Convergence detected, but we will keep iterating for now
            if ( (phi(i) >= 0) && (phi(i) <= u) )
                fprintf('\n SIR Convergence Detected\n');
            else
                %If the error is larger than the initial error we started
                %with, reset initial solution
                if ( phi(i) > phi(1) )
                    x = x0;
                    phi(i) = phi(1);
                    switchit = i;
                    cged = 0;
                    return;                  
                else
                    switchit = i;
                    cged = 0;
                    return;
                end
                
            end
        end
        
    end 
end


end