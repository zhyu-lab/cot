function [converged, decrease] = em_converged_m(loglik, previous_loglik, verbose1,threshold, check_increased)

if nargin < 3, verbose1 = 1; end
if nargin < 4, threshold = 1e-4; end
if nargin < 5, check_increased = 1; end

converged = 0;
decrease = 0;

if check_increased
    if loglik - previous_loglik < -1e-3 % allow for a little imprecision
        if verbose1 ==1
            fprintf(1, '******likelihood decreased from %6.4f to %6.4f!\n', previous_loglik, loglik);
        end
        decrease = 1;
    end
end

delta_loglik = abs(loglik - previous_loglik);
avg_loglik = (abs(loglik) + abs(previous_loglik) + eps)/2;
if (delta_loglik / avg_loglik) < threshold, converged = 1; end