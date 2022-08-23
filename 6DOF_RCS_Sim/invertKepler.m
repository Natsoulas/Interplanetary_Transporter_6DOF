function [E,counter] = invertKepler(M, e, E0)
 %invKepler - Newton-Raphson inversion of the Kepler Equation
 %      invKepler(M,e) calculates the eccentric anomalies (E)
 %      corresponding to mean anomalies (M) and eccentricites (e).
 %
 %      INPUTS
 %      M    n x 1   mean anomalies (in radians)
 %      e    n x 1   eccentricity values. if scalar, same eccentricity is
 %                   assumed for all values of M
 %      E0   n x 1   initial condition for iteration
 %      
 %      OUTPUTS
 %      E        n x 1 vector of eccentric anomalies (in radians)
 %      counter  n x 1 # of iterations taken 
 
 %condition inputs
 M = M(:);
 e = e(:);
 
 %initialize vars
 tol = eps(2*pi);
 counter = ones(size(M));
 maxcounter = 1000;
 del = ones(size(M));
 
 if numel(e) == 1 && e == 0
     E = M;
 else
     %initialize E
     if ~exist('E0','var') || isempty(E0)
         E = M./(1-e);
         inds = E > sqrt(6*(1-e)./e);
         if length(e) == 1
             einds = 1;
         else
             einds = inds;
         end
         E(inds) = (6*M(inds)./e(einds)).^(1/3);
     else
         E = E0;
     end
     
     %iterate
     while ((max(del) > tol) && (max(counter) < maxcounter))
         inds = del > tol;
         if numel(e) == 1
             E(inds) = E(inds) - (M(inds) - E(inds) + e*sin(E(inds)))./...
                 (e*cos(E(inds))-1);
             del(inds) = abs(M(inds) - (E(inds) - e*sin(E(inds))));
         else
             E(inds) = E(inds) - (M(inds) - E(inds) + ...
                 e(inds).*sin(E(inds)))./(e(inds).*cos(E(inds))-1);
             del(inds) = abs(M(inds) - (E(inds) - e(inds).*sin(E(inds))));
         end
         counter(inds) = counter(inds)+1;
     end
     if (max(counter) == maxcounter)
         warning('invKepler:overIteration',...
             'Maximum number of iterations exceeded');
     end
 end