function [sim, T] = model_NIMPLY_splicing(dose_ADintCZF1, dose_DsDedintN, z)
%Fig. 1H


%Retrieve parameters
[p, T] = Params;

%Inducible promoter
d.Rep1 = 21.1 / 200; %~4e9 gene copies

%Function for DsDed-ZF-mediated inhibition
function out = fDI1(A, DI)
    out = (p.b1 + p.w1E64 * A .* max(min(((((p.wr1E64 * DI) ...
        ./ (p.w1E64 * A + eps)) - p.l) * (1 - p.m1E64)) ...
        ./ (p.u - p.l) + p.m1E64, p.m1E64), 1)) ...
        ./ (1 + p.w1E64 * A + p.wr1E64 * DI);
end

%Initial values
IV = zeros(7, 1);


%Model
[~, sim] = ode15s(@(t, y, options)[ 
    
    % Y1 DsDed_intN RNA (INPUT)
      z(1)       * p.ktxEF1a * dose_DsDedintN ...
    - p.kdegR    * y(1)
    
    % Y2 DsDed_intN protein
      p.ktl      * y(1) ...
    - p.rec      * y(2) .* y(4) ...
    - p.kdegZFP  * y(2)
    
    % Y3 AD-intC-ZF1 RNA (INPUT)
      z(2)       * p.ktxEF1a * dose_ADintCZF1 ...     
    - p.kdegR    * y(3)
    
    % Y4 AD-intC-ZF1 protein
      p.ktl      * y(3) ...
    - p.rec      * y(2) .* y(4) ...
    - p.kdegintC * y(4)
    
    % Y5 DsDed-ZF1 reconstituted protein
      p.rec      * y(2) .* y(4) ...
    - p.kdegZFP  * y(5)  
    
    % Y6 Reporter RNA
      z(3)       * p.ktxZF * d.Rep1^0.5 * fDI1(y(4), y(5)) ...
    - p.kdegR    * y(6)
    
    % Y7 Reporter protein
      p.ktl      * y(6) ...
    - p.kdegRep  * y(7)
    
    ], T, IV);


end

