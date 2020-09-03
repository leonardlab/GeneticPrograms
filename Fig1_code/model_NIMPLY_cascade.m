function [sim, T] = model_NIMPLY_cascade(dose_DsDedZF10, dose_DsDedZF1, z)
%Fig. 1I


%Retrieve parameters
[p, T] = Params;

%Plasmid ng-equivalents of constitutive genes
d.ZF10a     = 21.1;
d.DsDedZF10 = 21.1;
d.ZF1a      = 2.11;

%Inducible promoters
d.Rep1         = 21.1 / 200; %~4e9 gene copies
d.DsDedZF1PEST = 21.1 / 200; %~4e9 gene copies

%Functions for DsDed-ZF-mediated inhibition
function out = f10(A, DI)
    out = (p.b10 + p.w10E64 * A .* max(min(((((p.wr10E64 * DI) ...
        ./ (p.w10E64 * A + eps)) - p.l) * (1 - p.m10E64)) ...
        ./ (p.u - p.l) + p.m10E64, p.m10E64), 1)) ...
        ./ (1 + p.w10E64 * A + p.wr10E64 * DI);
end
function out = f1(A, DI)
    out = (p.b1 + p.w1E64 * A .* max(min(((((p.wr1E64 * DI) ...
        ./ (p.w1E64 * A + eps)) - p.l) * (1 - p.m1E64)) ...
        ./ (p.u - p.l) + p.m1E64, p.m1E64), 1)) ...
        ./ (1 + p.w1E64 * A + p.wr1E64 * DI);
end

%Initial values
IV = zeros(12, 1);


%Model
[~, sim] = ode15s(@(t, y, options)[ 
    
    % Y1 ZF10a RNA 
      z(1)      * p.ktxEF1a * d.ZF10a ...
    - p.kdegR    * y(1)
    
    % Y2 ZF10a protein
      p.ktl      * y(1) ...
    - p.kdegZFP  * y(2)
    
    % Y3 DsDed-ZF10 RNA (INPUT)
      z(2)      * p.ktxEF1a * dose_DsDedZF10 ...     
    - p.kdegR    * y(3)
    
    % Y4 DsDed-ZF10 protein
      p.ktl      * y(3) ...
    - p.kdegZFP  * y(4)
    
    % Y5 ZF1a RNA
      z(3)       * p.ktxEF1a * d.ZF1a ...     
    - p.kdegR    * y(5)
    
    % Y6 ZF1a protein
      p.ktl      * y(5) ...
    - p.kdegZFP  * y(6)
    
    % Y7 DsDed-ZF1-PEST RNA
      z(4)       * p.ktxZF * d.DsDedZF1PEST^0.5 * f10(y(2), y(4)) ...
    - p.kdegR    * y(7)
    
    % Y8 DsDed-ZF1-PEST protein
      p.ktl          * y(7) ...
    - p.kdegZFP_PEST * y(8)
    
    % Y9 DsDed-ZF1 RNA (INPUT)
      z(5)      * p.ktxEF1a * dose_DsDedZF1 ...     
    - p.kdegR   * y(9)
    
    % Y10 DsDed-ZF1 protein
      p.ktl     * y(9) ...
    - p.kdegZFP * y(10)
    
    % Y11 Reporter RNA
      z(7)      * p.ktxZF * d.Rep1^0.5 * f1(y(6), y(8) + y(10)) ...
    - p.kdegR   * y(11)
    
    % Y12 Reporter protein
      p.ktl     * y(11) ...
    - p.kdegRep * y(12)
    
    ], T, IV);


end

