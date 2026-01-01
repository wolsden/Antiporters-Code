% clear workspace variable
clear all

% file name under which we ll save simulation data
% SaveName = 'Data/Figure1BCD.mat';

% Define parameters
BuildInvariants;
BuildCell;

% Defne some coarse grained parameter
eta = Model.Inv.F / Model.Inv.R / Model.Inv.T;

% Define the range of PMF/pHe that we try
psiH_vec = linspace(-0.24 , 0  , 101);
pHe_vec  = linspace(2  , 12 , 101);

% Define the range of captive charge we try 
% (alphaY := z_Y * [Y]i / Pie where Pie is total concentration of extracellular ions
zY_vec = [-0.1 0 0.1];

% Set total extracellular concentration of small ions (eg. := [Na]e + [Cl]e)
Model.Env.Pie = 200;

alphaY_vec = zY_vec/Model.Env.Pie;

% Set intracellular pH
Model.pHi = 7;

% start time counter
tic;

% Set two counters
counter = 0;
countertest = 0;
% Set loop over alphaY values (3 values in this case)
for i = 1 : length(alphaY_vec)
    alphaY = alphaY_vec(i);
    
    % Loop over PMF values
    for j = 1 : length(psiH_vec)
        psiH = psiH_vec(j);
        
        % Loop over pHe values
        for k = 1 : length(pHe_vec)
            Model.Env.pHe = pHe_vec(k);
            
            % Update the "countertest"
            countertest = countertest + 1;
            
            % Build environmental composition from pHe, Pie, upon some
            % assumption; see model for the extracellular environment in SI
            Model = BuildEnv(Model);
            
            % Here define the equation we wish to numerically solve for the limit Inf. of voltage
            % (assuming only cations are pumped, and their DG is such that
            % DG_cation -> -infinity (no constraint on the PMF yet)
            % since anions are not pumped: DG_anions = 0
            syms DVmin
            EQ0 = DVmin * Model.Inv.Cm * Model.Cell.S / Model.Inv.F / Model.Cell.Vol / Model.Env.Pie == ...
            alphaY - Model.A.alpha_e * exp(eta * DVmin);
            % Solve numerically the equation
            SolDVmin = vpasolve(EQ0);
            % MATLAB subtelty, here just change the type of the solution
            % from "sym type" to "double type"
            DV_Inf = double(SolDVmin);
            
            % Compute voltage necessary to maintain pHi = 7, given current
            % value of PMF and pHe
            DV = psiH - log(10) / eta * (Model.Env.pHe - Model.pHi);
            
            % here just record the value of the current loop
            Sim.alphaY(i)         = alphaY;
            Sim.pHe(k)            = Model.Env.pHe;
            Sim.psiH(j)           = psiH;
            
            % if voltage necessary is more than the inferior limit on
            % voltage (first accessibility criterion we did not talk about
            % last time)
            if DV > DV_Inf
            
                % then solve for the cationic motive force necessary to
                % achieve the voltage necessary to maintain pHi = 7 given
                % the current {PMF,pHe}
                
                % here define the equation to solve
                syms psiC
                EQ = DV * Model.Inv.Cm * Model.Cell.S / Model.Inv.F / Model.Cell.Vol / Model.Env.Pie == ...
                    alphaY + Model.C.alpha_e * exp(eta * (psiC - DV)) - Model.A.alpha_e * exp(eta * DV);
            
                % Solve the equation for Cationic motive force
                SolpsiC = vpasolve(EQ);
                % change the type to double
                Sim.psiC(j , k , i)   = double(SolpsiC);
                
                % Record the current voltage
                Sim.DV(j , k , i)     = DV;
            
            else
                % else if voltage required to maintain pHi given current
                % {PMF,pHe} is below the minimum voltage attainable
                % (assuming CMF goes to -infinity, effectively when
                % additioning the constraints we talked about last time the
                % effective inferior limit on the CMF will be finite,
                % typically in the range -500 mV to 0 mV.
                
                % Set voltage and CMF to NaN (non computable values, there
                % are no solutions).
                Sim.psiC(j , k , i)   = NaN;
                Sim.DV(j , k , i)     = NaN;
                
            end

            % Here just a praticality to tell me how much computational
            % time has elapsed as well as the number of simulations that
            % have been ran so far.
            if mod(countertest , 1e3) == 0
                sprintf('%0.3f , counter=%i' , toc , counter)
            end
            
        end
    end
end
