clear all;
BuildInvariants;
BuildCell;

% 1. Define the Pie values you want to iterate over
Pie_vec = [200, 100, 20, 2, 0.2, 0.02, 0.002];

% Define other static parameters
psiH_vec = linspace(-0.24 , 0  , 101);
pHe_vec  = [5.5, 7.5];
zY_vec   = [-10 0 10];
Model.pHi = 7.5;
eta = Model.Inv.F / Model.Inv.R / Model.Inv.T;

% Define symbolic variables ONCE outside the loops to save time
syms DVmin psiC

% 2. START THE OUTER LOOP
for m = 1:length(Pie_vec)
    % Update the current Pie value
    Model.Env.Pie = Pie_vec(m);
    alphaY_vec = zY_vec / Model.Env.Pie;
    
    fprintf('Starting simulation for Pie = %g...\n', Model.Env.Pie);
    tic;
    
    % Re-initialize Sim structure for each Pie value
    Sim = struct();
    countertest = 0;

    for i = 1 : length(alphaY_vec)
        alphaY = alphaY_vec(i);
        
        for j = 1 : length(psiH_vec)
            psiH = psiH_vec(j);
            
            for k = 1 : length(pHe_vec)
                Model.Env.pHe = pHe_vec(k);
                countertest = countertest + 1;
                
                % Build environment
                Model = BuildEnv(Model);
                
                % Solve for DVmin
                EQ0 = DVmin * Model.Inv.Cm * Model.Cell.S / Model.Inv.F / Model.Cell.Vol / Model.Env.Pie == ...
                      alphaY - Model.A.alpha_e * exp(eta * DVmin);
                SolDVmin = vpasolve(EQ0, DVmin);
                DV_Inf = double(SolDVmin);
                
                % Compute required voltage
                DV = psiH - log(10) / eta * (Model.Env.pHe - Model.pHi);
                
                % Record values
                Sim.alphaY(i) = alphaY;
                Sim.pHe(k)    = Model.Env.pHe;
                Sim.psiH(j)   = psiH;
                
                if DV > DV_Inf
                    % Solve for psiC
                    EQ = DV * Model.Inv.Cm * Model.Cell.S / Model.Inv.F / Model.Cell.Vol / Model.Env.Pie == ...
                         alphaY + Model.C.alpha_e * exp(eta * (psiC - DV)) - Model.A.alpha_e * exp(eta * DV);
                
                    SolpsiC = vpasolve(EQ, psiC);
                    Sim.psiC(j, k, i) = double(SolpsiC);
                    Sim.DV(j, k, i)   = DV;
                else
                    Sim.psiC(j, k, i) = NaN;
                    Sim.DV(j, k, i)   = NaN;
                end
                
                % Print progress every 500 iterations
                if mod(countertest, 500) == 0
                    fprintf('Pie: %g | Time: %0.2fs | Iter: %i\n', Model.Env.Pie, toc, countertest);
                end
            end
        end
    end
    
    % 3. SAVE DATA with unique filename
    filename = ['Sim_', num2str(Model.Env.Pie, '%g'), '_updated.mat'];
    save(filename, 'Sim');
    fprintf('Finished and saved: %s\n', filename);
end