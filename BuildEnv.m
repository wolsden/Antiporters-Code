function Model = BuildEnv(Model)

% Concentration of protons [mM]
Model.H.c_e  = 10.^(-Model.Env.pHe) * 1e3;

% COncentration of hydroxides [mM]
Model.OH.c_e = 10.^(Model.Env.pHe - Model.Inv.pKa_w) * 1e3;

% Concentration of ions outside [mM]
% Here it does not really matter because I ve always assumed (a bit of a stringent assumption
% that just aims at simplifying the problem) in my thesis
% there are no "captive molecules" outside, so the total concentration
% outside is identical to the total concentrations of small ions outside.
% Within this scope you may think of Model.Env.Pie_ion at total
% concentrations coming from addition of [NaCl] in the medium.
Model.Env.Pie_ion = Model.Env.Pie;

% Here it is describe in the draft SI, i compute the amount of eg [Na]e,
% [Cl]e depending on the pHe sought to set.
for i = 1 : length(Model.Env.pHe)
    % if pHe sought to be maintained is alkaline then 
    if Model.Env.pHe(i) > 7
        % compute the amount of "free" Cl (ie not complexed under the form
        % HCl)
        Model.A.c_e(i) = 1/2 * Model.Env.Pie_ion(i) / (1 + 10^(Model.Inv.pKa_AH - Model.Env.pHe(i)));
        % Compute the amount of [Na]e; indeed to reach alkaline pHe, we add
        % some NaOH, so [Na]e may deviate signifcantly from the amount of
        % [Na] we add to the medium via NaCl addition.
        Model.C.c_e(i) = Model.A.c_e(i) + Model.OH.c_e(i) - Model.H.c_e(i);
    % Similar to above, but in this time in the acidic pHe range (we add
    % HCl)
    elseif Model.Env.pHe(i) < 7
        Model.C.c_e(i) = 1/2 * Model.Env.Pie_ion(i) / (1 + 10^(Model.Inv.pKb_COH + Model.Env.pHe(i) - Model.Inv.pKa_w));
        Model.A.c_e(i) = Model.C.c_e(i) - Model.OH.c_e(i) + Model.H.c_e(i);
    % Simple case when pHe is neutral, we add neither HCl or NaOH, so their
    % concentrations are only set by the amount of NaCl we add (here called
    % Pie_ion)
    elseif Model.Env.pHe(i) == 7
        Model.C.c_e(i) = 1/2 * Model.Env.Pie_ion(i);
        Model.A.c_e(i) = 1/2 * Model.Env.Pie_ion(i);
    end 
end

% Here just compute some gratuitous variables we want to remember the
% values for latter.
Model.H.c_e  = Model.H.c_e;
Model.OH.c_e = Model.OH.c_e;
Model.C.c_e  = Model.C.c_e;
Model.A.c_e  = Model.A.c_e;

% fractional ionic composition
%Model.NDC.alpha_e  = Model.NDC.c_e  ./ Model.Env.Pie';
%Model.NDNC.alpha_e = Model.NDNC.c_e ./ Model.Env.Pie';
Model.H.alpha_e    = Model.H.c_e    ./ Model.Env.Pie';
Model.OH.alpha_e   = Model.OH.c_e   ./ Model.Env.Pie';
Model.C.alpha_e    = Model.C.c_e    ./ Model.Env.Pie';
Model.A.alpha_e    = Model.A.c_e    ./ Model.Env.Pie';

end