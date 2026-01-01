% Shape properties of the cell. Units: [meters]
Par.Inv.Length = 2.95 * 1e-6;
Par.Inv.Width  = 1.07 * 1e-6;
Par.Inv.Radius = Par.Inv.Width / 2;
% 
% Par.Inv.S = 2 * pi * Par.Inv.Radius^2 + 2 * pi * Par.Inv.Radius * Par.Inv.Length; % Surface area [meters^2]
% Par.Inv.V = pi * Par.Inv.Radius^2 * Par.Inv.Length;                               % Volume       [meters^3]

% Updated definition of surface and volume (for a spherocylinder)
Model.Cell.S   = pi * Par.Inv.Width * Par.Inv.Length;
Model.Cell.Vol = pi * (Par.Inv.Width / 2)^2 * (Par.Inv.Length - Par.Inv.Width / 3);
 
% Homeostatic targets
Model.Cell.pHi_Target  = 7.8;
Model.Cell.DPi_Target  = 0;
Model.Cell.Cost_Target = 0; 
Model.Cell.PMF_Target  = -0.170;

% Weight cost
Model.Cell.pHi_weight      = 1;
Model.Cell.DPi_weight      = 1;
% Model.Cell.Antiport_weight = 0;

% DV
Model.Cell.DV_Inf = [];     Model.Cell.DV_Sup = [];     Model.Cell.DV = [];

% DPi
Model.Cell.Cost_Inf = 0;    Model.Cell.Cost_Sup = 2e10 * 3.3;
%Model.Cell.DPi_Inf = [];     Model.Cell.DPi_Sup = [];

Model.Inv.F  = 96485;
Model.Inv.R  = 8.31;
Model.Inv.NA = 6.02 * 1e23;

Model.Inv.T = 298;

% Membrane specific capacitance. Units: [Farad/meters^2]
Model.Inv.Cm = 6.5e-3;

% pKa
Model.Inv.pKb_COH = 0.2;
Model.Inv.pKa_AH  = -6.3;
Model.Inv.pKa_w   = 14;

