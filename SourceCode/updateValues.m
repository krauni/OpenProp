% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function updateValues(hObject,ED)

%%

% Update values displayed in the Calculator panel

global SpecificationsValues Values;

N           = str2double(get(SpecificationsValues(2),'string'));  % propeller speed [RPM]
D           = str2double(get(SpecificationsValues(3),'string'));	% propeller diameter [m]
THRUST      = str2double(get(SpecificationsValues(4),'string')); 	% required thrust [N]
Vs          = str2double(get(SpecificationsValues(5),'string'));  % ship velocity [m/s]
Dhub        = str2double(get(SpecificationsValues(6),'string'));  % hub diameter [m]
rho         = str2double(get(SpecificationsValues(7),'string')); 	% water density [kg/m^3]

n       = N/60;    % ** propeller speed [rev/s] = Vs/(Js*D) = N/60
lambda 	= n*2*pi*(D/2)/Vs;
Js      = Vs/(n*D);              	% ** Js = Vs/(n*D) ,  advance coefficient
KT      = THRUST/(rho*n^2*D^4);          % ** KT = THRUST/(rho*n^2*D^4)

CT      = THRUST/(0.5*rho*Vs^2*pi*(D/2)^2);

set(Values(1),'string',num2str(Js));
set(Values(2),'string',num2str(lambda));
set(Values(3),'string',num2str(KT));
set(Values(4),'string',num2str(CT));

end
