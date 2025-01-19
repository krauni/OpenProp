% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% ================================ savedata ===============================
function savedata(hObject,ED)
    %
    % This subfunction saves all values presented in the OpenProp GUI.
    %

    global OpenPropDirectory SpecificationsValues DuctValues FlagValues FoilValues Filename...
           XR_in XCoD_in XCD_in VAI_in VTI_in ri_in Xt0oD_in skew0_in rake0_in...
           Meanline_cell Thickness_cell XCoD_values XCLmax_values; % CavValues

    %%


    filename   	= get(Filename,'string');                       % Filename prefix

    % -------------------------------------------------------------------------
    % Figure out what the current working directory is,
    % and change directories to OpenPropDirectory/filename
    rest = pwd;

    while ~isempty(rest)
        [CurrentDirectory,rest] = strtok(rest,'/');

        if strcmp(CurrentDirectory,OpenPropDirectory)

            if strcmp(rest(2:end),filename)
                % already in /OpenPropDirectory/filename/
                rest = [];

                addpath ../SourceCode

            elseif strcmp(rest(2:end),'SourceCode')

                mkdir(['../',filename])
                   cd(['../',filename])
                rest = [];

                addpath ../SourceCode

            elseif isempty(rest)

                % you are in /OpenPropDirectory/
                mkdir(['./',filename])
                   cd(['./',filename])

                addpath ../SourceCode

            else
                % you are in /OpenPropDirectory/wrongfolder
                disp('ERROR1: Must start OpenProp from the root directory.')
                return
            end
        end
    end
    % -------------------------------------------------------------------------


    Z           = str2double(get(SpecificationsValues(1),'string'));  % number of blades
    N           = str2double(get(SpecificationsValues(2),'string'));  % propeller speed [RPM]
    D           = str2double(get(SpecificationsValues(3),'string'));	% propeller diameter [m]
    THRUST      = str2double(get(SpecificationsValues(4),'string')); 	% required thrust [N]
    Vs          = str2double(get(SpecificationsValues(5),'string'));  % ship velocity [m/s]
    Dhub        = str2double(get(SpecificationsValues(6),'string'));  % hub diameter [m]
    rho         = str2double(get(SpecificationsValues(7),'string')); 	% water density [kg/m^3]
    Mp          = str2double(get(SpecificationsValues(8),'string')); 	% number of vortex panels over the radius
    Np          = str2double(get(SpecificationsValues(9),'string')); 	% Number of points over the chord [ ]

    TAU         = str2double(get(DuctValues(1),'string'));      % Thrust ratio
    CDd         = str2double(get(DuctValues(2),'string'));      % Duct section drag coefficient

    % H           = str2double(get(CavValues(1),'string'));       % Shaft centerline depth [m]
    % dV          = str2double(get(CavValues(2),'string'));       % Inflow variation [m/s]


    Propeller_flag	= get(FlagValues(1),'value');               % 0 == turbine, 1 == propeller
    Hub_flag  	= get(FlagValues(3),'value');                   % 0 == no hub, 1 == hub
    Duct_flag	= get(FlagValues(4),'value');                   % 0 == no duct, 1 == duct

    Chord_flag	= get(FlagValues(5),'value');                   % ** CHORD OPTIMIZATION FLAG **

    Viscous_flag	= get(FlagValues(6),'value');               % 0 == viscous forces off (CD = 0), 1 == viscous forces on
    Plot_flag       = get(FlagValues(7),'value');               % 0 == do not display plots, 1 == display plots

    Make2Dplot_flag = get(FlagValues(8),'value');               % 0 == do not make a 2D plot of the results, 1 == make plot
    Make3Dplot_flag = get(FlagValues(8),'value');               % 0 == do not make a 3D plot of the results, 1 == make plot

    Analyze_flag	= get(FlagValues(9),'value');

    % Cav_flag	= get(FlagValues(10),'value');                   % 0 == do not run cavitation mapping, 1 == run cavitation mapping



    Meanline_index	= get(FoilValues(1),'value');
    Meanline        = char(Meanline_cell(Meanline_index));       	% Meanline form

    Thickness_index	= get(FoilValues(2),'value');
    Thickness       = char(Thickness_cell(Thickness_index));    	% Thickness form


    XR        = str2double(get(XR_in,'string'));                % radius / propeller radius
    XCoD      = str2double(get(XCoD_in,'string'));              % chord / diameter
    XCD       = str2double(get(XCD_in,'string'));               % section drag coefficient

    ri        = str2double(get(ri_in, 'string'));
    VAI       = str2double(get(VAI_in,'string'));               % axial      inflow velocity / ship velocity
    VTI       = str2double(get(VTI_in,'string'));               % tangential inflow velocity / ship velocity

    % f0oc0     = str2double(get(f0oc_in,'String'));              % max section camber    / chord
    Xt0oD     = str2double(get(Xt0oD_in,'String'));             % max section thickness / chord
    skew0     = str2double(get(skew0_in,'String'));             % skew
    rake0     = str2double(get(rake0_in,'String'));             % rake


    save([filename,'_GUIsd'],'Z','N','D','THRUST','Vs','Dhub','rho','Mp','Np',...
        'TAU','CDd','Propeller_flag','Hub_flag','Duct_flag',...
        'Chord_flag','Viscous_flag','Plot_flag','Make2Dplot_flag','Make3Dplot_flag','Analyze_flag',...
        'Meanline','Meanline_index','Thickness',...
        'Thickness_index','filename','XR','XCoD','XCD','ri','VAI','VTI','Xt0oD','skew0',...
        'rake0','XCoD_values','XCLmax_values'); % 'Cav_flag','H','dV',

end

