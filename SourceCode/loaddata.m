% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% ================================ loaddata ===============================
function loaddata(hObject,ED)
    %
    % This subfunction loads all values presented in the OpenProp GUI.
    %

    global OpenPropDirectory SpecificationsValues DuctValues FlagValues FoilValues Filename...
           XR_in XCoD_in XCD_in VAI_in VTI_in ri_in Xt0oD_in skew0_in rake0_in; % CavValues

    %%

    % -------------------------------------------------------------------------
    % Figure out what the current working directory is,
    % and change directories to /OpenPropDirectory/
    rest = pwd;

    while ~isempty(rest)
        [CurrentDirectory,rest] = strtok(rest,'/');
    end

    if     strcmp(CurrentDirectory,OpenPropDirectory)    % OpenPropDirectory == 'OpenProp_vX.X.X';

        % stay in /OpenPropDirectory/
        uiload;
        cd(['./',filename])

    elseif strcmp(CurrentDirectory,'SourceCode')

        cd('../')
        uiload();
        cd(['./',filename])

    else
        % already in /OpenPropDirectory/filename
        uiload;
    end
    % -------------------------------------------------------------------------


    % ------------------------------------------------------------------------
    % ------------------------------------------------------------------------

    set(SpecificationsValues(1),'string',num2str(Z));         % number of blades
    set(SpecificationsValues(2),'string',num2str(N));         % propeller speed [RPM]
    set(SpecificationsValues(3),'string',num2str(D));         % propeller diameter [m]
    set(SpecificationsValues(4),'string',num2str(THRUST)); 	% required thrust [N]
    set(SpecificationsValues(5),'string',num2str(Vs));        % ship velocity [m/s]
    set(SpecificationsValues(6),'string',num2str(Dhub));      % hub diameter [m]
    set(SpecificationsValues(7),'string',num2str(rho));       % water density [kg/m^3]
    set(SpecificationsValues(8),'string',num2str(Mp));        % number of vortex panels over the radius
    set(SpecificationsValues(9),'string',num2str(Np));        % Number of points over the chord [ ]

    set(DuctValues(1),'string',num2str(TAU));           % Thrust ratio
    set(DuctValues(2),'string',num2str(CDd));           % Duct section drag coefficient

    % set(CavValues(1),'string',num2str(H));              % Shaft centerline depth [m]
    % set(CavValues(2),'string',num2str(dV));             % Inflow variation [m/s]

    set(FlagValues(1),'value',Propeller_flag);       	% 0 == turbine, 1 == propeller
    set(FlagValues(3),'value',Hub_flag);              	% 0 == no hub, 1 == hub
    set(FlagValues(4),'value',Duct_flag);              	% 0 == no duct, 1 == duct


    set(FlagValues(5),'value',Chord_flag);              % ** CHORD OPTIMIZATION FLAG **
    changeChord;       % Testing whether running the callback function on loading updates fields and avoids loading issue


    set(FlagValues(6),'value',Viscous_flag);          	% 0 == viscous forces off (CD = 0), 1 == viscous forces on
    set(FlagValues(7),'value',Plot_flag);         % 0 == do not display plots, 1 == display plots
    set(FlagValues(8),'value',Make3Dplot_flag);        	% 0 == do not display plots, 1 == display plots

    set(FlagValues(9),'value',Analyze_flag);

    % set(FlagValues(10),'value',Cav_flag);               	% 0 == do not run cavitation mapping, 1 == run cavitation mapping


    set(FoilValues(1),'value',Meanline_index);          	% Meanline form
    set(FoilValues(2),'value',Thickness_index);        	% Thickness form

    set(Filename,'string',filename);                	% Filename prefix

    % ----------------------------------------------
    % --- Loop to set new values for input table ---

    for index = 1 : length(XR);

        set(XR_in(index),'string',num2str(XR(index)));                    % radius / propeller radius
        set(XCoD_in(index),'string',num2str(XCoD(index)));                % chord / diameter
        set(XCD_in(index),'string',num2str(XCD(index)));                  % section drag coefficient

        % set(f0oc_in(index),'String',num2str(f0oc0(index)));             % max section camber    / chord
        set(Xt0oD_in(index),'String',num2str(Xt0oD(index)));              % max section thickness / chord
        set(skew0_in(index),'String',num2str(skew0(index)));              % skew
        set(rake0_in(index),'String',num2str(rake0(index)));              % rake

    end

    for index = 1 : length(ri)
        if isnan(ri(index))
            set(ri_in(index),'string','');
        else
            set(ri_in(index),'string',num2str(ri(index)));
        end

        if isnan(VAI(index))
            set(VAI_in(index),'string','');
        else
            set(VAI_in(index),'string',num2str(VAI(index)));
        end

        if isnan(VTI(index))
            set(VTI_in(index),'string','');
        else
            set(VTI_in(index),'string',num2str(VTI(index)));
        end
    end
    % ----------------------------------------------

end

