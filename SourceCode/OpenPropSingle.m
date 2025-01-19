% ----------------------------------------------------------------------- %
%                                                                         %
%                              0111000                                    %
%                           100 1 100 001                                 %
%                         10    1  1  1 00                                %
%                        01  1  1  1      0                               %
%                       0 1  1  1   1  1 1 0                              %
%                       0   1   1   1  1  1 0                             %
%                       0 1     1   1  1  1 0                             %
%                       0 1  1  1   1  0  1 0                             %
%                       0 1  1  1   0  1    0                             %
%                       01 1        1  1 1 0                              %
%                        0    0  1  0 1   0                               %
%                         0         1    0                                %
%                    10010 0 1101111110111                                %
%                  10 1 1  1111111111 11 11                               %
%                 0 1 1 1 11111111101011010111                            %
%                01 11    11111111 1  1    1 110                          %
%               011    1 1 111111110011  1 1 1 110                        %
%               0   11 1 1 1 111      0  1 1 1   10                       %
%               0 1   11  1  0         1 1 1 1 1 1 0                      %
%               1  11 1 1   11          0  1 1 1 1 11                     %
%                0     1 1  0           011  1 1 1 10                     %
%                10 1   1  0             0  1 1 1  11                     %
%                 10     01               01      10                      %
%                   10001                   001 100                       %
%                                             111                         %
%                                                                         %
%             ____                   _____                                %
%            / __ \                 |  __ \                               %
%           | |  | |_ __   ___ _ __ | |__) | __ ___  _ __                 %
%           | |  | | '_ \ / _ \ '_ \|  ___/ '__/ _ \| '_ \                %
%           | |__| | |_) |  __/ | | | |   | | | (_) | |_) |               %
%            \____/| .__/ \___|_| |_|_|   |_|  \___/| .__/                %
%                  | |                              | |                   %
%                  |_|                              |_|                   %
%                                                                         %
%             An integrated rotor design and analysis tool.               %
%                                                                         %
%                                                                         %
% Copyright (C) 2011, Brenden Epps.                                       %
%                                                                         %
% This program is free software; you can redistribute it and/or modify it %
% under the terms of the GNU General Public License as published by the   %
% Free Software Foundation.                                               %
%                                                                         %
% This program is distributed in the hope that it will be useful, but     %
% WITHOUT ANY WARRANTY; without even the implied warranty of              %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    %
% See the GNU General Public License for more details.                    %
%                                                                         %
% ----------------------------------------------------------------------- %


% =========================================================================
% OpenProp_v3.3.4
% Last modified: 10/21/2013 Brenden Epps
% =========================================================================
%
% This code runs the single-design GUI.
%
%--------------------------------------------------------------------------

% =========================================================================
% ========================= Initiate OpenProp ======================
function OpenPropSingle

    clear all;
    clear global;

    % warning off
    addpath ../SourceCode
    % addpath  ./SourceCode


    % =========================================================================
    % ------------------------- Setup global variables ------------------------

    global Fig_Main;    % Main GUI figure

    global Select;

    global OpenPropDirectory SpecificationsValues DuctValues FlagValues FoilValues Filename...
           XR_in XCoD_in XCD_in VAI_in VTI_in ri_in Xt0oD_in skew0_in rake0_in...
           Meanline_cell Thickness_cell Values; % CavValues


    % ------ Variable definitions ------
    %
    % Variables marked with ** are not taken from user input through the GUI.
    % Instead, these variables are computed through formulas or remain unused
    % in this version of the GUI but are still present to allocate space for
    % future use.
    %
    % %
    % % global Filename             % Filename prefix
    % %
    % %
    % % global SpecificationsValues;      % Contains the following variables:
    % % % General variables
    % % global Z_in;                % Number of blades
    % % global N_in;                % Propeller speed [RPM]
    % % global D_in;                % Propeller diameter [m]
    % % global T_in;                % Required thrust [N]
    % % global Vs_in;               % Ship speed [m/s]
    % % global Dhub_in;             % Hub diameter [m]
    % %
    % % global Js_in;               % ** Js = Vs/(n*D) ,  advance coefficient
    % % global KT_in;               % ** KT = THRUST/(rho*n^2*D^4)
    % % global n_in;                % ** propeller speed [rev/s] = Vs/(Js*D) = N/60
    % % % Advanced variables
    % % global rho_in;              % Water density
    % % global Mp_in;               % Number of vortex panels over the radius
    % % global Np_in;               % Number of points over the chord
    % % global ITER_in;             % Max iterations
    % % global Rhv_in;              % Hub vortex radius/hub radius
    % %     global SpecificationsText;
    % %
    % %
    % % global DuctValues;          % Contains the following variables:
    % % % Ducted Propeller variables
    % % global TAU;                 % Thrust ratio
    % % global CDd;                 % Duct section drag coefficient
    % %
    % %
    % % global CavValues;           % Contains the following variables:
    % % % Cavitation Analysis variables
    % % global H_in;                % Shaft centerline depth [m]
    % % global dV_in;               % Inflow variation [m/s]
    % %
    % %
    % % global FlagValues;          % Contains the following variables:
    % % % Flags
    % % global Propeller_flag;      % 0 == turbine, 1 == propeller
    % % global Hub_flag;            % 0 == no hub, 1 == hub
    % % global Duct_flag;           % 0 == no duct, 1 == duct
    % %
    % % global Chord_flag;          % ** CHORD OPTIMIZATION FLAG **
    % %
    % % global Viscous_flag;        % 0 == viscous forces off (CD = 0), 1 == viscous forces on
    % % global Plot_flag;           % 0 == do not display plots, 1 == display plots
    % %
    % % global FoilValues;          % Contains the following variables:
    % % % Airfoil type
    % % global Meanline;            % Meanline form
    % % global Thickness;           % Thickness form
    % %
    % % % Variables for geometry input table
    % % global XR_in;               % Radius ratio (r/R) at control point
    % % global XCoD_in;             % Camber to diameter ratio at XR
    % % global XCD_in;              % Coefficient of drag at XR
    % % global VAI_in;              % Va/Vs (wake profile at XR)
    % % global VTI_in;              % Vt/Vs (wake profile at XR)
    % % global Xt0oD_in;            % Thickness to camber ratio at XR
    % % global skew0_in;            % Skew at XR
    % % global rake0_in;            % Rake at XR

    global N_R0;                % Number of input radii

    global Col_Label;

    global XCoD_values XCLmax_values XCD_values;


    % --- Set GUI element variables as global ---

    % =========================================================================



    % =========================================================================
    % ========================== Initiate Main Figure =========================
     Meanline_cell   = {'NACA a=0.8' 'NACA a=0.8 (modified)' 'Parabolic'};
    Thickness_cell   = {'NACA 65A010' 'NACA 65A010 (modified)' 'Elliptical' 'Parabolic' 'NACA 66 (DTRC modified)'};

    % -------------------- Declare Default Input Variables --------------------

    % Variables for geometry input table
    XR_def          = [.2 .3 .4 .5 .6 .7 .8 .9 .95 1];          % Radius ratio (r/R)
    N_R0            = length(XR_def);                           % Number of input radii



    XCoD_def        = [0.1600 0.1812 0.2024 0.2196 0.2305 0.2311 0.2173 0.1807 0.1388 0.0010];     % chord to diameter ratio at XR
    XCD_def         = ones(1,N_R0).*0.008;                      % Coefficient of drag at XR

    % VAI_def         = ones(1,N_R0);                             % Va/Vs (wake profile at XR)
    % VTI_def         = zeros(1,N_R0);                            % Vt/Vs (wake profile at XR)

    % t0oc0_def       = [0.2055 0.1553 0.1180 0.09016 0.06960 0.05418 0.04206 0.03321 0.03228 0.03160];     % Thickness to chord ratio at XR

    % Xt0oD_def       = t0oc0_def .* XCoD_def;  % thickness / diameter at XR

      Xt0oD_def       = [0.0329    0.0281    0.0239    0.0198    0.0160    0.0125 0.0091    0.0060    0.0045  0];  % thickness / diameter at XR

    skew0_def       = zeros(N_R0);                              % Skew at XR
    rake0_def       = zeros(N_R0);                              % Rake at XR

    XCLmax_def      = 0.5 + (0.2-0.5)*(XR_def-XR_def(1))/(XR_def(end)-XR_def(1));  % CLmax distribution

    % --- Set defaults for change callbacks ---
      XCoD_values   =   XCoD_def;
    XCLmax_values   = XCLmax_def;
       XCD_values   =    XCD_def;
    % -----------------------------------------


    % General variables
    Z_def           = 3;                                        % Number of blades
    N_def           = 200;                                      % Propeller speed [RPM]
    D_def           = 2;                                        % Propeller diameter [m]
    T_def           = 25000;                                    % Required thrust [N]
    Vs_def          = 5;                                        % Ship speed [m/s]
    Dhub_def        = 0.4;                                      % Hub diameter [m]


    % Advanced variables
    rho_def         = 1000;                                     % Water density
    Mp_def          = 20;                                       % Number of vortex panels over the radius
    Np_def          = 20;                                       % Number of points over the chord


    n_def           = N_def/60;                                 % ** propeller speed [rev/s] = Vs/(Js*D) = N/60
    lambda_def      = n_def*2*pi*(D_def/2)/Vs_def;
    Js_def          = Vs_def/(n_def*D_def);                     % ** Js = Vs/(n*D) ,  advance coefficient
    KT_def          = T_def/(rho_def*n_def^2*D_def^4);          % ** KT = THRUST/(rho*n^2*D^4)

    CT_def          = T_def/(0.5*rho_def*Vs_def^2*pi*(D_def/2)^2);


    % Ducted Propeller variables
    TAU_def         = 1;                                        % Thrust ratio
    CDd_def         = 0.008;                                    % Duct section drag coefficient


    % Cavitation Analysis variables
    % H_def           = 3;                                        % Shaft centerline depth [m]
    % dV_def          = 0.3;                                      % Inflow variation [m/s]


    % % Flags
    % Propeller_flag  = 1;                                        % 0 == turbine, 1 == propeller
    % Hub_flag        = 1;                                        % 0 == no hub, 1 == hub
    % Duct_flag       = 0;                                        % 0 == no duct, 1 == duct
    %
    % Chord_flag      = 0;                                        % ** CHORD OPTIMIZATION FLAG **
    %
    % Viscous_flag    = 1;                                        % 0 == viscous forces off (CD = 0), 1 == viscous forces on
    % Plot_flag       = 0;                                        % 0 == do not display plots, 1 == display plots
    %

    % % Airfoil type
    % Thickness       = 'NACA66 (DTRC Modified)';                 % Thickness form
    % Meanline        = 'NACA a=0.8';                             % Meanline form

    filename        = 'DefaultPropeller';                           % Filename prefix

    % =========================================================================


    % =========================================================================
    % -------------------------- GUI Switching Check --------------------------
    if exist('OpenPropTempFile0307122010.mat','file')

        load('OpenPropTempFile0307122010.mat');

        Z_def       = ceil(mean([Zmin ZMax]));
        N_def       = mean([Nmin NMax]);
        D_def       = mean([Dmin DMax]);
        T_def       = THRUST;
        Vs_def      = Vs;
        Dhub_def    = Dhub;

        rho_def     = rho;
        Mp_def      = Mp;
        Np_def      = Np;

        XR_def      = XR;
        XCoD_def    = XCoD;
        XCD_def     = XCD;

        XCLmax_def  = 0.5 + (0.2-0.5)*(XR_def-XR_def(1))/(XR_def(end)-XR_def(1));  % CLmax distribution

    %     'VAI',
    %     'VTI',
    %     'ri'


        n_def           = N_def/60;                                 % ** propeller speed [rev/s] = Vs/(Js*D) = N/60
        lambda_def      = n_def*2*pi*(D_def/2)/Vs_def;
        Js_def          = Vs_def/(n_def*D_def);                     % ** Js = Vs/(n*D) ,  advance coefficient
        KT_def          = T_def/(rho_def*n_def^2*D_def^4);          % ** KT = THRUST/(rho*n^2*D^4)

        CT_def          = T_def/(0.5*rho_def*Vs_def^2*pi*(D_def/2)^2);


        % Ducted Propeller variables
        TAU_def         = 1;                                        % Thrust ratio
        CDd_def         = 0.008;                                    % Duct section drag coefficient

        delete('OpenPropTempFile0307122010.mat');

    end
    % -------------------------------------------------------------------------
    % =========================================================================


    % =========================================================================
    % -------------------------- GUI Layout Constants -------------------------
    %
    % This section presents the constant values for the dimensions of those
    % elements used in the GUI, as well as the construction formulas for the
    % different panels, including margins. The formulas are presented in a
    % linear reading order based on the actual order in the GUI.

    titlefontsize   = 25;
    subttlfontsize  = 20;
    panelfontsize   = 11;
    buttonfontsize  = 12;

    titleht         = 3;
    editbox         = 10;
    textbox         = 20;
    cavtextbox      = 22;
    filenametext    = 15;
    pushbox         = 10;
    runbox          = 20;

    selectboxht     = 2;
    selectbox       = 25;

    textht          = 1.5;
    editboxht       = 2;
    pushht          = 3;

    Specificationsboxht	= 1 + editboxht * 11 + 2;
    Specificationsbox 	= 2 + textbox + editbox + 2;

    Ductboxht           = 1 + editboxht * 3 + 2;
    Ductbox             = 2 + textbox + editbox + 2;

    BladeDesignboxht	= 1 + editboxht * 11 + 2;
    BladeDesignbox      = 2 + editbox * 6 + 2;

    Inflowboxht    	= 1 + editboxht * 11 + 2;
    Inflowbox      	= 2 + editbox * 3 + 2;

    % Cavboxht        = 1 + editboxht * 3 + 2;
    % Cavbox          = 2 + cavtextbox + editbox + 2;

    Valuesboxht        = 1 + editboxht * 3 + 2;
    Valuesbox          = 2 + cavtextbox + editbox + 2 + cavtextbox + 4 + editbox + 2;

    Flagboxht       = 2 + textht * 8 + 2;
    Flagbox         = 1 + textbox + 1;

    Foilboxht       = BladeDesignboxht - Flagboxht;
    Foilbox         = 1 + textbox + 1;

    Toolboxht       = 1 + pushht + 1 + editboxht + 2;
    Toolbox         = Specificationsbox + BladeDesignbox + Inflowbox + Flagbox - Ductbox - Valuesbox;

    filenamebox     = Toolbox - (2 + filenametext + 2);

    buttonspace     = (Toolbox-pushbox*2-runbox)/4;

    Windowht        = 1 + Ductboxht + Specificationsboxht + 1 + titleht + 1;
    Window          = 1 + Specificationsbox + BladeDesignbox + Inflowbox + Flagbox + 1;

    GUIselectionboxht   = 1 + selectboxht + 1;
    GUIselectionbox     = 1 + selectbox   + 1;
    % =========================================================================


    % =========================================================================
    % ------------------------- Create figure for GUI -------------------------
    close all;

    Fig_Main    = figure('units','characters','position',[5 55-Windowht Window Windowht],...
                         'numbertitle','off','name','OpenProp','menubar','none',...'toolbar','figure',...
                         'resize','off','color',[0.702 0.702 0.702]);

    % % -------------------------------
    % if strcmp(computer,'GLNX32') || strcmp(computer,'GLNXA64')

        set(Fig_Main,'resize','on');

    % end
    % % -------------------------------


    % -------------------------------------------------------------------------
    OpenPropDirectory = 'OpenProp_v3.3.4';
    OpenPropVersion   = 'OpenProp v3.3.4';

    Title       = uicontrol(Fig_Main,'style','text','fontsize',titlefontsize,...
                            'fontweight','bold','units','characters','position',...
                            [GUIselectionbox-12 Windowht-1-titleht Window-GUIselectionbox 3],'string',{OpenPropVersion});


    % -------------------------------------------------------------------------
    % --------- Setup panels --------

    % % GUIselection    = uibuttongroup('parent',Fig_Main,'fontsize',panelfontsize,...      %,'title',''
    % %                       'fontweight','bold','units','characters','position',...
    % %                       [1 1+Ductboxht+Specificationsboxht GUIselectionbox GUIselectionboxht],'clipping','on');

    Specifications	= uipanel('parent',Fig_Main,'title','Specifications','fontsize',panelfontsize,...
                          'fontweight','bold','units','characters','position',...
                          [1 1+Ductboxht Specificationsbox Specificationsboxht],'clipping','on');

    Duct            = uipanel('parent',Fig_Main,'title','    Ducted Propeller','fontsize',...
                          panelfontsize,'fontweight','bold','units','characters',...
                          'position',[1 1 Ductbox Ductboxht],'clipping','on');

    BladeDesign     = uipanel('parent',Fig_Main,'title','Blade Design Values','fontsize',...
                          panelfontsize,'fontweight','bold','units','characters',...
                          'position',[1+Ductbox 1+Ductboxht BladeDesignbox BladeDesignboxht],...
                          'clipping','on');

    Inflow          = uipanel('parent',Fig_Main,'title','Inflow Profile Values','fontsize',...
                          panelfontsize,'fontweight','bold','units','characters',...
                          'position',[1+Ductbox+BladeDesignbox 1+Ductboxht Inflowbox Inflowboxht],...
                          'clipping','on');

    % Cav             = uipanel('parent',Fig_Main,'title','    Cavitation Analysis','fontsize',...
    %                       panelfontsize,'fontweight','bold','units','characters',...
    %                       'position',[1+Ductbox 1 Cavbox Cavboxht],'clipping','on');

    Calculator      = uipanel('parent',Fig_Main,'title','Non-dimensional Parameters','fontsize',...
                          panelfontsize,'fontweight','bold','units','characters',...
                          'position',[1+Ductbox 1 Valuesbox Valuesboxht],'clipping','on');

    Flags           = uibuttongroup('parent',Fig_Main,'title','Options','fontsize',...
                          panelfontsize,'fontweight','bold','units','characters',...
                          'position',[1+Specificationsbox+BladeDesignbox+Inflowbox...
                          1+Ductboxht+Foilboxht Flagbox Flagboxht],'clipping','on',...
                          'selectionchangedfcn',@checkTurbine);

    Foil            = uipanel('parent',Fig_Main,'title','Airfoil type','fontsize',...
                          panelfontsize,'fontweight','bold','units','characters',...
                          'position',[1+Specificationsbox+BladeDesignbox+Inflowbox...
                          1+Ductboxht Foilbox Foilboxht],'clipping','on');

    Tools           = uipanel('parent',Fig_Main,'title','Tools','fontsize',...
                          panelfontsize,'fontweight','bold','units','characters',...
                          'position',[1+Ductbox+Valuesbox 1 Toolbox Toolboxht],...
                          'clipping','on');
    % -------------------------------------------------------------------------


    % -------------------------------------------------------------------------
    % ----------------------- Package Selection Elements ----------------------

    % % Select(1)     = uicontrol(Fig_Main,'units','characters','style',...
    % %                             'pushbutton','string','PS','position',...
    % %                             [1+1 1+Ductboxht+Specificationsboxht+2 selectbox selectboxht],...
    % %                             'horizontalalignment','left','callback','OpenPropParam',...
    % %                             'tooltipstring','Parametric Study');
    % %
    % % Select(2)     = uicontrol(Fig_Main,'units','characters','style',...
    % %                             'pushbutton','string','SP','position',...
    % %                             [1+1+selectbox 1+Ductboxht+Specificationsboxht+2 selectbox selectboxht],...
    % %                             'horizontalalignment','left','value',1,...
    % %                             'tooltipstring','Single Propeller Design');
    % %
    % % Select(3)     = uicontrol(Fig_Main,'units','characters','style',...
    % %                             'pushbutton','string','A','position',...
    % %                             [1+1+2*selectbox 1+Ductboxht+Specificationsboxht+2 selectbox selectboxht],...
    % %                             'horizontalalignment','left','callback','OpenPropAnalyze',...
    % %                             'tooltipstring','Analyze Propeller');
    % %

%     Selectcell      = {'Single Design','Parametric Study','Off-design Analysis'};
    Selectcell      = {'Single Design','Parametric Study'};

    Select          = uicontrol(Fig_Main,'units','characters','style','popupmenu',...
                                'position',[1+1 1+Ductboxht+Specificationsboxht+2 selectbox selectboxht],...
                                'backgroundcolor','w','string',Selectcell,'value',1,'callback',@Selectfn);
    % -------------------------------------------------------------------------


    % -------------------------------------------------------------------------
    % --------------------- Specifications Panel Elements ---------------------

    SpecificationsStrings   = {'Number of blades:'...
                               'Rotation speed (RPM):'...
                               'Rotor diameter (m):'...
                               'Required thrust (N):'...
                               'Ship speed (m/s):'...
                               'Hub diameter (m):'...
                               'Fluid density (kg/m^3):' ...
                               '# radial panels:'...
                               '# chordwise panels:'};

    SpecificationsValues_def      = [Z_def N_def D_def T_def Vs_def Dhub_def rho_def Mp_def Np_def];

    for index = 1 : length(SpecificationsValues_def)

        SpecificationsText(index)     = uicontrol(Specifications,'units','characters','style',...
                                            'text','string',SpecificationsStrings(index),...
                                            'position',[2 Specificationsboxht-4-editboxht*(index-1)...
                                            textbox textht],'horizontalalignment','left');

        SpecificationsValues(index)   = uicontrol(Specifications,'units','characters','style',...
                                            'edit','string',num2str(SpecificationsValues_def(index)),...
                                            'position',[2+textbox Specificationsboxht-4-editboxht*(index-1)...
                                            editbox editboxht],'backgroundcolor',[1 1 1]);
    end

    % Set callback for those Specs that affect Js and lambda

    set(SpecificationsValues(1),'callback',@checkBlades);
    set(SpecificationsValues(2),'callback',@updateValues);
    set(SpecificationsValues(3),'callback',@updateValues);
    set(SpecificationsValues(4),'callback',@updateValues);
    set(SpecificationsValues(5),'callback',@updateValues);
    set(SpecificationsValues(6),'callback',@updateValues);
    set(SpecificationsValues(7),'callback',@updateValues);
    % -------------------------------------------------------------------------


    % -------------------------- Duct Panel Elements --------------------------

    Ducttext(1)     = uicontrol(Duct,'units','characters','style',...
                                'text','string','Thrust Ratio:','position',...
                                [2 Ductboxht-4 textbox textht],...
                                'horizontalalignment','left');

    Ducttext(2)     = uicontrol(Duct,'units','characters','style',...
                                'text','string','Duct section drag (Cd):','position',...
                                [2 Ductboxht-4-editboxht textbox textht],...
                                'horizontalalignment','left');

    Ducttext(3)     = uicontrol(Duct,'units','characters','style',...
                                'text','string','duct D / prop D:','position',...
                                [2 Ductboxht-4-editboxht*2 textbox textht],...
                                'horizontalalignment','left');

    DuctValues(1)   = uicontrol(Duct,'units','characters','style',...
                                'edit','string',num2str(TAU_def),...
                                'position',[2+textbox Ductboxht-4,...
                                editbox editboxht],'backgroundcolor',[1 1 1],...
                               'enable','off');

    DuctValues(2)   = uicontrol(Duct,'units','characters','style',...
                                'edit','string',num2str(CDd_def),...
                                'position',[2+textbox Ductboxht-4-editboxht,...
                                editbox editboxht],'backgroundcolor',[1 1 1],...
                               'enable','off');

    DuctValuesoff(1)= uicontrol(Duct,'units','characters','style',...
                                'edit','string','1',...
                                'position',[2+textbox Ductboxht-4-editboxht*2,...
                                editbox editboxht],'backgroundcolor',[1 1 1],...
                                'enable','off');
    % -------------------------------------------------------------------------


    % ----------------------- Blade Design Panel Elements ---------------------

    ColName     = {'r/R' 'c/D' 'Cd' 't0/D' 'Skew' 'Xs/D'};

    for index = 1 : length(ColName)
        Col_Label(index) = uicontrol(BladeDesign,'style','edit','units','characters','FontSize',10,...
            'FontWeight','bold','position',[2+editbox*(index-1) BladeDesignboxht-4 editbox editboxht],...
            'string',ColName(index),'enable','inactive');
    end

    for index = 1 : N_R0

        XR_in(index)	= uicontrol(BladeDesign,'style','edit','units','characters','FontSize',10,...
                                   'backgroundcolor','w','string',num2str(XR_def(index)),'position',...
                                   [2 BladeDesignboxht-4-editboxht*index editbox editboxht]);

        XCoD_in(index)  = uicontrol(BladeDesign,'style','edit','units','characters','FontSize',10,...
                                    'backgroundcolor','w','string',num2str(XCoD_def(index)),'position',...
                                    [2+editbox BladeDesignboxht-4-editboxht*index editbox editboxht]);

        XCD_in(index)   = uicontrol(BladeDesign,'style','edit','units','characters','FontSize',10,...
                                    'backgroundcolor','w','string',num2str(XCD_def(index)),'position',...
                                    [2+editbox*2 BladeDesignboxht-4-editboxht*index editbox editboxht]);

    % 	f0oc_in(index)  = uicontrol(BladeDesign,'style','edit','units','characters','FontSize',10,...
    %                                 'backgroundcolor','w','string',num2str(f0oc_def(index)),'position',...
    %                                 [2+editbox*5 BladeDesignboxht-4-editboxht*index editbox editboxht]);

        Xt0oD_in(index)	= uicontrol(BladeDesign,'style','edit','units','characters','FontSize',10,...
                                     'backgroundcolor','w','string',num2str(Xt0oD_def(index)),'position',...
                                     [2+editbox*3 BladeDesignboxht-4-editboxht*index editbox editboxht]);

        skew0_in(index)	= uicontrol(BladeDesign,'style','edit','units','characters','FontSize',10,...
                                     'backgroundcolor','w','string',num2str(skew0_def(index)),'position',...
                                     [2+editbox*4 BladeDesignboxht-4-editboxht*index editbox editboxht]);

        rake0_in(index)	= uicontrol(BladeDesign,'style','edit','units','characters','FontSize',10,...
                                     'backgroundcolor','w','string',num2str(rake0_def(index)),'position',...
                                     [2+editbox*5 BladeDesignboxht-4-editboxht*index editbox editboxht]);

    end
    % -------------------------------------------------------------------------


    % ---------------------- Inflow Profile Panel Elements --------------------

    ColName2     = {'r' 'Va/Vs' 'Vt/Vs'};

    for index = 1 : length(ColName2)
        Col_Label2(index) = uicontrol(Inflow,'style','edit','units','characters','FontSize',10,...
            'FontWeight','bold','position',[2+editbox*(index-1) Inflowboxht-4 editbox editboxht],...
            'string',ColName2(index),'enable','inactive');
    end

    for index = 1 : N_R0

        ri_in(index)       = uicontrol(Inflow,'style','edit','units','characters','FontSize',10,...
                                   'backgroundcolor','w','string','','position',...
                                   [2 BladeDesignboxht-4-editboxht*index editbox editboxht]);

        VAI_in(index)   = uicontrol(Inflow,'style','edit','units','characters','FontSize',10,...
                                    'backgroundcolor','w','position',...
                                    [2+editbox BladeDesignboxht-4-editboxht*index editbox editboxht]);

        VTI_in(index)   = uicontrol(Inflow,'style','edit','units','characters','FontSize',10,...
                                    'backgroundcolor','w','position',...
                                    [2+editbox*2 BladeDesignboxht-4-editboxht*index editbox editboxht]);

    end

    % set(VAI_in(1),'string','1');
    % set(VTI_in(1),'string','0');
    % -------------------------------------------------------------------------


    % ------------------- Cavitation Analysis Panel Elements ------------------

    % Cavtext(1)     = uicontrol(Cav,'units','characters','style',...
    %                            'text','string','Shaft centerline depth (m):','position',...
    %                            [2 Cavboxht-4 cavtextbox textht],...
    %                            'horizontalalignment','left');
    %
    % Cavtext(2)     = uicontrol(Cav,'units','characters','style',...
    %                            'text','string','Inflow variation (m/s):','position',...
    %                            [2 Cavboxht-4-editboxht cavtextbox textht],...
    %                            'horizontalalignment','left');
    %
    % Cavtext(3)     = uicontrol(Cav,'units','characters','style',...
    %                            'text','string','Ideal angle of attackii:','position',...
    %                            [2 Cavboxht-4-editboxht*2 cavtextbox textht],...
    %                            'horizontalalignment','left');
    %
    % CavValues(1)   = uicontrol(Cav,'units','characters','style',...
    %                            'edit','string',num2str(TAU_def),...
    %                            'position',[2+cavtextbox Cavboxht-4,...
    %                            editbox editboxht],'backgroundcolor',[1 1 1],...
    %                            'enable','off');
    %
    % CavValues(2)   = uicontrol(Cav,'units','characters','style',...
    %                            'edit','string',num2str(CDd_def),...
    %                            'position',[2+cavtextbox Cavboxht-4-editboxht,...
    %                            editbox editboxht],'backgroundcolor',[1 1 1],...
    %                            'enable','off');
    %
    % CavValuesoff(1)= uicontrol(Cav,'units','characters','style',...
    %                            'edit','string','1',...
    %                            'position',[2+cavtextbox Cavboxht-4-editboxht*2,...
    %                            editbox editboxht],'backgroundcolor',[1 1 1],...
    %                            'enable','off');
    % -------------------------------------------------------------------------


    % -------------------------- Values Panel Elements ------------------------

    Valuestext(1)     = uicontrol(Calculator,'units','characters','style',...
                               'text','string','J = V/nD =','position',...
                               [2 Valuesboxht-4 cavtextbox textht],...
                               'horizontalalignment','left');

    Valuestext(2)     = uicontrol(Calculator,'units','characters','style',...
                               'text','string','L = omega*R/V =','position',...
                               [2 Valuesboxht-4-editboxht cavtextbox textht],...
                               'horizontalalignment','left');

    Valuestext(3)     = uicontrol(Calculator,'units','characters','style',...
                               'text','string','KT = T/(rho*n^2*D^4) =','position',...
                               [2+cavtextbox+editbox+2 Valuesboxht-4-editboxht cavtextbox+4 textht],...
                               'horizontalalignment','left');

    Valuestext(4)     = uicontrol(Calculator,'units','characters','style',...
                               'text','string','CT = T/(1/2*rho*V^2*pi*R^2) =','position',...
                               [2+cavtextbox+editbox+2 Valuesboxht-4 cavtextbox+4 textht],...
                               'horizontalalignment','left');

    Values(1)   = uicontrol(Calculator,'units','characters','style',...
                               'edit','value',Js_def,...
                               'position',[cavtextbox-4, Valuesboxht-4, editbox, editboxht],'backgroundcolor',[1 1 1],...
                               'enable','off');

    Values(2)   = uicontrol(Calculator,'units','characters','style',...
                               'edit','value',lambda_def,...
                               'position',[cavtextbox-4, Valuesboxht-4-editboxht, editbox, editboxht],'backgroundcolor',[1 1 1],...
                               'enable','off');

    Values(3)   = uicontrol(Calculator,'units','characters','style',...
                               'edit','value',KT_def,...
                               'position',[2+cavtextbox*2+editbox+2+4 Valuesboxht-4-editboxht,...
                               editbox editboxht],'backgroundcolor',[1 1 1],...
                               'enable','off');

    Values(4)   = uicontrol(Calculator,'units','characters','style',...
                               'edit','value',CT_def,...
                               'position',[2+cavtextbox*2+editbox+2+4 Valuesboxht-4,...
                               editbox editboxht],'backgroundcolor',[1 1 1],...
                               'enable','off');

    %Values(3)= uicontrol(Calculator,'units','characters','style',...
    %                            'edit','string','1',...
    %                            'position',[2+cavtextbox Valuesboxht-4-editboxht*2,...
    %                            editbox editboxht],'backgroundcolor',[1 1 1],...
    %                            'enable','off');

    % -------------------------------------------------------------------------


    % -------------------------------------------------------------------------
    % -------------------------- Flags Panel Elements --------------------------

    FlagValues(1)     = uicontrol(Flags,'units','characters','style','radiobutton',...
                                  'string','Propeller','value',1,'position',[1 Flagboxht-4 textbox textht]);

    FlagValues(2)     = uicontrol(Flags,'units','characters','style','radiobutton',...
                                  'string','Turbine','position',[1 Flagboxht-4-textht textbox textht]);

    FlagValues(3)     = uicontrol(Flags,'units','characters','style','checkbox',...
                                  'string','Hub','value',1,'position',[1 Flagboxht-4-textht*2 textbox textht]);

    % ---------------------------------------------------
    % ------ Ducted flag relocated to Ducted panel ------

    % FlagValues(4)     = uicontrol(Flags,'units','characters','style','checkbox',...
    %                               'string','Ducted','position',...
    %                               [1 Flagboxht-4-textht*3 textbox textht]);

    FlagValues(4)     = uicontrol(Duct,'units','characters','style','checkbox',...
                                  'string','','position',[1 Ductboxht-textht...
                                  textbox textht],'callback',@changeDuct);
    % ---------------------------------------------------

    FlagValues(5)     = uicontrol(Flags,'units','characters','style','checkbox',...
                                  'string','Chord optimization','position',[1 Flagboxht-4-textht*3 textbox textht],...
                                  'callback',@changeChord);

    FlagValues(6)     = uicontrol(Flags,'units','characters','style','checkbox',...
                                  'string','Viscous forces','value',1,'position',[1 Flagboxht-4-textht*4 textbox textht],...
                                  'callback',@changeViscous);

    FlagValues(7)     = uicontrol(Flags,'units','characters','style','checkbox',...
                                  'string','Optimization plots','position',[1 Flagboxht-4-textht*6 textbox textht]);

    FlagValues(8)     = uicontrol(Flags,'units','characters','style','checkbox',...
                                  'string','Geometry plots','value',1,'position',[1 Flagboxht-4-textht*7 textbox textht]);


    FlagValues(9)     = uicontrol(Flags,'units','characters','style','checkbox',...
                                  'string','Performance curve','value',0,'position',[1 Flagboxht-4-textht*8 textbox textht]);

    % ------ Added Cav_flag ------

    % FlagValues(10)     = uicontrol(Cav,'units','characters','style','checkbox',...
    %                               'string','','position',[1 Cavboxht-textht...
    %                               textbox textht],'callback',@changeCav);

    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------

    % ---------------------- Airfoil Type Panel Elements ----------------------
    FoilText(1)     = uicontrol(Foil,'units','characters','style','text',...
                                'string','Meanline type:','position',...
                                [1 Foilboxht-3.5 textbox textht]);

    FoilValues(1)	= uicontrol(Foil,'units','characters','style','popupmenu',...
                                'position',[1 Foilboxht-3.5-textht textbox textht],...
                                'backgroundcolor','w','string',...
                                Meanline_cell);

    FoilText(2)     = uicontrol(Foil,'units','characters','style','text',...
                                'string','Thickness type:','position',...
                                [1 Foilboxht-3.5-textht*2 textbox textht]);

    FoilValues(2)	= uicontrol(Foil,'units','characters','style','popupmenu',...
                                'position',[1 Foilboxht-3.5-textht*3 textbox textht],...
                                'backgroundcolor','w','string',...
                                Thickness_cell);
    % -------------------------------------------------------------------------


    % -------------------------- Tools Panel Elements -------------------------
    ToolText     = uicontrol(Tools,'units','characters','style','text',...
                                'string','Filename prefix:','horizontalalignment',...
                                'left','position',[2 Toolboxht-4 filenametext textht]);

    Filename     = uicontrol(Tools,'units','characters','style','edit',...
                                'backgroundcolor','w','string',filename,'position',...
                                [2+filenametext Toolboxht-4 filenamebox editboxht],...
                                'callback',@changeDir);

    LoadButton      = uicontrol(Tools,'units','characters','style','pushbutton',...
                                'string','Load','fontsize',buttonfontsize,...
                                'position',[buttonspace 1 pushbox pushht],...
                                'callback',@loaddata);

    SaveButton      = uicontrol(Tools,'units','characters','style','pushbutton',...
                                'string','Save','fontsize',buttonfontsize,...
                                'position',[buttonspace*2+pushbox 1 pushbox pushht],...
                                'callback',@savedata);

    RunButton       = uicontrol(Tools,'units','characters','style','pushbutton',...
                                'string','Run OpenProp','fontsize',buttonfontsize,...
                                'position',[buttonspace*3+pushbox*2 1 runbox pushht],...
                                'callback',@execute);
    % -------------------------------------------------------------------------

end

