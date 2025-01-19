%%
% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
% -------------------------------------------------------------------------
%               And the meek shall inherit the earth...                   %
% -------------------------------------------------------------------------
% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
%%
% ================================ execute ================================
function execute(hObject,ED)
    %
    % This subfunction runs OpenProp.
    %
    newPlots();

    global Plots PlotPanels Toggle OnDesignValues ConversionValues systemToggle;

    global OpenPropDirectory SpecificationsValues DuctValues FlagValues FoilValues Filename...
           XR_in XCoD_in XCD_in VAI_in ri_in VTI_in Xt0oD_in skew0_in rake0_in...
           Meanline_cell Thickness_cell; % CavValues

    global pt

    % -------------------------------------------------------------------------
    % Figure out what the current working directory is,
    % and change directories to OpenPropDirectory/filename
    filename = get(Filename,'string');                       % Filename prefix

    rest = pwd;

    while ~isempty(rest)
        [CurrentDirectory,rest] = strtok(rest,'/');

        if strcmp(CurrentDirectory,OpenPropDirectory)

            if isempty(rest)

                % you are in /OpenPropDirectory/
                mkdir(['./',filename])
                   cd(['./',filename])
                   addpath ../SourceCode

            elseif strcmp(rest(2:end),filename)
                % already in /OpenPropDirectory/filename/
                addpath ../SourceCode
                rest = [];

            elseif strcmp(rest(2:end),'SourceCode')

                mkdir(['../',filename])
                   cd(['../',filename])
                   addpath ../SourceCode
                rest = [];

            else
                % you are in /OpenPropDirectory/wrongfolder
                disp('ERROR2: Must start OpenProp from the root directory.')
                return
            end
        end
    end
    % -------------------------------------------------------------------------


    % --------------------------- Design parameters ---------------------------

    Z           = str2double(get(SpecificationsValues(1),'string'));  % number of blades
    N           = str2double(get(SpecificationsValues(2),'string'));  % propeller speed [RPM]
    D           = str2double(get(SpecificationsValues(3),'string'));	% propeller diameter [m]
    THRUST      = str2double(get(SpecificationsValues(4),'string')); 	% required thrust [N]
    Vs          = str2double(get(SpecificationsValues(5),'string'));  % ship velocity [m/s]
    Dhub        = str2double(get(SpecificationsValues(6),'string'));  % hub diameter [m]
    rho         = str2double(get(SpecificationsValues(7),'string')); 	% water density [kg/m^3]
    Mp          = str2double(get(SpecificationsValues(8),'string')); 	% number of vortex panels over the radius
    Np          = str2double(get(SpecificationsValues(9),'string')); 	% Number of points over the chord [ ]

    ITER        = 40;   % number of iterations in analysis
    Rhv         = 0.5;	% hub vortex radius / hub radius

    TAU         = str2double(get(DuctValues(1),'string'));      % Thrust ratio
    CDd         = str2double(get(DuctValues(2),'string'));      % Duct section drag coefficient

    % H           = str2double(get(CavValues(1),'string'));       % Shaft centerline depth [m]
    % dV          = str2double(get(CavValues(2),'string'));       % Inflow variation [m/s]


    % --------------------------------- Flags ---------------------------------

    Propeller_flag	= get(FlagValues(1),'value');               % 0 == turbine, 1 == propeller
    Hub_flag  	    = get(FlagValues(3),'value');                   % 0 == no hub, 1 == hub
    Duct_flag	    = get(FlagValues(4),'value');                   % 0 == no duct, 1 == duct

    Chord_flag	    = get(FlagValues(5),'value');                   % ** CHORD OPTIMIZATION FLAG **

    Viscous_flag	= get(FlagValues(6),'value');               % 0 == viscous forces off (CD = 0), 1 == viscous forces on
    Plot_flag       = get(FlagValues(7),'value');               % 0 == do not display plots, 1 == display plots



    Make2Dplot_flag = get(FlagValues(8),'value');               % 0 == do not make a 2D plot of the results, 1 == make plot
    Make3Dplot_flag = get(FlagValues(8),'value');               % 0 == do not make a 3D plot of the results, 1 == make plot

    Analyze_flag	= get(FlagValues(9),'value');

    % Cav_flag	= get(FlagValues(10),'value');                   % 0 == do not run cavitation mapping, 1 == run cavitation mapping


    % -------------------------------------------------------------------------
    % ---------------------- Blade 2D section properties ----------------------

    Meanline_index  = get(FoilValues(1),'value');
    Meanline        = char(Meanline_cell(Meanline_index));              % Meanline form

    Thickness_index	= get(FoilValues(2),'value');
    Thickness       = char(Thickness_cell(Thickness_index));            % Thickness form


    XR          = str2double(get(XR_in,  'string'));             	% radius / propeller radius

    XCoD        = str2double(get(XCoD_in,'string'));           	% chord / diameter
    XCLmax      = str2double(get(XCoD_in,'string'));            % maximum lift coefficient (for chord optimization)

    XCD     	= str2double(get(XCD_in, 'string'));            	% section drag coefficient

    ri          = str2double(get(ri_in, 'string'));
    VAI         = str2double(get(VAI_in,'string'));             % axial      inflow velocity / ship velocity
    VTI         = str2double(get(VTI_in,'string'));          	% tangential inflow velocity / ship velocity

    ri  = ri (~isnan(ri));
    VAI = VAI(~isnan(VAI));
    VTI = VTI(~isnan(VTI));


    % f0oc0     = str2double(get(f0oc_in,'String'));          	% max section camber    / chord
    Xt0oD       = str2double(get(Xt0oD_in,'String'));           % max section thickness / chord
    skew0       = str2double(get(skew0_in,'String'));          	% skew
    rake0       = str2double(get(rake0_in,'String'));         	% rake


    % ----------------------- Compute derived quantities ----------------------

    n           = N/60;                                         % ** propeller speed [rev/s] = Vs/(Js*D) = N/60
    Js          = Vs/(n*D);                                     % ** Js = Vs/(n*D) ,  advance coefficient
    KT          = THRUST/(rho*n^2*D^4);                        	% ** KT = THRUST/(rho*n^2*D^4)
    L           = pi/Js;                                        % tip speed ratio

    R           = D/2;                                          % propeller radius [m]
    Rhub        = Dhub/2;                                       % hub radius [m]
    Rhub_oR     = Rhub/R;

    CTDES       = THRUST/(0.5*rho*Vs^2*pi*R^2);                 % CT thrust coefficient required

    % dVs         = dV/Vs;                                        % axial inflow variation / Vs



    % *************************************************************************
    % *************************************************************************
    input.part1      = '------ Performance inputs ------';
    input.Z          = Z;           % [1 x 1], [ ] number of blades
    input.N          = N;           % propeller speed [RPM]
    input.D          = D;           % propeller diameter [m]
    input.Vs         = Vs;          % [1 x 1], [m/s] ship speed
    input.Js         = Js;          % [1 x 1], [ ] advance coefficient, Js = Vs/nD = pi/L
    input.L          = L;           % [1 x 1], [ ] tip speed ratio, L = omega*R/V
    input.THRUST     = THRUST;      % required thrust [N]
    input.CTDES      = CTDES;       % [1 x 1], [ ] desired thrust coefficient
    input.TAU        = TAU;          % Thrust ratio

    input.part2      = '------ Geometry inputs ------';
    input.Mp         = Mp;          % [1 x 1], [ ] number of blade sections
    input.Np         = Np;          % [1 x 1], [ ] number of points along the chord
    input.R          = R;           % [1 x 1], [m] propeller radius
    input.Rhub       = Rhub;        % [1 x 1], [m] hub radius
    input.XR         = XR;          % [length(XR) x 1], [ ] input radius/propeller radiusat XR
    input.XCD        = XCD;         % [length(XR) x 1], [ ] input drag coefficient       at XR
    input.XCoD       = XCoD;        % [length(XR) x 1], [ ] input chord / diameter       at XR
    input.Xt0oD      = Xt0oD;       % [length(XR) x 1], [ ] input thickness / chord      at XR
    input.skew0      = skew0;       % [length(XR) x 1], [ ] input skew  [deg]      at XR
    input.rake0      = rake0;       % [length(XR) x 1], [ ] input rake X/D       at XR
    input.Meanline   = Meanline;    % 2D section meanline  flag
    input.Thickness  = Thickness;   % 2D section thickness flag
    input.XCLmax     = XCLmax;

    if ~isempty(ri)  , input.ri  = ri;  end
    if ~isempty(VAI) , input.VAI = VAI; end        % [length(XR) x 1], [ ] input axial inflow velocity  at XR
    if ~isempty(VTI) , input.VTI = VTI; end        % [length(XR) x 1], [ ] input swirl inflow velocity


    input.Rduct     = R;
    input.Cduct     = D/2;
    input.CDd       = CDd;

    input.part3      = '------ Computational inputs ------';
    input.Propeller_flag  = Propeller_flag; % 0 == turbine, 1 == propeller
    input.Viscous_flag    = Viscous_flag;   % 0 == viscous forces off (CD = 0), 1 == viscous forces on
    input.Hub_flag        = Hub_flag;       % 0 == no hub, 1 == hub
    input.Duct_flag       = Duct_flag;      % 0 == no duct, 1 == duct
    input.Plot_flag       = Plot_flag;      % 0 == do not display plots, 1 == display plots
    input.Chord_flag      = Chord_flag;     % 0 == do not optimize chord lengths, 1 == optimize chord lengths

    input.Make2Dplot_flag = Make2Dplot_flag;
    input.Make3Dplot_flag = Make3Dplot_flag;
    % input.Make_Rhino_flag = Make_Rhino_flag;
    input.ITER            = ITER;           % [ ] number of iterations
    input.Rhv              = Rhv;         % [1 x 1], [ ] hub vortex radius / hub radius

    input.part4      = '------ Cavitation inputs ------';
    input.rho        = rho;         % [1 x 1], [kg/m^3] fluid density
    % input.dVs        = dVs;         % [1 x 1], [ ] ship speed variation / ship speed
    % input.H          = H;           % [1 x 1]

    input.part5      = '------ Duct inputs ------';


    % ---------------------------- Pack up propeller/turbine data structure, pt
    pt.filename = filename; % (string) propeller/turbine name
    pt.date     = date;     % (string) date created
  % pt.notes    = ' ';    % (string or cell matrix)   notes
    pt.input    = input;    % (struct) input parameters
    pt.design   = [];       % (struct) design conditions
    pt.geometry = [];       % (struct) design geometry
    pt.states	= [];       % (struct) off-design state analysis

    pt.input.GUI_flag = 1;

    % *************************************************************************
    % *************************************************************************



    % =========================================================================
    % ============================ execution script ===========================

    % ---------------------------------------------------------------------
    % Plot from input:

    % Expanded Blade, only if Chord_flag = 0
    if Chord_flag == 0

        set(0,'CurrentFigure',Plots);
        h = axes('parent',PlotPanels(1));
        axes(h);
        hold on;

        XXR   = XR(1) + (XR(end)-XR(1))*(sin((0:60)*pi/(2*60)));
        XXCoD = InterpolateChord(XR,XCoD,XXR);

        plot(XXR, XXCoD,'b','LineWidth',2);
        plot(XXR,-XXCoD,'b','LineWidth',2);

        plot(XR, XCoD,'.b','MarkerSize',16)
        plot(XR,-XCoD,'.b','MarkerSize',16)

        xlabel('r/R','Fontsize',16,'FontName','Times');
        ylabel('c/R','Fontsize',16,'FontName','Times');
        set(gca,'Fontsize',14,'FontName','Times')
        grid on, box on,

    else
        set(Toggle(1),'enable','off');
    end

    % Thickness profile:
        set(0,'CurrentFigure',Plots);
        h = axes('parent',PlotPanels(2));
        axes(h);
        hold on;

        XXR    = XR(1) + (XR(end)-XR(1))*(sin((0:60)*pi/(2*60)));
        XXt0oD = pchip(XR,Xt0oD,XXR);

        plot(XXR, XXt0oD,'b','LineWidth',2);
        plot(XXR,-XXt0oD,'b','LineWidth',2);

        plot(XR, Xt0oD,'.b','MarkerSize',16)
        plot(XR,-Xt0oD,'.b','MarkerSize',16)

        xlabel('r/R','FontSize',16,'FontName','Times');
        ylabel('t0/D','FontSize',16,'FontName','Times');
        set(gca,     'FontSize',14,'FontName','Times');
        grid on, box on,

    % Inflow profile:
        set(0,'CurrentFigure',Plots);
        h = axes('parent',PlotPanels(3));
        axes(h);
        hold on;

        if ~isempty(VAI)

            plot(VAI,ri,'b','LineWidth',2);
            plot(VTI,ri,'r','LineWidth',2);

            xlim([-0.1 1.1*max(VAI)])
        else
            plot( ones(size(XR)),XR,'b','LineWidth',2);
            plot(zeros(size(XR)),XR,'r','LineWidth',2);

            xlim([-0.1 1.1])
        end

        xlabel('VA / Vs (blue),   VT / Vs (red)','FontSize',16,'FontName','Times');
        ylabel('r / R','FontSize',16,'FontName','Times');
        set(gca,     'FontSize',14,'FontName','Times');
        grid on, box on,


    % ---------------------------------------------------------------------
    % Perform design optimization
    pt.design   = EppsOptimizer(input);
    % ---------------------------------------------------------------------


    % ---------------------------------------------------------------------
    % Set On Design Performance values

    pt.design.Q = pt.design.CQ * 0.5 * rho * Vs^2 * pi*D^2/4 * D/2; % [Nm]  torque

    omega = 2*pi*n; % [rad/s]

    pt.design.P = pt.design.Q * omega;

    if Propeller_flag == 1
        set(OnDesignValues(1),'string',num2str(pt.design.Js));
        set(OnDesignValues(2),'string',num2str(pt.design.KT));
        set(OnDesignValues(3),'string',num2str(pt.design.KQ));
        set(OnDesignValues(4),'string',num2str(pt.design.EFFY));
    else
        set(OnDesignValues(1),'string',num2str(pi/pt.design.L));
        set(OnDesignValues(2),'string',' ');
        set(OnDesignValues(3),'string',' ');
        set(OnDesignValues(4),'string',' ');
    end
%     set(OnDesignValues(2),'string',num2str(pt.design.KT));
%     set(OnDesignValues(3),'string',num2str(pt.design.KQ));
%     set(OnDesignValues(4),'string',num2str(pt.design.EFFY));

    if Propeller_flag == 1
        set(OnDesignValues(5),'string',num2str(pt.design.ADEFFY));
    end

    set(OnDesignValues(6),'string',num2str(pt.design.CT));
    set(OnDesignValues(7),'string',num2str(pt.design.CQ));
    set(OnDesignValues(8),'string',num2str(pt.design.CP));

    set(ConversionValues(1),'string',num2str(pt.input.Vs));
    set(ConversionValues(2),'string',num2str(pt.input.N));
    set(ConversionValues(3),'string',num2str(pt.input.D));
    set(ConversionValues(4),'string',num2str(pt.input.THRUST));
    set(ConversionValues(5),'string',num2str(pt.design.Q));
    set(ConversionValues(6),'string',num2str(pt.design.P));

    set(systemToggle,'enable','on');
    % ---------------------------------------------------------------------



    % ---------------------------------------------------------------------
    disp(' ')
    disp('Creating graphical and text reports')
    disp(' ')

    % Create graphical and text reports
    Make_Reports(pt);
    % ---------------------------------------------------------------------


    % ---------------------------------------------------------------------
    % Plot t0/D vs RC
        set(0,'CurrentFigure',Plots);
        h = axes('parent',PlotPanels(11));
        axes(h);
        hold on;

        plot([Rhub,pt.design.RC,1],interp1(pt.design.RC, pt.design.t0oD,[Rhub,pt.design.RC,1],'spline','extrap'),'b','LineWidth',2);
        plot([Rhub,pt.design.RC,1],interp1(pt.design.RC,-pt.design.t0oD,[Rhub,pt.design.RC,1],'spline','extrap'),'b','LineWidth',2);

        plot(pt.design.RC, pt.design.t0oD,'b.','MarkerSize',16);
        plot(pt.design.RC,-pt.design.t0oD,'b.','MarkerSize',16);

        xlabel('r/R','FontSize',16,'FontName','Times');
        ylabel('t0/D','FontSize',16,'FontName','Times');
        set(gca,     'FontSize',14,'FontName','Times');
        grid on; box on,


    % ---------------------------------------------------------------------
    % Plot CL vs RC
        set(0,'CurrentFigure',Plots);
        h = axes('parent',PlotPanels(12));
        axes(h);
        hold on;

        plot([Rhub,pt.design.RC,1],interp1(pt.design.RC, pt.design.CL,[Rhub,pt.design.RC,1],'spline','extrap'),'b','LineWidth',2);

        plot(pt.design.RC,pt.design.CL,'b.','MarkerSize',16);

        xlabel('r/R','FontSize',16,'FontName','Times');
        ylabel('CL','FontSize',16,'FontName','Times');
        set(gca,     'FontSize',14,'FontName','Times');
        grid on; box on,


            % Set lower y limit to 0
            ylimits = get(gca,'Ylim');

            if abs(ylimits(2)) > abs(ylimits(1))

                set(gca,'Ylim',[0 ylimits(2)]);
            else
                set(gca,'Ylim',[ylimits(1) 0]);
            end
    % ---------------------------------------------------------------------


    % ---------------------------------------------------------------------
    % Determine propeller geometry
    if Make2Dplot_flag == 1 | Make3Dplot_flag == 1
        pt.geometry = Geometry(pt);
    end
    % ---------------------------------------------------------------------

    % ---------------------------------------------------------------------
    % if Cav_flag
    %
    %     % Pefrorm cavitation analysis
    %     Cav_CavitationMap(pt);
    %
    %     VLMbucket
    %
    % end
    % ---------------------------------------------------------------------


    % ---------------------------------------------------------------------
    if Analyze_flag == 1
        % % Analyze off-design states
        % Js_all      = [1.05:-0.1:0.55];     % advance coefficient
        % LAMBDAall   = pi./Js_all;           % tip-speed ratio
        % pt.states	= Analyze(pt,LAMBDAall)

        pt.states      = AnalyzeAuto(pt);

        % ---------------------------------------------------------------------
        set(0,'CurrentFigure',Plots);
        h = axes('parent',PlotPanels(15));
        axes(h);
        hold on;

        if Propeller_flag == 1

            VMIV = pt.design.VMIV;

            Js_curve = linspace(min(pt.states.Js),max(pt.states.Js),100);

            EFFY_curve = (Js_curve/(2*pi)) * VMIV .* pchip(pt.states.Js,pt.states.KT,Js_curve)./pchip(pt.states.Js,pt.states.KQ,Js_curve);

            % Efficiency (green squares)
                    plot(Js_curve,EFFY_curve,'-','LineWidth',2,'Color',[0 0.8 0])
            Heffy = plot(pt.states.Js,pt.states.EFFY,'sk','MarkerSize',5,'LineWidth',1,'MarkerFaceColor',[0 0.8 0]);

            % Thrust coefficient (blue diamonds)
                    plot(pt.states.Js,pt.states.KT,'b-','LineWidth',2)
            Hkt   = plot(pt.states.Js,pt.states.KT,'dk','MarkerSize',5,'LineWidth',1,'MarkerFaceColor','b');

            % Torque coefficient (red circles)
                  plot(pt.states.Js,10*pt.states.KQ,'r-','LineWidth',2)
            Hkq = plot(pt.states.Js,10*pt.states.KQ,'ok','MarkerSize',5,'LineWidth',1,'MarkerFaceColor','r');

            % Design point
            plot(pt.design.Js*[1 1],[0 2],'k--','LineWidth',1);

            xlabel('Js','FontSize',16,'FontName','Times'),
            ylabel('KT, 10*KQ, EFFY','FontSize',16,'FontName','Times')
            axis([min(pt.states.Js) max(pt.states.Js) 0 0.9])
            set(gca,     'FontSize',14,'FontName','Times');
            box on, grid on,


            ylimits = get(gca,'Ylim');
            set(gca,'Ylim', [0  max(1,ylimits(2)) ] );

        else

            % Power coefficient (blue dots)
            plot(pt.states.L,-pt.states.CP,'b.-','LineWidth',2,'MarkerSize',12)

            % Design point
            plot(pt.design.L*[1 1],[0 0.6],'k--','LineWidth',1);

            % % Betz limit
            % plot([0 ceil(max(pt.states.L))],(16/27)*[1 1],'k--','LineWidth',1);

            xlabel('L','FontSize',16,'FontName','Times'),
            ylabel('CP','FontSize',16,'FontName','Times'),

            set(gca,'Ytick',[0:0.1:0.6])

            axis([0 ceil(max(pt.states.L))  0 0.6])
            set(gca,     'FontSize',14,'FontName','Times');
            box on, grid on,

        end

    end
    % ---------------------------------------------------------------------



    % ---------------------------------------------------------------------
    % Plot parametric study results
    if exist([filename '.mat'],'file')

                CLR = [     1       0       0;      ... % (1) Red
                            0       0.9     0;      ... % (2) Green
                            0       0       1;      ... % (3) Blue
                            0.75    0       0.75;   ... % (4) Purple
                            1       0.5     0;      ... % (5) Orange
                            0       1       1;      ... % (6) Cyan
                            1       0       1;      ... % (7) Magenta
                            0.75    0.5     0.25;   ... % (8) Brown
                            0.25    0.25    0.75;   ... % (9) Navy blue
                            0.25    0.5     0.75;   ... % (10) Steel blue
                            0.75    0.75    0];         % (11) Burnt Yellow

        temp    = pt;

        load([filename '.mat']);

        if isfield(pt,'paroutput')

            paroutput   = pt.paroutput;
            parinput    = pt.parinput;

            % --- For EFFY vs N ---

            set(0,'CurrentFigure',Plots);
            h = axes('parent',PlotPanels(4));
            axes(h);

            hold on, box on, grid on,
            ylim([0 1]),

            tempEFFY = paroutput.EFFY(1,:,1);

            plot(paroutput.N,tempEFFY(:),'-','color',CLR( (mod(0,11)+1) ,:));

            str_legend ={[num2str(parinput.D),' m  ']};


            legend(str_legend,'location','southwest');

            xlabel('Rotation Speed (RPM)  ');
            ylabel('Efficiency');
            title(['Number of Blades: ',num2str(Z),'  '])

            % --- For EFFY vs D ---

            set(0,'CurrentFigure',Plots);
            h = axes('parent',PlotPanels(5));
            axes(h);

            hold on, box on, grid on,
            ylim([0 1]),

            tempEFFY = paroutput.EFFY(1,1,:);

            plot(paroutput.D,tempEFFY(:),'-','color',CLR( (mod(0,11)+1) ,:));

            str_legendD ={[num2str(parinput.N),' RPM  ']};

            legend(str_legendD,'location','southwest');

            xlabel('Propeller Diameter (m)  ');
            ylabel('Efficiency');
            title(['Number of Blades: ',num2str(Z),'  '])

            % =========
            % =========

        else

            set(Toggle(4),'enable','off');
            set(Toggle(5),'enable','off');

        end

        pt      = temp;

    else

        set(Toggle(4),'enable','off');
        set(Toggle(5),'enable','off');

    end
    % ---------------------------------------------------------------------


    % ---------------------------------------------------------------------
    % pt overwrite sequence:

    if exist([filename '.mat'],'file')

        disp(['Found original file:  ',filename,'.mat']);

        temp  = pt;

        load(filename)

        disp(['Overwriting file:  ',filename,'.mat']);

        pt.filename = temp.filename;
        pt.date     = temp.date;
        pt.input	= temp.input;
        pt.design   = temp.design;
        pt.geometry = temp.geometry;
        pt.states   = temp.states;

    end

    pt

    save(filename,'pt');

    save([filename,'_GUIsd'],'Z','N','D','THRUST','Vs','Dhub','rho','Mp','Np',...
        'TAU','CDd','Propeller_flag','Hub_flag','Duct_flag',...
        'Chord_flag','Viscous_flag','Plot_flag','Make2Dplot_flag','Make3Dplot_flag','Analyze_flag',...
        'Meanline','Meanline_index','Thickness',...
        'Thickness_index','filename','XR','XCoD','XCD','ri','VAI','VTI','Xt0oD','skew0',...
        'rake0');
    % ---------------------------------------------------------------------

end

