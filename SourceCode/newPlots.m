% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function newPlots

    global PlotPanels Toggle OnDesignValues systemToggle ConversionValues UnitsText;

    % ---------------------------------------------------------------------
    % -------------------------- GUI Layout Constants -------------------------
    %
    % This section presents the constant values for the dimensions of those
    % elements used in the GUI, as well as the construction formulas for the
    % different panels, including margins. The formulas are presented in a
    % linear reading order based on the actual order in the GUI.

    titlefontsize   = 25;
    panelfontsize   = 11;

    titleht         = 3;
    textht          = 1.5;
    editboxht       = 2;
    pushht          = 2; % 1.75;

    textbox         = 20;
    ondesignbox     = 15;
    editbox         = 10;

    Windowht        = 44; % 40;           % 1 + Ductboxht + Specificationsboxht + 1 + titleht + 1;
    Window          = 142;          % 1 + Specificationsbox + Inputbox + Flagbox + 1;

    Togglebox       = textbox + 2;
    Toggleboxht     = Windowht - 1 - titleht - 1;

    Displaybox      = Window - 1 - Togglebox - 1;
    Displayboxht    = Windowht - 1 - titleht - 1;

    Conversionbox   = 2 + textbox + ondesignbox + editbox + 2;
    Conversionboxht	= pushht * 8 * 1.25 + 1;

    % =========================================================================
    % ---------------------- Create figure for Plots GUI ----------------------

    Plots       = figure('units','characters','position',[25 20 Window Windowht],...
                         'name','Plots','numbertitle','off','toolbar','figure',... 'menubar','none'
                         'color',[0.702 0.702 0.702],'resize','off');
    % --- --- ---

    Title       = uicontrol(Plots,'style','text','fontsize',titlefontsize,...
                            'fontweight','bold','units','characters','position',...
                            [0 Windowht-1-titleht Window 3],'string',{'OpenProp v3.3.4'});

    % =========================================================================
    % --------- Setup panels --------

    Togglebar           = uibuttongroup('parent',Plots,'title','Figures and plots',...
                                        'fontsize',panelfontsize,'fontweight','bold',...
                                        'units','characters','position',[1 1 Togglebox...
                                        Toggleboxht],'clipping','on','selectionchangedfcn',...
                                        @togglefn);

    % ---
    %
    % The panels defined from here on are set up using a (chronological) order
    % criteria based on the different design stages, as presented in the following
    % list:
    %
    % Inputs:
    % 1 Expanded Blade
    % 2 Thickness profile
    % 3 Inflow profile
    %
    % Parametric:
    % 4 Efficiency vs Diameter
    % 5 Efficiency vs Rotational Speed (N)
    %
    %
    % Design results:
    %
    % 6 On Design Perf. (default)
    %
    % 7 Circulation distribution
    % 8 Induced Velocity
    % 9 Beta & Beta i
    %
    % 10 Expanded Blade (after)
    % 11 t0/D vs Radius at control points
    % 12 CL vs Radius at control points
    %
    % 13 2D Geometry
    % 14 3D Geometry
    %
    % 15 Performance curves
    %
    % ---


    PlotPanels(1)       = uipanel('parent',Plots,'title','Expanded Blade (input blade)',  	'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(2)       = uipanel('parent',Plots,'title','Blade Thickness (input blade)', 	'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(3)       = uipanel('parent',Plots,'title','Inflow Profile',                	'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(4)       = uipanel('parent',Plots,'title','Efficiency vs Diameter',          'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(5)       = uipanel('parent',Plots,'title','Efficiency vs Rotation Speed',    'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(6)       = uipanel('parent',Plots,'title','On-design Performance',           'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','on');

    PlotPanels(7)       = uipanel('parent',Plots,'title','Circulation Distribution',        'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(8)       = uipanel('parent',Plots,'title','Induced Velocity',                'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(9)       = uipanel('parent',Plots,'title','Inflow Angle',                    'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(10)       = uipanel('parent',Plots,'title','Expanded Blade (as designed)',   'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(11)       = uipanel('parent',Plots,'title','Blade Thickness (as designed)',  'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(12)       = uipanel('parent',Plots,'title','Lift Coefficient',               'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(13)       = uipanel('parent',Plots,'title','2D Geometry',                    'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(14)       = uipanel('parent',Plots,'title','3D Geometry',                    'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(15)       = uipanel('parent',Plots,'title','Performance Curves',             'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    % ---------------------------------------------------------------------

    % =========================================================================
    % ------------------------ Togglebar Panel Elements -----------------------

    ToggleText(1)	= uicontrol(Togglebar,'units','characters','style',...
                                'text','string','From Inputs:','position',[1 Toggleboxht-2-pushht textbox textht]);

    Toggle(1)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Expanded Blade','position',[1 Toggleboxht-2-pushht*2 textbox pushht]);

    Toggle(2)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Blade Thickness','position',[1 Toggleboxht-2-pushht*3 textbox pushht]);

    Toggle(3)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Inflow Profile','position',[1 Toggleboxht-2-pushht*4 textbox pushht]);

    ToggleText(2)	= uicontrol(Togglebar,'units','characters','style',...
                                'text','string','Parametric:','position',[1 Toggleboxht-2-pushht*5 textbox textht]);

    Toggle(4)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Effy vs D','position',[1 Toggleboxht-2-pushht*6 textbox pushht]);

    Toggle(5)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Effy vs N','position',[1 Toggleboxht-2-pushht*7 textbox pushht]);

    ToggleText(3)	= uicontrol(Togglebar,'units','characters','style',...
                                'text','string','Design results:','position',[1 Toggleboxht-2-pushht*8 textbox textht]);

    Toggle(6)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Design Performance','position',[1 Toggleboxht-2-pushht*9 textbox pushht],'value',1);

    Toggle(7)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Circulation Distribution','position',[1 Toggleboxht-2-pushht*10 textbox pushht]);

    Toggle(8)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Induced Velocity','position',[1 Toggleboxht-2-pushht*11 textbox pushht]);

    Toggle(9)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Inflow Angle','position',[1 Toggleboxht-2-pushht*12 textbox pushht]);

    Toggle(10)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Expanded Blade','position',[1 Toggleboxht-2-pushht*13 textbox pushht]);

    Toggle(11)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Blade Thickness','position',[1 Toggleboxht-2-pushht*14 textbox pushht]);

    Toggle(12)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Lift Coefficient','position',[1 Toggleboxht-2-pushht*15 textbox pushht]);

    Toggle(13)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','2D Geometry','position',[1 Toggleboxht-2-pushht*16 textbox pushht]);

    Toggle(14)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','3D Geometry','position',[1 Toggleboxht-2-pushht*17 textbox pushht]);

    Toggle(15)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Performance Curves','position',[1 Toggleboxht-2-pushht*18 textbox pushht]);

    % ---------------------------------------------------------------------


    % =========================================================================
    % ------------------ On Design Performance Panel Elements -----------------

    ODtextsrc           = {'J =' 'KT =' 'KQ =' 'EFFY =' 'ADEFFY =' 'CT =' 'CQ =' 'CP ='};

    for index = 1 : length(ODtextsrc)

        OnDesignText(index)     = uicontrol(PlotPanels(6),'units','characters','style',...
                                'text','string',ODtextsrc(index),'position',...
                                [2 Displayboxht-4-pushht*index*1.25 ondesignbox textht],...
                                'horizontalalignment','left');

        OnDesignValues(index)	= uicontrol(PlotPanels(6),'units','characters','style',...
                                'edit','string','','position',...
                                [2+ondesignbox Displayboxht-4-pushht*index*1.25 editbox editboxht],...
                                'enable','off');

    end


    % Create a subpanel to group those values that may be converted between SI and English

    ConversionPanel     = uibuttongroup(PlotPanels(6),'units','characters','position',...
                                        [2+ondesignbox+editbox+6 Displayboxht-4-1-pushht*8*1.25...
                                        Conversionbox Conversionboxht],'clipping','on',...
                                        'selectionchangedfcn',@convertfn);

    systemToggle(1)     = uicontrol(ConversionPanel,'style','radiobutton','units','characters',...
                                    'position',[2 Conversionboxht-2.5 editbox textht],...
                                    'string','SI','horizontalalignment','left','enable','off');

    systemToggle(2)     = uicontrol(ConversionPanel,'style','radiobutton','units','characters',...
                                    'position',[2+editbox Conversionboxht-2.5 editbox textht],...
                                    'string','English','horizontalalignment','left','enable','off');


    Conversiontextsrc   = {'Ship speed (Vs) =' 'Rotation speed (N) =' 'Diameter (D) ='...
                           'Thrust (T) =' 'Torque (Q) =' 'Power (P) ='};

    unitssrc            = {'m/s' 'RPM' 'm' 'N' 'Nm' 'W'};

    for index = 1 : length(Conversiontextsrc)

        ConversionText(index)	= uicontrol(ConversionPanel,'units','characters','style',...
                                'text','string',Conversiontextsrc(index),'position',...
                                [2 Conversionboxht-2.5-pushht*(index)*1.25 textbox textht],...
                                'horizontalalignment','left');

        ConversionValues(index)	= uicontrol(ConversionPanel,'units','characters','style',...
                                'edit','string','','position',...
                                [2+textbox Conversionboxht-2.5-pushht*(index)*1.25 ondesignbox editboxht],...
                                'enable','off');

        UnitsText(index)	= uicontrol(ConversionPanel,'units','characters','style',...
                                'text','string',unitssrc(index),'position',...
                                [2+textbox+ondesignbox+1 Conversionboxht-2.5-pushht*(index)*1.25 editbox textht],...
                                'horizontalalignment','left');



    end


end


