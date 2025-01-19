% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function convertfn(hObject,ED)

    global systemToggle ConversionValues UnitsText;

    global pt;

    if get(systemToggle(1),'value')

        set(ConversionValues(1),'string',num2str(pt.input.Vs));
        set(ConversionValues(2),'string',num2str(pt.input.N));
        set(ConversionValues(3),'string',num2str(pt.input.D));
        set(ConversionValues(4),'string',num2str(pt.input.THRUST));
        set(ConversionValues(5),'string',num2str(pt.design.Q));
        set(ConversionValues(6),'string',num2str(pt.design.P));

        unitssrc            = {'m/s' 'RPM' 'm' 'N' 'Nm' 'W'};

        for index = 1 : length(unitssrc)

            set(UnitsText(index),'string',unitssrc(index));

        end

    else

        set(ConversionValues(1),'string',num2str(pt.input.Vs*1.94384449));     % Convert to knots
        set(ConversionValues(2),'string',num2str(pt.input.N));
        set(ConversionValues(3),'string',num2str(pt.input.D*3.2808399));       % Convert to feet
        set(ConversionValues(4),'string',num2str(pt.input.THRUST*0.224808943));     % Convert to lbf
        set(ConversionValues(5),'string',num2str(pt.design.Q*0.737562149277));  % Convert to lb ft
        set(ConversionValues(6),'string',num2str(pt.design.P*0.00134102209));   % Convert to Hp

        unitssrc            = {'knots' 'RPM' 'ft' 'lb' 'lb/ft' 'HP'};

        for index = 1 : length(unitssrc)

            set(UnitsText(index),'string',unitssrc(index));

        end

    end

end
