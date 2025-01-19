% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function changeChord(hObject,ED)
    % global SpecificationsValues DuctValues CavValues FlagValues FoilValues Filename...
    %        XR_in XCoD_in XCD_in VAI_in VTI_in Xt0oD_in skew0_in rake0_in...
    %        Meanline_cell Thickness_cell;

    global FlagValues XCoD_in;

    global Col_Label XCoD_values XCLmax_values;

    if get(FlagValues(5),'value')

        XCoD_values     = str2double(get(XCoD_in,'string'));

        set(Col_Label(2),'string','XCLmax');

        for i = 1:length(XCoD_in)
            set(XCoD_in(i),'string',num2str(XCLmax_values(i)));
        end
    else

        XCLmax_values  = str2double(get(XCoD_in,'string'));


        set(Col_Label(2),'string','c/D');

        for i = 1:length(XCoD_in)
            set(XCoD_in(i),'string',num2str(XCoD_values(i)));
        end
    end

end


