% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function changeViscous(hObject,ED)

    global FlagValues XCD_in XCD_values;

    if get(FlagValues(6),'value')

        set(XCD_in,'enable','on');

        for index = 1 : length(XCD_values)

            set(XCD_in(index),'string',num2str(XCD_values(index)));

        end

    else

        XCD_values  = str2double(get(XCD_in,'string'));

        set(XCD_in,'enable','off','string','0');

    end

end

