 =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function changeDir(hObject,ED)

    global Filename OpenPropDirectory;

    filename    = get(Filename,'string');                       % Filename prefix

    % -------------------------------------------------------------------------
    % Figure out what the current working directory is,
    % and change directories to OpenPropDirectory/filename
    rest    = pwd;
    root  = '';

    while ~isempty(rest)
        [CurrentDirectory,rest] = strtok(rest,'/');

        if strcmp(CurrentDirectory,OpenPropDirectory)

            if isempty(rest)

                % you are in /OpenPropDirectory/
                %             addpath ../SourceCode

            elseif strcmp(rest(2:end),filename)
                % already in /OpenPropDirectory/filename/
                %             addpath ../SourceCode
                rest = [];

            elseif strcmp(rest(2:end),'SourceCode')

                %             addpath ../SourceCode
                rest = [];

            else
                % you are in /OpenPropDirectory/wrongfolder
    %             disp('NOTICE: Must start OpenProp from the root directory.')

                cd([root '/' OpenPropDirectory]);

    %             disp('Changed to root directory.')

            end
        end

        root  = [root '/' CurrentDirectory];

    end
    % -------------------------------------------------------------------------
end

