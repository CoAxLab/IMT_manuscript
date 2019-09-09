function mywd = cdtodrive(whichdrive)

% change directory to drive of interest
% if it is not connected, change directory to the next on the list
% thus anything done on TEK external can be replicated on ProjectDrive by any user

if nargin < 1, whichdrive = 'ProjectDrive'; end

% also now accompanies different users accessing ProjectDrive via the linux machine
% ADD YOUR DRIVE HERE (but keep '/Volumes/ProjectDrive' in the bottom spot)
ProjectDriveDirs = ...
    {'/home/tek31/ProjectDrive' ... 
    '/home/chk65/ProjectDrive' ... 
    '/home/jrd177/ProjectDrive' ...
    '/Volumes/ProjectDrive'};

fprintf('\n=======================================\n')

switch whichdrive
    case {'tekexternal'}
        drive='/Volumes/THOMAS/pitt/';
        if exist(drive, 'dir')
            cd(drive)
            mywd = pwd;
        else
            fprintf('Could not find TEK external drive');
        end
    otherwise
        for p = 1:length(ProjectDriveDirs)
            if exist(ProjectDriveDirs{p}, 'dir')
                cd(ProjectDriveDirs{p})
                mywd = pwd;
                fprintf('CHANGED DIRECTORY TO:\n%s', ProjectDriveDirs{p})
                fprintf('\n=======================================\n')  
                return;
            end
        end
end

% If none of the drives can be found, report an error.
fprintf('Could not find any of the listed drives')

%
fprintf('\n=======================================\n')              
%             
%                 
% if ~exist(drive,'dir')
%     drive='/Volumes/ProjectDrive/Thomas';
% end
% if ~exist(drive,'dir')
%    
%    return;
% end
% cd(drive);
% fprintf('\n\nCHANGED DIRECTORY TO:\n%s\n\n',drive)
% mywd = pwd;