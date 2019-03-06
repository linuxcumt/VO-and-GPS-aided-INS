function pathsetup()
%% pathsetup
% 
% Adds folders and subfolders from the include directory to the matlab
% path.
% 
% @author: Matt Marti
% @date: 2019-03-05

clear, clc, clear global

recursiveadd('include');

end


function recursiveadd(str)
% Recursively add all directories in str to the path
addpath(str);
dirstruct = dir(str);
for i = 3:length(dirstruct)
    d = dirstruct(i);
    dirname = [str '/' d.name];
    if exist(dirname, 'dir')
        recursiveadd( dirname );
    end
end

end