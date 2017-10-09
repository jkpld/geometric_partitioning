function setup

if verLessThan('matlab','9.1')
    error('geometricPartitioning:setup','The declumpNuclei code requires at least Matlab 2016b because implicit expansion is heavily used.')
end

fprintf('\nStarting setup...\n')

% Get the path to the current folder -------------------------------------
fileLocation = mfilename('fullpath');
path = fileparts(fileLocation);
fprintf('...found path\n')

% Compile the needed C files into .mex files -----------------------------
% Compile the needed C files into .mex files -----------------------------
have_compiler = true;
try 
    evalc('mex(''-setup'',''c'')');
catch % ME
    have_compiler = false;
    warning('geometricPartitioning:setup','There is no supported C compiler installed.\nThe code will still run; however, it could be slower for 2D data without the compiled mex functions.')
%     rethrow(ME)
end

if have_compiler
    try
%         if ~exist(fullfile(path,'private',['nakeinterp1.' mexext]),'file')
            mex(fullfile(path,'private','nakeinterp1.c'),'-outdir',fullfile(path,'private'),'-silent')
%         end
        fprintf('...Compiled nakeinterp1.c\n')
    catch % ME
    %     rethrow(ME)
        warning('geometricPartitioning:setup','There was an error compiling the required C function ''nakeinterp1.c''. Make sure that the function ''mex'' is coorectly setup to compile C code.\nThe code will still run; however, it could be slower without the compiled mex functions.')
    end
end

% Copy histcountsmex into the private folder ---------------------------
folder = fullfile(path,'private');
d = dir(folder);
names = {d.name};
% if ~any(strncmp('geoPart_histcountsmex',names,16))
    folder = fullfile(matlabroot,'toolbox','matlab','datafun','private');
    d = dir(folder);
    names = {d.name};
    nameIdx = strncmp('histcountsmex',names,13);

    if ~any(nameIdx)
        error('geometricPartitioning:setup','File ''histcountsmex'' was not found in\n %s\nThis file is required. Try manually locating it; if found copy into the private folder and rename it to ''geoPart_histcountsmex.(extension)''.',folder)
    end

    copyfile(fullfile(folder,names{nameIdx}),fullfile(path,'private',['geoPart_' names{nameIdx}]))
% end
fprintf('...added histcountsmex\n')

% Copy pdistmex into the private folder --------------------------------
% folder = fullfile(path,'private');
% d = dir(folder);
% names = {d.name};
% if ~any(strncmp('geoPart_pdistmex',names,11))
%     folder = fullfile(matlabroot,'toolbox','stats','stats','private');
%     d = dir(folder);
%     names = {d.name};
%     nameIdx = strncmp('pdistmex',names,8);
% 
%     if ~any(nameIdx)
%         error('declumpNuclei:setup','File ''pdistmex'' was not found in\n %s\nThis file is required. Try manually locating it; if found copy into the private folder and rename it to ''geoPart_pdistmex.(extension)''.',folder)
%     end
% 
%     copyfile(fullfile(folder,names{nameIdx}),fullfile(path,'private',['geoPart_' names{nameIdx}]))
% end
% fprintf('...added pdistmex\n')

% Add subfolders of current location to path -----------------------------
pths = split(string(genpath(path)),';');
toIgnore = {'.git','docs','unused'};
pths = pths(~pths.contains(toIgnore)).join(';');
addpath(pths.char());
fprintf('...Added subfolders to path\n')

fprintf('Setup finished!\n')
    
end
