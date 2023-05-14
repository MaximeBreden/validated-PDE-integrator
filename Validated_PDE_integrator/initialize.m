% this is an initialization file
% please adapt the path appropriately

intlabpath = [];

% If you have intlab provide the path to it below.
% Comment out the next line if you don't have intlab

intlabpath = '~/matlab/intlab/Intlab_V12/';

if ~exist('intval','file') && ~isempty(intlabpath)  
    dir = pwd;
    cd(intlabpath)
    startintlab;
    cd(dir)
end

clearvars;

global intervalarithmeticavailable
if exist('intval','file')
    intervalarithmeticavailable = true;
    setround(0);
else
    % the code will run without interval arithmetic
    % but will then only produce a "proof" modulo rounding errors
    intervalarithmeticavailable = false;
end


