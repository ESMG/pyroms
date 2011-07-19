% makefile for building mxlpsolve for lpsolve VERSION 5.5
%
% Usage: Makefile({what {, debug {, wait}}})
%  what: '' or not provided: compile build all source code
%        else any combination of lpsolve.c, matlab.c hash.c
%  debug: 0 or not provided: no debug info
%         1: argument checking and print extra debug info.
%  wait: 1 or not provided: ask to press a key before starting compiling
%        0: immediately start compiling

function Makefile(what, debug, wait)

windir=getenv('windir');
if isempty(windir) | length(windir) == 0
	lpsolvepath='../../../lib';
	objext='.o';
%	 libs=['-cxx -L' lpsolvepath '/lpsolve55 -llpsolve55 -lm'];
        libs=['-cxx -L' lpsolvepath '/lpsolve55 -lm'];
else
	lpsolvepath='..\..\..';
	objext='.obj';
	libs='lp_explicit.lib';
end

	lpsolvelibpath='';

if nargin == 0 | isempty(what)
        what = 'lpsolve.c matlab.c hash.c';
end

if nargin < 2
        debug = 0;
end

if debug == 0
        debug = '';
else
        debug = '-argcheck -DDEBUG'
end

if nargin < 3
        wait = 1;
end

disp('Automatic compilation of Matlab MEX interface for lp_solve 5.5');
disp(' ');
disp(['We assume that lp_solve 5.5 is installed in: ' lpsolvepath]);
disp('If that is not correct path for lp_solve 5.5, modify this file accordingly.');
disp(' ');
disp('***************************************************************');
disp('*  Old version of MEX lp_solve files will be overwritten !!!  *');
disp('***************************************************************');
disp(' ');
if wait == 0
else
	disp('Press any key to continue');
	disp(' ');
	pause;
end
disp('Compiling in progress. Please wait...');

% compile lp that uses lp_solve 5.5
eval(['mex ' debug ' -D_WINDOWS -DMATLAB -DINLINE=static -I' lpsolvepath ' ' lpsolvelibpath ' -c ' what]);
eval(['mex ' ' lpsolve' objext ' matlab' objext ' hash' objext ' ' libs ' -output mxlpsolve']);

disp('Compiling finished.');
disp('You can now run ex and lpdemo on matlab command line to test.');
