function rootPath = banksRootPath()
% Return the path to the root of the Banks1991 directory
%
% This function must reside in the directory at the base of the
% Banks code directory structure.  It is used to determine the location of various
% sub-directories.
% 
% Example:
%   fullfile(banksRootPath,'scripts')

rootPath=which('banksRootPath');

rootPath=fileparts(rootPath);

return
