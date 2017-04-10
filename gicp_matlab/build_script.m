close all
clear all
clc

save_cd = pwd(); 
cd(fileparts(mfilename('fullpath'))); 
mex -setup:mex_C++_mingw-w64.xml C++

include_str = ['-I', getenv('MINGWROOT'), '\include '];
out_dir_str = '-outdir "..\" ';
file_name_str = '-output gicp_alg ';
build_command = ['mex ', include_str, out_dir_str, file_name_str, '*.cpp *.c'];
   
disp(build_command)
eval(build_command)

cd(save_cd)