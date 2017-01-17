% Part of the code used in:
% Weitz et al. Lysis, Lysogeny, and Virus-Microbe Ratios
% 
% From https://github.com/WeitzGroup/VMR-Lysis-Lysogeny-v3
% MIT License

function psprintc(filename)
% psprint('filename')
% prints current figure with 
% the following flags
% -depsc2
% .ps will be appended
%
% print in color
filenameps = sprintf('%s.ps',filename);
print(filenameps,'-depsc2');
filenamjpg = sprintf('%s.jpg',filename);
print(filenamjpg,'-djpeg80');
end
