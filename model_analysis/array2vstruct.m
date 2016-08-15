% Part of the code used in:
% Weitz et al. Lysis, Lysogeny, and Virus-Microbe Ratios
% 
% From https://github.com/WeitzGroup/VMR-Lysis-Lysogeny-v2
% MIT License

function s = array2vstruct(x_array,fnames)
% function s = array2vstruct(x_array,fnames)
% This is a function which takes an array of values
% and places them in order into a particular structure
% in the order of the fieldnames in fnames.

if (length(fnames)~=length(x_array))
  error('Wrong number of elements to unpack in array2vstruct');
  exit(-1);
else
  tmpc=cell(1,length(x_array));
  for i=1:length(x_array),
    tmpc{i}=x_array(i);
  end
  s = cell2struct(tmpc,fnames,2);
end
end
