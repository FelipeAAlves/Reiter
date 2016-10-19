function dyntable(title,headers,labels,values,label_width,val_width, ...
                  val_precis)

% Copyright (C) 2002-2010 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global options_

if options_.noprint
    return
end

label_width = max(size(deblank(char(headers(1,:),labels)),2)+2, ...
                  label_width);
val_width = max(size(deblank(headers(2:end,:)),2)+2,val_width);
label_fmt = sprintf('%%-%ds',label_width);
header_fmt = sprintf('%%-%ds',val_width);
val_fmt = sprintf('%%%d.%df',val_width,val_precis);
if length(title) > 0
    disp(sprintf('\n\n%s\n',title));
end
if length(headers) > 0
    hh = sprintf(label_fmt,headers(1,:));
    hh = [hh char(32*ones(1,floor(val_width/4)))];
    for i=2:size(headers,1)
        hh = [hh sprintf(header_fmt,headers(i,:))];
    end
    disp(hh);
end
for i=1:size(values,1)
    disp([sprintf(label_fmt,labels(i,:)) sprintf(val_fmt,values(i,:))]);
end

% 10/30/02 MJ