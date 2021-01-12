function x = proj(a)
% turns a state into a projector

% Author: Christina Giarmatzi
% Requires: tensor
% License: GPL3

%     Copyright (C) 2020  Dr Christina Giarmatzi, christina.giar@gmail.com
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.


x = tensor(a, a');