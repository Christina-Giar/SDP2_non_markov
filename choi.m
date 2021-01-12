 function x = choi(U)
% Takes unitary U of any size and outputs it's Choi-Jamiolkowski representation,
% without the transpose

% Author: Christina Giarmatzi
% Requires: proj
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

U = reshape(U',1,size(U,1)*size(U,2));
x = proj(U);
