function x = evolutionW(H,t)
% Evolution of unitary H for time t. Exp is a matrix function.

% Author: Christina Giarmatzi
% Requires: TrX, choi, syspermute
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

k = funm(-(1i*H.*t), @exp);
W = TrX(choi(k), 4, [2 2 2 2]);
x = syspermute(W, [2 1 3], [2 2 2]); % So that W has AIAOBI ordering.