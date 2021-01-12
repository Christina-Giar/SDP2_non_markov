function x = traceoutcom(W, sys, dim)
% W is the matrix, sys the system whose complementary will be traced out, dim is the dimension [d_in d_out d_in d_out etc.]
% output is normalized by dividing the dimension of the output systems that were traced out.

% Author: Christina Giarmatzi
% Requires: TrX
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

sys = sort(sys);
s = 1;
i = 1;
product = 1;
sysT = zeros(1, length(dim) - length(sys));
while s <= length(dim);
    for o = 1:length(sys);
        if s == sys(o);
            s = s+1;
        end
    end
    if s <= length(dim);
        sysT(i) = s;
        if ~mod(s,2); % If s is even, then s is an output system and its dimension adds to the normalization.
        product = product*dim(s);
        end
        i = i+1;
        s = s+1;
        
    end
end
x = TrX(W, sysT, dim)/product;% TrX just traces out normalizing by the product of the dimensions of output systems