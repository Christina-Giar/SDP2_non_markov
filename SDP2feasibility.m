function x = SDP2feasibility(W)
% SDP for checking entanglement in the 2-timestep input W, seen as a bipartite state. 
% Called the feasibility problem, using PPT criterion to an extension of W^{AiAoBi} with Ai.
% The input is W^{AiAoBi} with dimensions [2 2 2] (size of W is 8x8).
% The ouput x is 0 if the state is entangled, i.e. quantum memory has been
% detected. Else, x = 1. Uses the code SDP2feasibility_W.m which has examples
% of input W.

% Author: Christina Giarmatzi
% Requires: tensor, sysexchange, syspermute, PartialTranspose, opt_args, 
%           PermuteSystems, cprintf
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
    
    
dlist = [2 2 2]; % Dimensions of the W subsystems A_i A_o B_i
dim = dlist;
N = 2;  % Number of parties in the process W.
% Constructing the dimensions of W_ext (the extension of W or rho) 
% [input channel channel ... output] such that we have rho^ABC... 
% to perfrom the swap, traceout, and partial transpose operations.
% Note that the last party has only input.
n_arr = N-1; % Parties are labelled such that arrows are from A to B, B to C etc.
dim_interm = zeros(1,n_arr);
for i = 1:n_arr
    dim_interm(i) = dim(2*i)*dim(2*i+1);
end
% Here the extension is the subsystem Ai.
dim_W_ext = [dim(1) dim_interm dim(1)];
d_ext = size(W,1)*dim(1); % Dimension of the W_ext.

% The SDP. Write cvx_begin quiet sdp to suppress the output "Calling SDP3..."
cvx_begin quiet sdp
    variable W_ext(d_ext,d_ext) semidefinite;
    minimize 0
    subject to 
    % W_ext is a symmetric extension of W iff:
    TrX(W_ext, 3, dim_W_ext) == W; %  Trace out the last system and end up with initial W.
    W_ext == sysexchange(W_ext, [1,3], dim_W_ext);% W_ext is symmetric under swap of systems 1 and 3.
    % rho and the partial transposes are positive
    blkdiag(W_ext,PartialTranspose(W_ext, 1, dim_W_ext),PartialTranspose(W_ext, 2, dim_W_ext)) >= 0 ;
cvx_end

% Output on command window the result
if strcmp(cvx_status,'Infeasible')
    x=0;
    cprintf('*blue', 'Quantum memory detected :).');
    fprintf('\n')
elseif strcmp(cvx_status,'Solved')
    x=1;
     cprintf('blue', 'Quantum memory not detected :(.');
     fprintf('\n')
end


