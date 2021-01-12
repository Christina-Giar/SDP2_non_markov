% SDP for checking entanglement in the 2-timestep input W, seen as a bipartite state. 
% Called the feasibility problem, using PPT criterion to an extension of W^{AiAoBi} with Ai.
% The input is W^{AiAoBi} with dimensions [2 2 2] (size of W is 8x8).
% The ouput x is 0 if the state is entangled, i.e. quantum memory has been
% detected. Else, x = 1.

% Author: Christina Giarmatzi
% Requires: evolutionW, choi, TrX, proj, tensor, sysexchange, syspermute, 
%           PartialTranspose, opt_args, PermuteSystems, cprintf
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

clear all


% Example 1: Create an entangled and separable state.
% Pauli basis
% I2 = eye(2)/(2);
% sigX = [0, 1 ; 1, 0]/(2);
% sigY = [0, -1j ; 1j, 0]/(2);
% sigZ = [1, 0 ; 0, -1]/(2);
% ket0 = [1;0];
% ket1 = [0;1];
% norm = 1/sqrt(2);
% ket_phi = norm*(tensor(ket0,ket0) - tensor(ket1,ket1));
% phi = ket_phi*ket_phi';
% rho_ent = phi;
% rho_sep = tensor(ket0*ket0',ket0*ket0');
% % Use rho_ent or rho_sep to construct W.
% % W = syspermute(tensor(rho_ent, eye(2)), [1 3 2], [2 2 2]);
% W = syspermute(tensor(rho_sep, eye(2)), [1 3 2], [2 2 2]);
% dmat = [2 2; 2 1];
% N = 2;
% dlist = [2 2 2]; %input the dimensions of the W subsystems A_i A_o B_i
% dim = dlist;


% Example 2: Create an example (the example of the paper)
I2 = eye(2)/(2);
sigX = [0, 1 ; 1, 0]/(2);
sigY = [0, -1j ; 1j, 0]/(2);
sigZ = [1, 0 ; 0, -1]/(2);
H1 = tensor(sigZ, sigZ);
H2 = -(tensor(sigX, sigX) + tensor(sigY, sigY) + tensor(sigZ, sigZ));
t = 1;% For t=0 the process is Markovian. t=1 has quantum memory.
W = evolutionW(H2,t); % evolution(hamiltonian, time)
dlist = [2 2 2]; %input the dimensions of the W subsystems A_i A_o B_i
dim = dlist;
N = (size(dim,2)+1)/2;% Number of parties.
epsilon = 10e-30;


% Constructing the dimensions of W_ext (the extension of W/rho) 
% [input channel channel ... output] to be used with the transpose function
n_arr = N-1;% Parties are labelled such that arrows are from A to B, B to C etc.
dim_interm = zeros(1,n_arr);
for i = 1:n_arr;
    dim_interm(i) = dim(2*i)*dim(2*i+1);
end
% Assuming the extension is a party with input output dimension of 2
dim_W_ext = [dim(1) dim_interm dim(1)];
n = size(W,1)*dim(1);% The size of the W_ext


% The SDP. Write cvx_begin quiet sdp to suppress the output "Calling SDP3..."
cvx_begin sdp
    variable W_ext(n,n) semidefinite;
    minimize 0
    subject to
    % W_ext is a symmetric extension of W iff:
    TrX(W_ext, 3, dim_W_ext) == W; %  Trace out sys 3 and end up with initial W.
    W_ext == sysexchange(W_ext, [1,3], dim_W_ext);% W_ext is symmetric under swap of systems 1 and 3.
    % rho and the partial transposes are positive
    blkdiag(W_ext, PartialTranspose(W_ext, 1, dim_W_ext),PartialTranspose(W_ext, 2, dim_W_ext)) >= 0;   
cvx_end

if strcmp(cvx_status,'Infeasible')
    cprintf('*blue', 'Quantum memory detected :).');
    fprintf('\n')
elseif strcmp(cvx_status,'Solved')
     cprintf('blue', 'Quantum memory not detected :(.');
     fprintf('\n')
end