% Dual SDP to find a witness for a quantum non-Markovian process.
% There are several ways to create the extension of W (W_ext), depending on the partition of W.
% Here we assume that subsystem A is A_I and this is what we extend the input W. i.e. tensor(W,A_I)

% Author: Christina Giarmatzi
% Requires: evolutionW, choi, TrX, proj, tensor, sysexchange, syspermute, PartialTranspose, 
%           opt_args, PermuteSystems, traceoutcom, cprintf
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
sigX = [0, 1 ; 1, 0]/(2);
sigY = [0, -1j ; 1j, 0]/(2);
sigZ = [1, 0 ; 0, -1]/(2);
tag = 'evolution';
t = 1;% For t=0 the process is Markovian. t=1 has quantum memory.
H1 = tensor(sigZ, sigZ);
H2 = -(tensor(sigX, sigX) + tensor(sigY, sigY) + tensor(sigZ, sigZ));
W = evolutionW(H2,t); % evolution(hamiltonian, time)
dlist = [2 2 2]; %input the dimensions of the W subsystems A_i A_o B_i
dim = dlist;
N = (size(dim,2)+1)/2;% Number of parties.
epsilon = 10e-31;




% Pauli basis
I2 = eye(2)/(2);
sigX = [0, 1 ; 1, 0]/(2);
sigY = [0, -1j ; 1j, 0]/(2);
sigZ = [1, 0 ; 0, -1]/(2);
sigs = {I2, sigX, sigY, sigZ};

% Constructing the dimensions [input channel channel ... output] to be used with the PartialTranspose function
% Note that the last party has only input.
dim_interm = zeros(1,N-1);
for i = 1:N-1
    dim_interm(i) = dim(2*i)*dim(2*i+1);% i.e. Ao Bi
end
% Assuming the extension is a party with input output dimension of 2
dim_ext = [dim(1) dim_interm dim(1)];% i.e. [Ai AoBi Ai]
d = size(W,1)*2; % the dimension of the W_ext=dim(W)*(d of the system we extend) here A_I

% Define the G0 term through W.


% Note that this is for a bipartite state. For a tripartite we would need i,j,k,l, etc..
presum1 = zeros(d, d, 16);
for i = 1:4
    for j = 1:4
        W_1ij = 2*trace(W*(tensor(sigs{1},sigs{i},sigs{j})));
        presum1(:,:,(i-1)*4+j) = W_1ij*tensor(sigs{1}, sigs{i}, sigs{j}, sigs{1});
    end
end
sum1 = sum(presum1, 3); % for efficiency, this could be done in the loop
sum2 = zeros(d, d);
for k = 2:4
    for i = 1:4
        for j = 1:4
            W_kij = 2*trace(W*tensor(sigs{k}, sigs{i}, sigs{j}));
            sum2 = sum2 + W_kij*(tensor(sigs{k}, sigs{i}, sigs{j}, sigs{1}) + tensor(sigs{1}, sigs{i}, sigs{j}, sigs{k}));
        end
    end
end
G0 = sum1 + sum2;

F0 = blkdiag(G0,PartialTranspose(G0, 1, dim_ext),PartialTranspose(G0, 2, dim_ext));

tic
cvx_begin sdp
    variable Z(3*d,3*d) hermitian semidefinite
    variable Z0(d,d) hermitian semidefinite
    variable Z1p(d,d) hermitian semidefinite
    variable Z2p(d,d) hermitian semidefinite
    maximize -real(trace(F0*Z))
%     maximize -real(trace(G0*(Z0+Z1+Z2)))
    subject to
    Z == blkdiag(Z0,Z1p,Z2p);
    % Z1p = PartialTranspose(Z1, 1, dim_ext);
    % Z2p = PartialTranspose(Z2, 2, dim_ext);
for i = 2:4
    for j = 1:4
        for n = 1:4
        for k = 2:4
%             trace(tensor(sigs{i1}, sigs{j1}, sigs{n1}, sigs{k1})*(Z0 + Z1 + Z2)); % For some reason this does not work.
            term1 = trace(tensor(sigs{i}, sigs{j}, sigs{n}, sigs{k})*Z0);
            term2 = trace(PartialTranspose(tensor(sigs{i}, sigs{j}, sigs{n}, sigs{k}), 1, dim_ext)*Z1p);
            term3 = trace(PartialTranspose(tensor(sigs{i}, sigs{j}, sigs{n} , sigs{k}), 2, dim_ext)*Z2p);
            term1 + term2 + term3 == 0;
        end
        end
    end
end

z = Z0 + PartialTranspose(Z1p, 1, dim_ext) + PartialTranspose(Z2p, 2, dim_ext);
% z = Z0 + Z1 + Z2;
% z = Z0 + Z1p + Z2p;
di = 2;% This is the dimension of the system we used to extend the W, here Ai.
la1 = 1/di*TrX(z, 3, dim_ext);
la2 = 1/di*TrX(sysexchange(z, [1 3], dim_ext), 3, dim_ext);
la3 = (1/di^2)*tensor(eye(di), TrX(z, [1 3], dim_ext));
wit = la1 + la2 - la3;
trace(wit)==1;

cvx_end
toc

Tr_W_input = real(trace(wit*W))

if Tr_W_input < -epsilon
    cprintf('*blue', 'Quantum memory detected in the process #%s',tag);
    fprintf('\n')
else
     cprintf('blue', 'Quantum memory not detected in the process #%s',tag);
     fprintf('\n')
end
    



        