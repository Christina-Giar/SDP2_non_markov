function [x, dist, witness] = SDP2witness(W)
% SDP for detecting quantum memory in the input 2-timestep process matrix W.
% The input is W^{AiAoBi} with dimensions [2 2 2] (size of W is 8x8).
% Uses the PPT criterion on an extension of input W(A_iA_oB_i) with Ai. 
% Uses the code of SDP2witness_W.m which has examples of input W.
% Outputs:
% dist is the norm of (W - Wmarkovian) as a measure of non-Markovianity.
% x is real( Tr(W Witness))  which is another measure of non-Markovianity.
% witness is the Witness of qunatum memory for the input W.

% Author: Christina Giarmatzi
% Requires: TrX, tensor, sysexchange, syspermute, PartialTranspose, 
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

dlist = [2 2 2]; % Dimensions of the W subsystems A_i A_o B_i.
dim = dlist;
N = 2;
epsilon = 10e-31; % Accuracy on the output value real(trace('witness'*'process')).

% Pauli basis.
I2 = eye(2)/(2);
sigX = [0, 1 ; 1, 0]/(2);
sigY = [0, -1i ; 1i, 0]/(2);
sigZ = [1, 0 ; 0, -1]/(2);
sigs = {I2, sigX, sigY, sigZ};

% Constructing the dimensions of W_ext (the extension of W or rho) 
% [input channel channel ... output] such that we have rho^ABC... 
% to perfrom the swap, traceout, and partial transpose operations.
% Note that the last party has only input.
dim_interm = zeros(1,N-1);
for i = 1:N-1
    dim_interm(i) = dim(2*i)*dim(2*i+1);% i.e. Ao*Bi
end
% Here the extension is the subsystem Ai.
dim_ext = [dim(1) dim_interm dim(1)];% i.e. [Ai AoBi Ai]
d = size(W,1)*2; % the dimension of the W_ext=dim(W)*(d of the system we extend) here A_I

% Define the G0 term through W. G is the extension of W and lives on AiAoBiAi.
presum1 = zeros(d, d, 16);
for i = 1:4
    for j = 1:4
        W_1ij = 2*trace(W*(tensor(sigs{1},sigs{i},sigs{j})));
        presum1(:,:,(i-1)*4+j) = W_1ij*tensor(sigs{1}, sigs{i}, sigs{j}, sigs{1}); % AiAoBiAi
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

% The SDP. Write cvx_begin quiet sdp to suppress the output "Calling SDP3..."
% tic % Uncomment if you want to time the SDP.
cvx_begin sdp
    variable Z(3*d,3*d) hermitian semidefinite
    variable Z0(d,d) hermitian semidefinite
    variable Z1p(d,d) hermitian semidefinite
    variable Z2p(d,d) hermitian semidefinite
    maximize -real(trace(F0*Z))
    subject to
    
    Z == blkdiag(Z0,Z1p,Z2p);

    for i = 2:4
        for j = 1:4
            for n = 1:4
                for k = 2:4
                    term1 = trace(tensor(sigs{i}, sigs{j}, sigs{n}, sigs{k})*Z0);
                    term2 = trace(PartialTranspose(tensor(sigs{i}, sigs{j}, sigs{n}, sigs{k}), 1, dim_ext)*Z1p);
                    term3 = trace(PartialTranspose(tensor(sigs{i}, sigs{j}, sigs{n} , sigs{k}), 2, dim_ext)*Z2p);
                    term1 + term2 + term3 == 0;
                end
            end
        end
    end

    z = Z0 + PartialTranspose(Z1p, 1, dim_ext) + PartialTranspose(Z2p, 2, dim_ext);
    di = 2;% Dimension of the system we used to extend the W, here Ai.
    la1 = 1/di*TrX(z, 3, dim_ext);
    la2 = 1/di*TrX(sysexchange(z, [1 3], dim_ext), 3, dim_ext);
    la3 = (1/di^2)*tensor(eye(di), TrX(z, [1 3], dim_ext));
    wit = la1 + la2 - la3;
    trace(wit)==1;

% An extra constraint on the witness if anyone wants: Alice measures psi and prepares psi.
% wit == syspermute(wit, [2 1 3], [2 2 2])

cvx_end
% toc % Uncommment to time the SDP.

% Assign the output variables.
Tr_W_input = real(trace(wit*W));
x = Tr_W_input;
witness = wit;
% For the distance to the Markovian one.
% All trace functions (which are based on TrX) rely on the fact that 
% dim = [input output input output etc] and do not deal well with output = 1.
% Hence I add identity on Bo of dimension 2
d_bo = 2;
Bo = eye(d_bo)/2;
Wnew = tensor(W, Bo*2); % Now Wnew is with Ai Ao Bi Bo
dimnew = [dim d_bo];
Wmarkov = tensor(traceoutcom(Wnew, 1, dimnew), traceoutcom(Wnew, [2 3], dimnew), Bo);
dist = norm(Wnew - Wmarkov);

% Output on command window the result
if Tr_W_input < -epsilon
    cprintf('*blue', 'Quantum memory detected :).');
    fprintf('\n')
else
     cprintf('blue', 'Quantum memory not detected :(.');
     fprintf('\n')
end
    