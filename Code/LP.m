clear;
close all;
clc;

%% Parameter Definition
c = [1 2 3 4 5];
Vertices = eye(5);
numVertices = size(Vertices,1);

% Add Explicit Constraints (Ax <= b)
A = [ 1 1  0 0 0
    0 0  1 1 1
    -1 0 -1 0 0];

b = [1; 2; 0];

% Filter Vertices Based on Constraints
ValidVertices    = Vertices(all(A*Vertices' <= b, 1), :);
NumValidVertices = size(ValidVertices,1);

% Update Matrices Based on Valid Vertices
Q = diag([1 zeros(1,NumValidVertices-1)]);
N = diag(ValidVertices*c');

% Initial Orthogonal Matrix Theta
Theta = orth(randn(NumValidVertices));
H     = Theta'*Q*Theta;

%% Dynamical Parameters
StepSize   = 0.001;
Iterations = 5000;

%% ODE Solve
for i = 1:Iterations
    Hdot = H*H*N - 2*H*N*H + N*H*H;
    H    = H + StepSize*Hdot;
end

%% Results
% Optimal Vertex
[~, idx]      = max(diag(H));
OptimalVertex = ValidVertices(idx, :);
OptimalValue  = c * OptimalVertex';

fprintf('Optimal Vertex: [%g, %g, %g, %g, %g]\n', OptimalVertex);
fprintf('Optimal Value: %g\n', OptimalValue);
