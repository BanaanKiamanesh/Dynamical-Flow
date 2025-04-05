clear;
close all;
clc;

%% Parameter Definition
ASize = 50;

% % Create A Matrix as invT * A * T
% A = diag(randperm(100, ASize));
% T = orth(rand(ASize));
% A = T' * A * T;

% Create a Randomized A Matrix
A = randn(ASize);
A = (A + A')/2;

N = diag(1:ASize);

%% ODE Solve
H0 = A(:);
tSpan = linspace(0, 3, 10000);
opt = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, "OutputFcn", @(t, y, flag) ProgressBar(t, y, flag, tSpan));
[~, Hsol] = ode45(@(t, Hvec) ODEFun(t,Hvec,ASize,N), tSpan, H0, opt);

%% Final Diagonalized Matrix
Hfinal = reshape(Hsol(end, :), ASize, ASize);

%% Visualization
videoFile = VideoWriter('H_Evolution.mp4','MPEG-4');
videoFile.FrameRate = 30;
open(videoFile);

figure("Name", "Plot", "Units", "normalized", "OuterPosition", [0, 0, 1, 1]);
h = imagesc;
axis equal
axis tight
axis off
colormap('jet')

% Record the Transient with Higher FPS
Start = 1;
End   = (size(Hsol, 1)/100);
Step  = 1;
for i = Start:Step:End
    Data = reshape(Hsol(i, :), ASize, ASize);
    h.CData = fliplr(Data);

    drawnow limitrate

    frame = getframe(gcf);
    writeVideo(videoFile, frame);

    disp(i)
end

% Record the Convergence with Lower FPS
Start = round(size(Hsol, 1)/100);
End   = size(Hsol, 1);
Step  = 20;
for i = Start:Step:End
    Data = reshape(Hsol(i, :), ASize, ASize);
    h.CData = fliplr(Data);

    drawnow limitrate

    frame = getframe(gcf);
    writeVideo(videoFile, frame);

    disp(i)
end

close(videoFile);

%% Result Check
fprintf('Original Matrix:\n');
disp(A);

fprintf('Diagonalized Matrix:\n');
disp(Hfinal);

fprintf('Eigenvalues of Diagonalized Matrix:\n');
disp(sort(eig(Hfinal))');

fprintf('Eigenvalues of Original Matrix:\n');
disp(sort(eig(A))');

%% ODE Function Definition
function dHvec = ODEFun(~, Hvec, ASize, N)
    H = reshape(Hvec, ASize, ASize);
    Hdot = H^2*N - 2*H*N*H + N*H^2;
    dHvec = Hdot(:);
end

function status = ProgressBar(t, ~, flag, tSpan)
    status = 0;
    persistent wb

    switch flag
        case 'init'
            wb = waitbar(0, "Initiating the ODE Solve");
        case ''
            % Calc Percent
            Percent = t(1) / tSpan(end);

            % Set to the Wb
            waitbar(Percent, wb, "ODE Solve in Progress.");
        case 'done'
            waitbar(1, wb, "ODE Solve Terminated!");

            pause(1)
            close(wb)
    end
end
