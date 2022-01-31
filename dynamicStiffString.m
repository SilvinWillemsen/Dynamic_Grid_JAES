%{
    Implementation of the dynamic stiff string accompanying the paper: 
    S. Willemsen, S. Bilbao, M. Ducceschi, and S. Serafin, "The Dynamic 
    Grid: Time-Varying Parameters for Musical Instrument Simulations based 
    on Finite-Difference Time-Domain Schemes," submitted to J. Audio Eng.
    Soc. (JAES), 2021.

    Equations referred to from this code can be found in the paper.

    The code generates the sound example "stringOctaveVibraMarimba.wav"
    found via https://tinyurl.com/DynSoundEx. Parameters are changed to 
    model the following:

    Stiff string -> stiff string with half length -> vibraphone -> marimba

    CC 4.0 Silvin Willemsen 2021.
%}

%% Reset workspace
clear all;
close all;
clc;

%% Initialise variables
fs = 44100;             % Sample rate (in Hz)
k = 1/fs;               % Time step (in s)

% The number of grid points from the right boundary the split happens (> 0)
numFromRightBound = 1;

if numFromRightBound < 1
    error("numFromRightBound is undefined")
end

%% Generate parameter vectors used for the simulation

%{
    The eventual simulation will be divided into chunks of static 
    parameters and chunks time-varying parameters. Each chunk is one second 
    long.

    One row in the 'params' matrix contains parameter values for
        - Length 'L' (in m)
        - Material density 'rho' (in kg/m^3)
        - Radius 'r' (in m) 
        - Tension 'T' (in N)
        - Young's Modulus 'E' (in Pa)
        - Frequency independent damping 'sig0' (in 1/s)
        - Frequency dependent damping 'sig1' (in m^2/s)

    Add rows for more changes (and more chunks). 

    Note that changes can not result in a change in grid configuration of 
    more than 2 grid points per sample.

%}
params = [1, 7850, 5e-4, 555, 2e11, 1, 0.005;
          0.5, 7850, 5e-4, 555, 2e11, 1, 0.005;
          0.5, 7850, 5e-4, 0, 7e13, 1, 0.005;
          0.5, 7850, 5e-4, 0, 7e13*1/4, 1.5, 0.05;]';

numParams = size(params, 1);    % number of parameters
numChanges = size(params, 2);   % number of changes between parameter settings

lengthSound = numChanges*2 * fs;  % length of the resulting simulation (in samples)
chunkSize = fs;                   % chunck length (one second)

% Create collection of vectors with parameter values. This will end up 
% being of size numparams x lengthSound. We start with a static chunk.
paramVecs = params(:, 1) * ones(1, chunkSize);

for i = 1:numChanges-1
    %{
        Generate chunk with time-varying parameters
    
        As changes in some parameters have a nonlinear effect on the 
        resulting pitch of the output sound, the change from one value to
        the next is altered such that the change in pitch is as linear as 
        possible
    %}
    paramMat = [];
    for j = 1:numParams
        if j == 1
            paramMat = [paramMat; (1./linspace(1/(params(j, i)), 1/(params(j, i+1)), chunkSize))];
        elseif j == 2
            paramMat = [paramMat; 1./linspace(1/sqrt(params(j, i)), 1/sqrt(params(j, i+1)), chunkSize).^2];
        elseif j == 4 || j == 5
            paramMat = [paramMat; linspace(sqrt(params(j, i)), sqrt(params(j, i+1)), chunkSize).^2];
        else
            paramMat = [paramMat; linspace(params(j, i), params(j, i+1), chunkSize)];
        end
    end
    % Add time-varying chunk to the collection of parameter vectors
    paramVecs = [paramVecs, paramMat];
      
    % Add static chunk 
    paramVecs = [paramVecs, params(:, i+1) .* ones(numParams, chunkSize)];
end

% End with extra static chunk
paramVecs = [paramVecs, params(:, end) .* ones(numParams, chunkSize)];

% Retrieve individual parameter vectors
Lvec = paramVecs(1, :);
rhoVec = paramVecs(2, :);
rVec = paramVecs(3, :);
Tvec = paramVecs(4, :);
Evec = paramVecs(5, :);
sig0Vec = paramVecs(6, :);
sig1Vec = paramVecs(7, :);

% Calculate wave speed and stiffness vectors used in the simulation
cSqVec = Tvec ./ (rhoVec .* pi .* rVec.^2);
kappaSqVec = (Evec .* pi / 4 .* rVec.^4) ./ (rhoVec .* pi .* rVec.^2);

% Grid spacing (in m) calculated by satisfying the stability condition for 
% the stiff string with equality
hVec = sqrt((cSqVec * k^2 + 4 * sig1Vec * k + ...
       sqrt((cSqVec * k^2 + 4 * sig1Vec * k).^2 + 16 * kappaSqVec * k^2))/2);

% (Fractional) number of intervals
Nfrac = Lvec(1) / hVec(1);
N = floor(Nfrac);
NPrev = N;

%% Initialise Mv and Mw (number of intervals between grid points for left and right system respectively)
Mv = N - numFromRightBound;
Mw = numFromRightBound;   

%% Initialise state vectors (excluding outer boundaries due to Dirichlet conditions)
uNext = zeros (N, 1);
u = zeros (N, 1);
uPrev = u;

%% Initialise matrices approximating spatial derivatives 
% (without the interpolation introduced by the dynamic grid)

% Dxx
DxxNoInterpolation = zeros(Mv + Mw);

ev = ones(Mv, 1);
ew = ones(Mw, 1);
DxxNoInterpolation(1:Mv, 1:Mv) = spdiags([ev -2*ev ev], -1:1, Mv, Mv);         
DxxNoInterpolation(Mv+1:end, Mv+1:end) = spdiags([ew -2*ew ew], -1:1, Mw, Mw);

% Dxxxx
DxxxxNoInterpolation = DxxNoInterpolation * DxxNoInterpolation;

%% Initialise output vector
out = zeros(floor(lengthSound), 1);


%% Main loop
for n = 1:lengthSound  

    % Excite system every half second with a raised cosing scaled by the 
    % current grid spacing
    if mod(n, floor(chunkSize / 2)) == 1
        u(1:5) = u(1:5) + 1/hVec(n) * hann(5);
    end
    
    % Retrieve new fractional number of grid points
    Nfrac = Lvec(n) / hVec(n);
    N = floor(Nfrac);
    
    % Calculate alpha: the fractional part of Nfrac (Eq. (19))
    alpha = Nfrac - N;
    
    % Can only add/remove one point at a time
    if abs(N - NPrev) > 1
        error('Can only add/remove one grid point at a time')
    end
    
    %% Check whether to add or remove points
    if N ~= NPrev
        % Note that uNext will not be altered, as it will be overwritten by
        % the update equation
        if N > NPrev
            %% Add point if N^n > N^{n-1}
        
            % Create cubic interpolator (Eq. (25))
            I3 = [alpha * (alpha + 1) / -((alpha + 2) * (alpha + 3)); ...
                        2 * alpha / (alpha + 2); ...
                        2 / (alpha + 2); ...
                        2 * alpha / -((alpha + 3) * (alpha + 2))]';
             
            v = u(1:Mv);
            w = u(Mv+1:end);
            
            vPrev = uPrev(1:Mv);
            wPrev = uPrev(Mv+1:end);
            
            % Add grid points to v (Eq. (24))
            if numFromRightBound == 1
                v = [v; I3(1:3) * [v(Mv-1:Mv); w(1)]];
                vPrev = [vPrev; I3(1:3) * [vPrev(Mv-1:Mv); wPrev(1)]];
            else
                v = [v; I3 * [v(Mv-1:Mv); w(1:2)]];
                vPrev = [vPrev; I3 * [vPrev(Mv-1:Mv); wPrev(1:2)]];
            end
            
            % Concatenate v and w to yield full state vector u (Eq. (21))
            u = [v; w];
            uPrev = [vPrev; wPrev];
            
            disp("Added Point")

        else
            %% Remove point if N^n < N^{n-1}

            % Remove grid point from system (Eq. (26))
            u(Mv) = [];
            uPrev(Mv) = [];

            disp("Removed Point")
            
        end
        
        %% Refresh Mv and Mw
        Mv = N-numFromRightBound;
        Mw = numFromRightBound;
        
        %% Refresh matrices
        
        % Dxx
        DxxNoInterpolation = zeros(Mv + Mw);
    
        ev = ones(Mv, 1);
        ew = ones(Mw, 1);
        DxxNoInterpolation(1:Mv, 1:Mv) = spdiags([ev -2*ev ev], -1:1, Mv, Mv);
        DxxNoInterpolation(Mv+1:end, Mv+1:end) = spdiags([ew -2*ew ew], -1:1, Mw, Mw);
        
        % Dxxxx
        DxxxxNoInterpolation = DxxNoInterpolation * DxxNoInterpolation;
        
    end

    %% Reinitialise matrices 
    Dxx = DxxNoInterpolation;
    Dxxxx = DxxxxNoInterpolation;

    % Calculate (quadratic) interpolator (Eqs. (17), (18))
    ip = [-(alpha - 1) / (alpha + 1), 1, (alpha - 1) / (alpha + 1)];
    
    % Add the effect of the interpolation to the matrices (Eq. (23))
    if numFromRightBound == 1
        Dxx(Mv, Mv:(Mv+1)) = Dxx(Mv, Mv:(Mv+1)) + fliplr(ip(2:end));
        Dxx(Mv+1, (Mv-1):(Mv+1)) = Dxx(Mv+1, (Mv-1):(Mv+1)) + ip;
        Dxxxx = Dxx * Dxx;
    else
        Dxx(Mv, Mv:(Mv+2)) = Dxx(Mv, Mv:(Mv+2)) + fliplr(ip);
        Dxx(Mv+1, (Mv-1):(Mv+1)) = Dxx(Mv+1, (Mv-1):(Mv+1)) + ip;
        Dxxxx = Dxx * Dxx;
    end

    % Create matrices for update equation (Eq. (32))
    A = (1 + sig0Vec(n) * k);
    B = 2 * eye(N) + cSqVec(n) * k^2 / hVec(n)^2 * Dxx ...
        - kappaSqVec(n) * k^2 / hVec(n)^4 * Dxxxx ...
        + 2 * sig1Vec(n) * k  / hVec(n)^2 * Dxx;
    C = -(1 - sig0Vec(n) * k) * eye(N) - 2 * sig1Vec(n) * k  / hVec(n)^2 * Dxx;
    
    %% Calculate update equation (Eq. (20)) 
    % (Note that as A is a scalar, we don't need a matrix division)
    uNext = (B * u + C * uPrev) / A;

    % Save output
    out(n) = u(5);

    % Update states    
    uPrev = u;
    u = uNext;

    % Update N^{n-1}
    NPrev = N;
    
end

%% Plot Spectrogram
figure ('Position', [440, 378, 516, 250]);
spectrogram(out, 512, 64, 512, fs, 'yaxis');
set(gca, 'Fontsize', 16, 'TickLabelInterpreter', 'latex', ...
    'Position', [0.095 0.18 0.9 0.81]);

labelX = get(gca, 'xLabel');
labelX.Interpreter = 'latex';
labelY = get(gca, 'yLabel');
labelY.Interpreter = 'latex';

set(gcf, 'color', 'w')

%% Play sound (be sure to use soundsc(), and not sound()!!)
% soundsc (out, fs)