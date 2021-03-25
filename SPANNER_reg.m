function Xout = SPANNER_reg(M,rf,Xin,tolObj,maxItrProj,maxItrCG,sigma,reg)
 
%% Calculate the 0th iteration
nItrCG = 0;                % Set the CG iteration counter to 0
[m,n]  = size(Xin);        % Get the size of the image
r      = rf - M * Xin(:);  % Calculate the residual 
rho    = r' * r ;          % Calculate the objective function (residual squared- L2 norm)
rhoDes = rho * tolObj;     % Set the target objective function value

%% Calculate the fixed step sizes for the FGP algorithm
s = zeros(100,1);
t = zeros(101,1); t(1) = 1;

for k = 1:100
    t(k+1) = 0.5 * (1 + sqrt(1 + 4 * t(k) ^ 2)); % Step parameter
    s(k)   = (t(k) - 1) / t(k+1);                % Step size
end

%% Main loop - Run until convergence or maximal # of iterations
while rho > rhoDes && nItrCG < maxItrCG
    nItrCG = nItrCG + 1; % Update the CG iteration counter
    g      = -(M' * r) + reg * Xin(:);  % Calculate the gradient
    gamma  = g' * g;     % Calculate the gradient squared- L2 norm
    if nItrCG == 1 
        delta = 0;
        beta  = 1;
        d     = -g; % Gradient descent
    else 
        delta   = g' * d;
        epslion = g' * gPrv;
        beta    = (gamma - max(0, epslion * gamma / gammaPrv)) / max(deltaAux + delta, gammaPrv); % Jian, Han and Jiang hybrid coefficient
        d       = -g / beta + d; % Update the direction
    end
    q        = beta * (M * d);                  % Calculate the conjugate direction
    theta    = q' * q + reg * beta^2 * (d' * d);  % Conjugate direction squared- L2 norm
    deltaAux = gamma - beta * delta;            % Update the delta parameter
    alpha    = beta * deltaAux / theta;         % CG Step size
    Xin      = Xin + reshape(alpha * d, [m,n]); % Update the (unconstrained) estimate
    
    %% Fast Gradient Projection loop
    sigma  = sigma*0.975; % Superiorization parameter
    lambda = 8*sigma;     % The inverse of an upper bound on the Lipschitz constant
    DxPrv = zeros(m,n); DyPrv = DxPrv;
    Xout   = max(Xin, 0); % Project the unconstaruiend estaimte into the feasble domain

    for k =1:maxItrProj % take nIterProj steps

        Dx = DxPrv * lambda + conv2(Xout,[-1 1]','same'); % Taking a step towards minus of the gradient (x-direction)
        Dy = DyPrv * lambda + conv2(Xout,[-1 1] ,'same'); % Taking a step towards minus of the gradient (y-direction)
         
        % Anisotropic (l1-norm) total-variation projection of the gradient
        Dx = Dx ./ max(lambda, abs(Dx));
        Dy = Dy ./ max(lambda, abs(Dy));
        
        DxPrv = Dx + s(k) * (Dx - DxPrv); %
        DyPrv = Dy + s(k) * (Dy - DyPrv); %
      
        Lx = conv2(DxPrv,[0 1 -1]','same');
        Ly = conv2(DyPrv,[0 1 -1] ,'same');

%         Xout = max(Xin - sigma * (Lx+Ly),0);                   % Projected estimate
        Xout = max(Xin - sigma * (Lx+Ly),0);                   % Projected estimate

    end
    
    Xin      = Xout;            % Save the previous estimate 
    gPrv     = g;               % Save the previous gradient
    gammaPrv = gamma;           % Save the previous gradient squared- L2 norm 
    r        = rf - M * Xin(:); % Calculate the residual 
    rho      = r' * r;          % Calculate the objective function
end
