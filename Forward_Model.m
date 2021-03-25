function [M,Grid,SensMat,Weights] = Forward_Model(Grid,Elem,Medium,doPlot,Verbosity)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Idan Steinberg, Ph.D 7/30/18
% Multimodality Molecular Imaging Lab (MMIL)
% Stanford University School of Medicine \ Department of Radiology
% 318 Campus Drive, Clark Center,  Room E150
% Stanford, CA 94305-5427
% Mobile: +1(650)469-6437, Work: +1(650)724-3624, Fax: +1(650)723-8649
%
% This function generates the forward sparse matrix which maps each point of
% the initial pressure on a grid to a detector and time. Uses on Bi-linear interpolation,
% includes the effects of detector's directivity, impulse response, sensativity, finate
% aperature and location in a 3D space w.r.t the imaging plane.
%
% Last update - 03/22/2020
% Improved memory mangment and detection of edge points
% Added Grid.Offset
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs -
%     Grid - A structure which includes the following fields:
%       Grid.BLCorner - [2,1], Bottom-left corner x,y coordinates [m]
%       Grid.Res      - [2,1], Grid x,y resolution [m]
%       Grid.Size     - [2,1], Number of points in x,y directions [#]
%       Grid.dt       - Scalar, Sampling interval [sec]
%       Grid.tBounds  - [2,1], Minimal and maximal time points [#], For default values see Output.
%       Grid.SensCut  - Scalar, Sensitivity cut-off - will remove from the model Matrix nodes where
%                       sensitivity is lower than SensCut [#], Default is not none removed.
%                       sensitivity is lower than SensCut [#], Default is not none removed.
%       Grid.Offset   - Constant offset in RF samples [Samples]. Default is zero offset.
%       Grid.DimOrder - Scalar. 1- Column first order (Matlab format), 2- row first order (C format). 
%                       Default is Matlab format.
%
%     Elem - An array of structures, each includes the following fields per Detector element:
%       Elem(i).Position    - [3,1], Position of Detector element on grid [m]
%       Elem(i).Direction   - [2,1], Detector element's polar and azimuthal rotation angles [rad]
%       Elem(i).Size        - [2,1], Detector element width and height [m], Default is eps (point Detector element)
%       Elem(i).DirParam    - [2,1], Polar (elevational) and azimuthal (lateral) directivity
%                                   directivity parameter [#] Default is Size*Fc/Vs
%       Elem(i).DirCut      - [2,1], Polar and azimuthal directivity cutoff [Rad] Default is pi/2
%       Elem(i).Impulse     - Vector, Impulse response function - it is best practice to keep as short a possible.
%                             Default is the Kronecker delta (i.e. [1]);
%       Elem(i).Sensitivity - Scalar, Relative sensitivity (0-1) [#] Default is 1
%       Elem(i).CentFreq    - Scalar, Detector element's central frequency [Hz]
%       Elem(i).ElementNum  - Scalar, Element number [#] Default is i
%
%     Medium - A structure which includes acoustic properties of the medium
%       Medium.Vs    - Scalar, Speed of sound [m/s], Default is 1540 m/s
%       Medium.Grun  - Scalar, Gruneisen coefficient [#] Default is 1
%       Medium.Alpha - Scalar, Acoustic attenuation coefficient [dB/cm/Hz] Default is 0
%
%     doPlot -  Scalar, logical value. If true, will plot every step. This will make the
%               calculation very slow and should only be used for
%               debugging. Default value is 0.
%
%     Verbosity - Scalar, integer value. 0 - Don't display anything
%                                        1 - Display only errors
%                                        2 - Display errors and warnings
%                                        3 - Display errors, warnings and computation stage
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs -
%
%     M - A sparse matrix [Grid.NtD, Grid.Nxy] which maps each point of the initial pressure ona grid to
%         a Detector element and time
%
%     Grid - An updated grid structure which includes the following additional fields:
%         Grid.TRCorner - [2,1], Top-right corner coordinates [m]
%         Grid.tBounds  - [2,1], Minimal and maximal time points [#], If the user supplied these values,
%                         they will be kept for all elements and will not be updated. Otherwise, the values
%                         used buy the algorithm will be returned.
%         Grid.xAxis    - [Grid.Size(1),1], x-Axis node coordinates [m]
%         Grid.yAxis    - [Grid.Size(2),1], y-Axis node coordinates [m]
%         Grid.xMesh    - [Grid.Size(1),Grid.Size(1)], Meshgrid of node coordinates, x component [m]
%         Grid.yMesh    - [Grid.Size(1),Grid.Size(1)], Meshgrid of node coordinates, y component [m]
%         Grid.tAxis    - [max(Grid.tBounds(2,:))-min(Grid.tBounds(1,:))+1,1], time axis points [sec]
%         Grid.Nxy      - Number of image pixels [#]
%         Grid.Rxy      - Surface area of each pixel [m^2]
%         Grid.NtD      - Number of time points times number of detectors [#]
%         Grid.Nt       - Number of time points
%         Grid.GridIdx  - [Grid.Nxy, 1] - If Grid.SensCut was spesified, Grid.GridIdx will contain logical
%                         values for each pixel on the grid that indicate wheter it is included or not.
%
%     SensMat - A matrix [Nd, Grid.Nxy] with values corresponding to the sensitivity of each
%               pixel on the grid
%
%     Weights - Vector, length Grid.Nxy. If the user asks for the weights, the model matrix M will be
%     column-normalized and weights are the normalization weights.
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Code starts here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Anonymous function definitions

% Distance to range function
D2R = @(x,xMin,xMax) max([xMin(:)-x(:),zeros(size(x(:))),x(:)-xMax(:)],[],2); 

% Distance from p0 to all points pt
D2P = @(p0,pt) sqrt(sum((bsxfun(@minus,p0,pt)).^2));

% distance from p0 along a line with angle ang to the edges
Dint = @(p0,pBL,pTR,ang) [(pBL(1)-p0(1))/cos(ang)*(D2R(p0(2)+(pBL(1)-p0(1))*tan(ang),pBL(1),pTR(1)) == 0);
                          (pBL(2)-p0(2))/sin(ang)*(D2R(p0(1)+(pBL(2)-p0(2))/tan(ang),pBL(2),pTR(2)) == 0);
                          (pTR(1)-p0(1))/cos(ang)*(D2R(p0(2)+(pTR(1)-p0(1))*tan(ang),pBL(1),pTR(1)) == 0);
                          (pTR(2)-p0(2))/sin(ang)*(D2R(p0(1)+(pTR(2)-p0(2))/tan(ang),pBL(2),pTR(2)) == 0)];
    
% Closest point on the grid (defined by pBL and pTR) to p0
CPG = @(p0,pBL,pTR) [min(pTR(1),max(p0(1),pBL(1)));min(pTR(2),max(p0(2),pBL(2)));0]; 

% Wrap to Pi around x0, ignore complex angles
wrap2PiAround = @(x,x0) (imag(x)==0).*mod(real(x)+pi-x0,2*pi)-pi+x0+(imag(x)~=0).*x; 

% Tests if ang is in the range between AngMin and AngMax
AngBetween = @(ang,angMin,angMax) (sin((ang-angMin)/2).*sin((angMax-ang)/2))>0;  

% -1 if p0 is below angMin, 0 if p0 is in the sector between angMin angMax and 1 if p0 is above angMax
inSector = @(p0,angMin,angMax) sign(sin(angMin)*p0(1)-cos(angMin)*p0(2))+sign(sin(angMax)*p0(1)-cos(angMax)*p0(2)); 

%% Input Processing
% profile on
M = [];

% General input processing
if exist('doPlot','var')==0 || isempty(doPlot)         % Check for plotting
    doPlot = 0;                                        % Default - no plotting
end

if exist('Verbosity','var')==0 || isempty(Verbosity)   % Check for verbosity level
    Verbosity = 0;                                     % Default - verbosity level 0
end

if isfield(Grid,'SensCut')==0 || isempty(Grid.SensCut) % Check for sensitivity cut-off
    Grid.SensCut = -eps;                               % Default - no cut-off
end

if isfield(Grid,'DimOrder')==0 || isempty(Grid.DimOrder) % Check for sensitivity cut-off
    Grid.DimOrder = 1;                                   % Default dimension 1 first
end

if isfield(Grid,'Offset')==0 || isempty(Grid.Offset) % Check for time axis offset
    Grid.Offset = 0;                                 % Default Offset is 0
end

if isfield(Grid,'tBounds')==0 || isempty(Grid.tBounds) % Check for time bounderies
    UserBnds = false;                                  % Default - bounderies T.B.D
else
    UserBnds = true;
end

if exist('Medium','var')==0 || isempty(Medium)         % Check for medium
    Medium = [];                                       % Generate medium structure
end

if isfield(Medium,'Alpha')==0 || isempty(Medium.Alpha) % Check for acoustic attenuation
    Medium.Alpha = 0;                                  % Default - no acoustic attenuation
end

if isfield(Medium,'Vs')==0 || isempty(Medium.Vs)       % Check for speed of sound
    Medium.Vs = 1540;                                  % Default speed of sound
end
if isfield(Medium,'Grun')==0 || isempty(Medium.Grun)   % Check for Gruneisen coefficient
    Medium.Grun = 0.15;                                % Default Grunheisen parameter
end

% Per element input processing
NElem = length(Elem);                                                          % Num Elements
for el=1:NElem                                                                 % Per Element input check
    if isfield(Elem(el),'Size')==0 || isempty(Elem(el).Size)                   % Check for Element size
        Elem(el).Size = [eps; eps];                                            % Default Element size
    end
    if isfield(Elem(el),'DirParam')==0 || isempty(Elem(el).DirParam)              % Check for directivity parameter
        Elem(el).DirParam = Elem(el).Size(:)*Elem(el).CentFreq/Medium.Vs.*[1; 1]; % Default directivity parameter
    end
    if isfield(Elem(el),'DirCut')==0 || isempty(Elem(el).DirCut)               % Check for directivity cut-off
        Elem(el).DirCut = [pi/2; pi/2];                                        % Default directivity cut-off
    end
    if isfield(Elem(el),'Sensitivity')==0 || isempty(Elem(el).Sensitivity)     % Check for Sensitivity
        Elem(el).Sensitivity = 1;                                              % Default sensitivity
    end
    if isfield(Elem(el),'Impulse')==0 || isempty(Elem(el).Impulse)             % Check for impulse response
        Elem(el).Sensitivity = 1;                                              % Default impulse response
    end
    if isfield(Elem(el),'DetNum')==0 || isempty(Elem(el).DetNum)               % Check for detector number
        Elem(el).DetNum = el;                                                  % Default detector number
    end
end

%% Initialization
if Verbosity > 0, tic; end
if Verbosity > 2, disp('Starting initialization.'); toc; end

% General initialization
Medium.A_np = -2e-7*log10(exp(1))*Medium.Alpha; % Attenuation in nepers/Meter/Hz
CardinalDir = pi*(1:0.5:2.5)';                  % Cardinal directions on the imaging plane
Filter = [1; 0; -1]/(2*Grid.dt);                % 1st order central dervitative

% Grid initialization
Grid.Nxy      = prod(Grid.Size);                                   % Number of pixels in the image [#]
Grid.Rxy      = prod(Grid.Res);                                    % DeltaX * DeltaY
Grid.TRCorner = Grid.BLCorner+Grid.Res.*(Grid.Size-1);             % Top-right corner coordinates [m]
Grid.xAxis    = Grid.BLCorner(1):Grid.Res(1):Grid.TRCorner(1);     % x-Axis node coordinates [m]
Grid.yAxis    = Grid.BLCorner(2):Grid.Res(2):Grid.TRCorner(2);     % y-Axis node coordinates [m]
Grid.Conv     = [1./Grid.Res, (Grid.Res-Grid.BLCorner)./Grid.Res]; % Conversion factors from value to index
[Grid.xMesh,Grid.yMesh] = meshgrid(Grid.xAxis,Grid.yAxis);         % Mesh grid of node coordinates [m]
if UserBnds
    NtMin = Grid.tBounds(1);
    NtMax = Grid.tBounds(2);
else
    NtMin = 0;
    NtMax = Inf;
end
OfstTime = Grid.Offset*Grid.dt;

% Some important points in 3D
P1    = [Grid.BLCorner(:);0];  P2 = [Grid.TRCorner(:);0];
Pedge = combvec([P1(1),P2(1)],[P1(2),P2(2)]); Pedge(3,:) = 0;

% Per element initialization
Dets = unique([Elem(:).DetNum]); % Detector identifiers
Nd = length(Dets);               % Number of detectors

DirCutSum = 0;
[~, ind] = sort([Elem.DetNum]);
Elem = Elem(ind);
MaxImpLen = 0;

if ~UserBnds
    Grid.tBounds(1,1:Nd) = Inf;
    Grid.tBounds(2,1:Nd) = 0;
end

for el = 1:NElem
    if Elem(el).Sensitivity > 0 % Ignore inactive Elements
        
        Elem(el).DetNum = find(Dets==Elem(el).DetNum);
        
        % Determine the last element for every detector
        if ~exist('CurrDet','var')
            CurrDet = Elem(el).DetNum;
        end
        if Elem(el).DetNum ~= CurrDet
            Elem(PrevInd).LastInDet = true;
            CurrDet = Elem(el).DetNum;
        else
            Elem(el).LastInDet = false;
        end
        PrevInd = el;
        
        Elem(el).Impulse   = conv(Elem(el).Impulse, Filter);                            % Apply derivative filter
        if length(Elem(el).Impulse) > MaxImpLen
            MaxImpLen = length(Elem(el).Impulse);
        end
        Elem(el).Direction = wrapToPi(Elem(el).Direction);                              % Element's rotation angle [rad] should be within [-pi pi] range
        Elem(el).Normal    = [sin(Elem(el).Direction(1))*cos(Elem(el).Direction(2));...
            sin(Elem(el).Direction(1))*sin(Elem(el).Direction(2));...
            cos(Elem(el).Direction(1))];                                                % Element's normal
        DirCutSum          = DirCutSum + mean(Elem(el).DirCut);                         % Use the average directivity cutoff for estimating the number of values to be calculated
        Elem(el).EffSize   = pi*Elem(el).Size*Elem(el).CentFreq/Medium.Vs;
        
        % Now we need to find the minimal and maximal time sample for each element
        z_elem = abs(Elem(el).Position(3));
        if z_elem<eps % Element is on the imaging plane
            isParaxial = AngBetween(pi/2,Elem(el).Direction(1)-Elem(el).DirCut(1),Elem(el).Direction(1)+Elem(el).DirCut(1)); % Test to see if the imaging plane is within the element's elevational directivity
            if isParaxial  %  Paraxial case
                Dmax_Dir = Inf;
                Dmin_Dir = 0;
                Elem(el).Face = 0;
            else
                Elem(el).Sensitivity = 0; % this element won't see a thing
            end
        else % Element is outside the imaging plane
            Elem(el).Face = sign(Elem(el).Normal(3)); % Check whether the element is facing up (1), down (-1) or parallel to the imaging plane (0)
            
            % Minimal and maximal distance on the grid based on elevational cut-off
            Dmin_Dir = max(Elem(el).Face*z_elem/cos(Elem(el).Direction(1)-Elem(el).Face*Elem(el).DirCut(1)),0);
            Dmax_Dir = 1./max(Elem(el).Face*cos(Elem(el).Direction(1)+Elem(el).Face*Elem(el).DirCut(1))/z_elem,0);
        end
        
        % Minimal and maximal time step for this element
        ClosestPt = CPG(Elem(el).Position,P1,P2);
        RelPnt = Elem(el).Position(:)-ClosestPt;
        Dmin_Dis = norm(RelPnt);
        if norm(RelPnt) > eps 
            t = inSector(RelPnt,Elem(el).Direction(2)-Elem(el).DirCut(2),Elem(el).Direction(2)+Elem(el).DirCut(2));
            if t < 0
                D = Dint(Elem(el).Position,P1,P2,Elem(el).Direction(2)-Elem(el).DirCut(2));
                Dmin_Dis = min(D(D>0));
            elseif t > 0
                D = Dint(Elem(el).Position,P1,P2,Elem(el).Direction(2)+Elem(el).DirCut(2));
                Dmin_Dis = min(D(D>0));
            end
            if isempty(Dmin_Dis)
                Elem(el).Sensitivity = 0; % this element won't see a thing
            end
        end
        
        Elem(el).MinSamp = max(ceil(max(Dmin_Dir,Dmin_Dis)/(Medium.Vs*Grid.dt))+1,NtMin);
        Dmax_Dis         = max(D2P(Elem(el).Position(:),Pedge));
        Elem(el).MaxSamp = min(floor(min(Dmax_Dir,Dmax_Dis)/(Medium.Vs*Grid.dt)),NtMax);
        if ~UserBnds
            Grid.tBounds(1,Elem(el).DetNum) = min(Grid.tBounds(1,Elem(el).DetNum),Elem(el).MinSamp);
            Grid.tBounds(2,Elem(el).DetNum) = max(Grid.tBounds(2,Elem(el).DetNum),Elem(el).MaxSamp);
        end
    end
end
Elem(PrevInd).LastInDet = true;

%% Grid time initialization
NtMin = min(Grid.tBounds(1,:));
NtMax = max(Grid.tBounds(2,:));
Grid.tAxis = Grid.dt*(NtMin:NtMax); % t axis [sec]
Grid.NtD = (NtMax-NtMin+1)*Nd; % Product of the number of time samples and detectors

%% Preallocate memory
% Preallocate memory for the sensitivity matrix
if nargout > 2
    SensMat = zeros([Nd,Grid.Nxy]);
end

% Define the detector buffer
% NumValsEst = 80*(2^30)/(8*2); % 140 Gb of buffer for double precision values
NumValsEst   = ceil(10*MaxImpLen*Grid.Nxy*Grid.Rxy/(Medium.Vs*Grid.dt)^2*DirCutSum/pi*NElem/Nd); % This is an empiric estimate for the number of values needed. If we go over that algorithm will still work but will be much slower
DetBuffer    = zeros(NumValsEst,2);
DetBufferInd = 1;

% Matrix buffer -
FullMatSize   = Grid.Nxy*Grid.NtD;
NumMatValsEst = ceil(FullMatSize*0.02); % set to about 2% of the final full matrix size
MatBuffer     = zeros(NumMatValsEst,2);
MatBufferInd  = 1;

%% Main loop - compute the matrix values for each active Element and relevant time point
ConstFactor = Medium.Grun/(4*pi*Medium.Vs*Grid.Rxy);
if Verbosity > 2, disp('Finished Initialization. Main loop started'); toc; end
for el = 1:NElem
    if Elem(el).Sensitivity > 0
        %% Arc structure initialization
        Arc.Origin = Elem(el).Position(1:2);                           % Coordinates of the origin [m]
        Arc.Angles = Elem(el).Direction(2)+Elem(el).DirCut(2)*[-1; 1]; % Minimal and maximal angles [rad]
        
        %% For each time point
        for l = Elem(el).MinSamp:Elem(el).MaxSamp
            
            %% Complete arc initialization
            R          = Medium.Vs*(Grid.tAxis(l)+OfstTime);       % Total radius for the current time point
            if R < 0, continue; end
            Arc.RadSq  = R^2-Elem(el).Position(3)^2;               % Squared radius of arc on the imaging plane [m]
            Arc.Radius = sqrt(Arc.RadSq);                          % Radius of arc
            Arc.xSpan  = Arc.Radius*cos(Arc.Angles)+Arc.Origin(1); % x-Coordinates of the boundary points
            Arc.ySpan  = Arc.Radius*sin(Arc.Angles)+Arc.Origin(2); % y-Coordinates of the boundary points
            
            % Calculate the factor
            Arc.factor = Elem(el).Sensitivity*ConstFactor*exp(Medium.A_np*R*Elem(el).CentFreq)*Arc.Radius/R;
            
            %% Check point 1
            if Arc.RadSq < 0
                if Verbosity > 0, error('Error - R < z'); end
                return
            end
            
            %% Edge points - compute all possible 10 and then remove irrelevant ones
            EdgePoints = zeros(4,10); % Allocate memory
            
            CxMin = (Grid.BLCorner(1)-Arc.Origin(1))/Arc.Radius; % Cosine of the angle with the xMin edge
            CxMax = (Grid.TRCorner(1)-Arc.Origin(1))/Arc.Radius; % Cosine of the angle with the xMax edge
            SyMin = (Grid.BLCorner(2)-Arc.Origin(2))/Arc.Radius; % Sine of the angle with the yMin edge
            SyMax = (Grid.TRCorner(2)-Arc.Origin(2))/Arc.Radius; % Sine of the angle with the yMax edge
            
            EdgePoints(3,2:2:10) = wrap2PiAround([acos(CxMin),asin(SyMin),acos(CxMax),...
                asin(SyMax),Arc.Angles(2)],Elem(el).Direction(2));                        % Angles
            EdgePoints(3,1:2:9)  = wrap2PiAround([Arc.Angles(1),-EdgePoints(3,2),pi-EdgePoints(3,4),...
                -EdgePoints(3,6),pi-EdgePoints(3,8)],Elem(el).Direction(2));               % Conjugate Angles
            EdgePoints(1,:) = [Arc.xSpan(1),...
                Grid.BLCorner(1),Grid.BLCorner(1),...                                      % x-coordinate
                Arc.Radius*cos(EdgePoints(3,4))+Arc.Origin(1),...
                Arc.Radius*cos(EdgePoints(3,5))+Arc.Origin(1),...
                Grid.TRCorner(1),Grid.TRCorner(1),...
                Arc.Radius*cos(EdgePoints(3,8))+Arc.Origin(1),...
                Arc.Radius*cos(EdgePoints(3,9))+Arc.Origin(1),...
                Arc.xSpan(2)];
            EdgePoints(2,:) = [Arc.ySpan(1),...
                Arc.Radius*sin(EdgePoints(3,2))+Arc.Origin(2),...                           % y-coordinate
                Arc.Radius*sin(EdgePoints(3,3))+Arc.Origin(2),...
                Grid.BLCorner(2),Grid.BLCorner(2),...
                Arc.Radius*sin(EdgePoints(3,6))+Arc.Origin(2),...
                Arc.Radius*sin(EdgePoints(3,7))+Arc.Origin(2),...
                Grid.TRCorner(2),Grid.TRCorner(2),...
                Arc.ySpan(2)];
            EdgePoints(4,:) = [3,1,1,2,2,1,1,2,2,3];
            
            ind = logical((EdgePoints(1,:)>=Grid.BLCorner(1)).*(EdgePoints(1,:)<=Grid.TRCorner(1)).*(imag(EdgePoints(1,:))==0).*...
                (EdgePoints(2,:)>=Grid.BLCorner(2)).*(EdgePoints(2,:)<=Grid.TRCorner(2)).*(imag(EdgePoints(2,:))==0).*...
                [1, AngBetween(EdgePoints(3,2:9),Arc.Angles(1),Arc.Angles(2)),1]);
            EdgePoints = sortrows(EdgePoints(:,ind)',3)'; % Remove unphysical points and sort the rest
            
            % Remove zero length sections
            if EdgePoints(3,1)     == EdgePoints(3,2),   EdgePoints(:,1) = [];   end
            if EdgePoints(3,end-1) == EdgePoints(3,end), EdgePoints(:,end) = []; end
            
            NumSections = size(EdgePoints,2)/2; % Number of arc section is half the number of edge points
            if NumSections < 1
                if Verbosity > 1,  warning(['Warning, zero arc sections on time step #' num2str(l-Elem(el).MinSamp+1)]); toc; end
                break % if the arc is completely outside the grid - no more time points for this Element
            end
            
            %% Check point 2
            if floor(NumSections)<NumSections   % Sanity check
                if Verbosity > 0,  error('Error - Number of arc sections is not integer'); end
                return
            end
            
            %% Divide to sections and initialize ArcSegments
            InRange = zeros(NumSections,4);
            
            for n=1:NumSections
                InRange(n,:) = AngBetween(CardinalDir,EdgePoints(3,2*n-1),EdgePoints(3,2*n)); % Find out which of the cardinal directions is in range
            end
            NumSegments = ones(NumSections,1); % Number of arc segments for each section
            
            %% Now get all the points on each arc segment
            
            for n=1:NumSections
                
                %% Find the bounding box of the arc section
                if InRange(n,1)
                    BLCorner(1) = max(Arc.Origin(1)-Arc.Radius,Grid.BLCorner(1)); % x min = xc-R or x min
                else
                    BLCorner(1) = min(EdgePoints(1,2*n-1),EdgePoints(1,2*n))+eps; % x min is determined by the edge points
                end
                if InRange(n,2)
                    BLCorner(2) = max(Arc.Origin(2)-Arc.Radius,Grid.BLCorner(2)); % y min = yc-R or y min
                else
                    BLCorner(2) = min(EdgePoints(2,2*n-1),EdgePoints(2,2*n))+eps; % y min is determined by the edge points
                end
                if InRange(n,3)
                    TRCorner(1) = min(Arc.Origin(1)+Arc.Radius,Grid.TRCorner(1)); % x max = xc+R or x max
                else
                    TRCorner(1) = max(EdgePoints(1,2*n-1),EdgePoints(1,2*n))-eps; % x max is determined by the edge points
                end
                if InRange(n,4)
                    TRCorner(2) = min(Arc.Origin(2)+Arc.Radius,Grid.TRCorner(2)); % y max = yc+R or y max
                else
                    TRCorner(2) = max(EdgePoints(2,2*n-1),EdgePoints(2,2*n))-eps; % y max is determined by the edge points
                end
                
                NumPntsEst = 2*(1+ceil((TRCorner(1)-BLCorner(1))/Grid.Res(1))+ceil((TRCorner(2)-BLCorner(2))/Grid.Res(2)));
                Pnts = zeros(6,NumPntsEst);
                
                %% Add the edge points
                Pnts(1:4,1)          = EdgePoints(:,2*n-1);
                Pnts(1:4,NumPntsEst) = EdgePoints(:,2*n);
                
                %% Get all the vertical points
                indx = logical((Grid.xAxis>BLCorner(1)).*(Grid.xAxis<TRCorner(1)));                % Remove x values outside bounding box
                Nv = sum(indx); xc = zeros(1,2*Nv); yc = zeros(1,2*Nv);                            % Preallocate memory
                xc(1:Nv) = Grid.xAxis(indx)-Arc.Origin(1); yc(1:Nv) = sqrt(Arc.RadSq-xc(1:Nv).^2); % Convert to physical values around the center
                xc(Nv+1:2*Nv) = xc(1:Nv);                  yc(Nv+1:2*Nv) = -yc(1:Nv);              % Double side the values to cover all options
                indy = logical((yc>=BLCorner(2)-Arc.Origin(2)).*(yc<=TRCorner(2)-Arc.Origin(2)));    % Remove y values outside bounding box
                Nv = sum(indy);
                Pnts(:,2:Nv+1) = [xc(indy)+Arc.Origin(1); yc(indy)+Arc.Origin(2);...
                    atan2(yc(indy),xc(indy)); ones(1,Nv); zeros(2,Nv)];              % Generate the points structure
                
                %% Get all the horizontal points
                indy = logical((Grid.yAxis>BLCorner(2)).*(Grid.yAxis<TRCorner(2)));                % Remove y values outside bounding box
                Nh = sum(indy); yc = zeros(1,2*Nh); xc = zeros(1,2*Nh);                            % Preallocate memory
                yc(1:Nh) = Grid.yAxis(indy)-Arc.Origin(2); xc(1:Nh) = sqrt(Arc.RadSq-yc(1:Nh).^2); % Convert to physical values around the center
                yc(Nh+1:2*Nh) = yc(1:Nh);                  xc(Nh+1:2*Nh) = -xc(1:Nh);              % Double side the values to cover all options
                indx = logical((xc>=BLCorner(1)-Arc.Origin(1)).*(xc<=TRCorner(1)-Arc.Origin(1)));    % Remove x values outside bounding box
                Nh = sum(indx);
                Pnts(:,1+Nv+(1:Nh)) = [xc(indx)+Arc.Origin(1);yc(indx)+Arc.Origin(2);...
                    atan2(yc(indx),xc(indx));2*ones(1,Nh);zeros(2,Nh)];         % Generate the points structure
                
                %% Remove unused points and sort
                Pnts(:,Nv+Nh+2:NumPntsEst-1) = [];                                  % Remove
                Pnts(3,:) = wrapToPi(Pnts(3,:)-Elem(el).Direction(2));              % Angle w.r.t the element's normal
                NumSegments(n)            = Nv+Nh+1;                                % Get the number of segments (number of points is NumSegments + 1
                Pnts(:,2:NumSegments(n))  = sortrows(Pnts(:,2:NumSegments(n))',3)'; % Sort all points according to angle. Edge points don't need to be sorted
                
                %% Finally, compute the matrix values for each segments
                ind1 = 1:NumSegments(n); ind2 = ind1+1;
                
                % Get the coordinates
                Theta1 = Pnts(3,ind1); Theta2 = Pnts(3,ind2);
                MTheta = (Theta1+Theta2)/2;
                Pnts(5,ind1) = floor((Arc.Origin(1)+Arc.Radius*cos(MTheta+Elem(el).Direction(2)))*Grid.Conv(1,1)+Grid.Conv(1,2)); % x indices
                Pnts(6,ind1) = floor((Arc.Origin(2)+Arc.Radius*sin(MTheta+Elem(el).Direction(2)))*Grid.Conv(2,1)+Grid.Conv(2,2)); % y indices
                
                xBL = Grid.xAxis(Pnts(5,ind1));      yBL = Grid.yAxis(Pnts(6,ind1));
                xTR = xBL+Grid.Res(1);               yTR = yBL+Grid.Res(2);
                Mxn = (Pnts(1,ind2)+Pnts(1,ind1))/2; Myn = (Pnts(2,ind2)+Pnts(2,ind1))/2;
                
                % Differences of all the relevant coordinates
                Dx2x1  = Pnts(1,ind2)-Pnts(1,ind1); Dy2y1  = Pnts(2,ind2)-Pnts(2,ind1);
                DxTRxd = xTR-Arc.Origin(1);         DyTRyd = yTR-Arc.Origin(2);
                DxBLxd = xBL-Arc.Origin(1);         DyBLyd = yBL-Arc.Origin(2);
                DxBLMx = xBL-Mxn;                   DyBLMy = yBL-Myn;
                DxTRMx = xTR-Mxn;                   DyTRMy = yTR-Myn;
                DTheta = wrapTo2Pi(Theta2-Theta1);
                
                % Matrix indices
                indxy = zeros(4,NumSegments(n));
                indt  = l+(NtMax-NtMin+1)*(Elem(el).DetNum-1);
                
                if Grid.DimOrder == 1 % X -first dim, Y - second dim
                    indxy(1,:) = Pnts(5,ind1)+Grid.Size(1)*(Pnts(6,ind1)-1);
                    indxy(2,:) = indxy(1,:) + Grid.Size(1);                
                    indxy(3,:) = indxy(1,:) + 1;
                    indxy(4,:) = indxy(2,:) + 1;
                else %  Y - first Dim, X - Second Dim
                    indxy(1,:) = Pnts(6,ind1) + Grid.Size(2)*(Pnts(5,ind1)-1);
                    indxy(2,:) = indxy(1,:) + 1;
                    indxy(3,:) = indxy(1,:) + Grid.Size(2);
                    indxy(4,:) = indxy(3,:) + 1;
                end
                
                % Calculate the directivity function
                ind = abs(MTheta)>1e-10;
                Arg = sin(MTheta(ind)*Elem(el).DirParam(2))*Elem(el).EffSize(1);
                DirAz = ones(1,NumSegments(n));
                DirAz(ind) = abs(cos(MTheta(ind)*Elem(el).DirParam(2)).*sin(Arg)./Arg);
                switch Elem(el).Face
                    case 1
                        Phi = Elem(el).Direction(1)-pi+atan2(Arc.Radius*cos(MTheta),Elem(el).Position(3));
                        if Phi == 0
                            DirEl = 1;
                        else
                            Arg   = sin(Phi*Elem(el).DirParam(1))*Elem(el).EffSize(2);
                            DirEl = abs(cos(Phi*Elem(el).DirParam(1)).*sin(Arg)./Arg);
                        end
                    case 0
                        DirEl = 1;
                    case -1
                        Phi = Elem(el).Direction(1)-atan2(Arc.Radius*cos(MTheta),Elem(el).Position(3));
                        if Phi == 0
                            DirEl = 1;
                        else
                            Arg   = sin(Phi*Elem(el).DirParam(1))*Elem(el).EffSize(2);
                            DirEl = abs(cos(Phi*Elem(el).DirParam(2)).*sin(Arg)./Arg);
                        end
                end
                
                % Calculate the matrix elements
                Psi = Arc.factor*bsxfun(@times,DirAz.*DirEl,...
                    bsxfun(@times,Dx2x1/2,[DxTRxd+DxTRMx;  -DxTRxd-DxTRMx;  -DxBLxd-DxBLMx;  DxBLxd+DxBLMx])-...
                    bsxfun(@times,Dy2y1/2,[DyTRyd+DyTRMy;  -DyBLyd-DyBLMy;  -DyTRyd-DyTRMy;  DyBLyd+DyBLMy])+...
                    bsxfun(@times,DTheta, [DxTRxd.*DyTRyd; -DxTRxd.*DyBLyd; -DxBLxd.*DyTRyd; DxBLxd.*DyBLyd]));
                
                %% That's it. Now just add the values to the lists of values and subscripts
                NumVals = numel(Psi);
                DetBuffer(DetBufferInd:DetBufferInd+NumVals-1,:)=[Psi(:),indt+Grid.NtD*(indxy(:)-1)];
                DetBufferInd = DetBufferInd+NumVals;
                
                if nargout > 2
                    Sens = reshape(bsxfun(@mtimes,DirEl.*DirAz*Elem(el).Sensitivity*exp(Medium.A_np*R*Elem(el).CentFreq),ones(4,1)),[1, NumSegments(n)*4]);
                    SensMat(Elem(el).DetNum,indxy(:)) = SensMat(Elem(el).DetNum,indxy(:))+Sens.^2;
                end
            end
            
            if doPlot
                plot(Grid.xMesh(:),Grid.yMesh(:),'+')
                title(['Element # ' num2str(el) ', time step # ' num2str(l)]);
                hold all
                axis equal
                Theta = linspace(Arc.Angles(1),Arc.Angles(2),100);
                ArcPnts(1,:) = Arc.Origin(1)+Arc.Radius*cos(Theta);
                ArcPnts(2,:) = Arc.Origin(2)+Arc.Radius*sin(Theta);
                plot(ArcPnts(1,:),ArcPnts(2,:),'k')
                N = [Arc.Origin(:) Arc.Origin(:)+Elem(el).Normal(1:2)*5*Grid.Res(1)];
                plot(N(1,:),N(2,:));
                plot(Arc.Origin(1),Arc.Origin(2),'sm')
                
                ind = Pnts(4,:) == 1;
                plot(Pnts(1,ind),Pnts(2,ind),'or')
                ind = Pnts(4,:) == 2;
                plot(Pnts(1,ind),Pnts(2,ind),'ok')
                ind = Pnts(4,:) == 3;
                plot(Pnts(1,ind),Pnts(2,ind),'*g')
                hold off
                drawnow
            end
        end
    end
    
    %% Accumulate data for the detector, convolve and save to matrix buffer
    if Elem(el).LastInDet
        if Verbosity > 2
            disp(['Detector #' num2str(Elem(el).DetNum) ' finished calculating. Detector buffer ' num2str(round(100*DetBufferInd/NumValsEst)) '% full']); toc;
            disp('Data accumulation and convolution starting');
        end
        
        % Data accumulation
        [IndDet,~,ValsDet] = find(accumarray(DetBuffer(1:DetBufferInd-1,2),DetBuffer(1:DetBufferInd-1,1),[FullMatSize,1],[],0,1));
        
        % Sparse convolution of the impulse response with the detector data
        n = numel(Elem(el).Impulse);
        ind = repmat(0:n-1,numel(IndDet),1)+repmat(IndDet,1,n);
        Vals = (Elem(el).Impulse(:)*ValsDet(:).').';
        [ind,~,Vals] = find(accumarray(ind(:),Vals(:),[FullMatSize,1],[],0,1));
        
        % Remove values that overflow beyond the scope of the current detector
        [IndT,~] = ind2sub([Grid.NtD, Grid.Nxy],ind);
        indOut = IndT>(NtMax*Elem(el).DetNum);
        ind(indOut) = []; Vals(indOut) = [];
        
        % Update the Mat buffer
        NumVals = numel(Vals);
        MatBuffer(MatBufferInd:MatBufferInd+NumVals-1,:)=[Vals(:),ind(:)];
        MatBufferInd = MatBufferInd+NumVals;
        
        % Reset detector buffer and clear variables
        DetBuffer = DetBuffer*0; % zero out the buffer
        DetBufferInd = 1;
        clear ValsDet IndDet Vals ind idxBrep idxArep IndT indOut
        
        if Verbosity > 2, disp(['Convolution finished. Matrix buffer ' num2str(round(100*MatBufferInd/NumMatValsEst)) '% full']); toc; end
    end
end

clear DetBuffer; % free up memeory

%% final step - generate the sparse model matrix using array accumulation
if Verbosity > 2, disp('Using array accumulation to generate the model matrix'); end
sz = [Grid.NtD, Grid.Nxy];
MatBuffer = MatBuffer(1:MatBufferInd-1,:);
[S1,S2] = ind2sub(sz, MatBuffer(:,2));
M = accumarray([S1,S2],MatBuffer(:,1),sz,[],0,true);

%% Post processing if necessary

% Removing pixels with low-sensitivity
if Grid.SensCut > 0
    if Verbosity > 2, disp('Removing pixels with low sensitivity'); toc; end
    W = sqrt(sum(SensMat)');
    Grid.GridIdx = (W./max(W)) >= Grid.SensCut;
    M(:,~Grid.GridIdx) = [];
    Grid.Nxy = sum(Grid.GridIdx);
end

% Matrix normalization
if nargout == 4
    if Verbosity > 2, disp('Normalizing the model matrix columns'); toc; end
    Weights = zeros(Grid.Nxy,1);
    [row,col,Val] = find(M);
    nnz = length(Val);
    clear M;
    IdxStart = 1;
    IdxEnd = IdxStart;
    for k = 1:Grid.Nxy
        C = col(IdxStart);
        while 1
            IdxEnd = IdxEnd +1;
            if IdxEnd > nnz
                IdxEnd = nnz;
                break
            elseif col(IdxEnd) > C
                break
            end
        end
        
        idx = IdxStart:IdxEnd-1;
        V = Val(idx);
        N = norm(V);
        if N > 0
            Val(idx) = V ./ N;
            Weights(C) = N;
        end
        IdxStart = IdxEnd;
    end
    M = sparse(row,col,Val,Grid.NtD,Grid.Nxy);
end

if Verbosity > 2, disp('Forward model done'); toc; end

% profile viewer
end