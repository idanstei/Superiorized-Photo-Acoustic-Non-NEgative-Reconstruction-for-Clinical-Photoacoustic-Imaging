%% Start fresh
close all
clear
clc

%% Define parameters
% Note (!!!) unlike ultrasound, here x is azimuth, y is depth and z is (out of plane) elevation
% Grid parameters
Grid.Size     = [246;300];                       % grid size is 246 by 300 pixels
Grid.Res      = 120e-6*[1;1];                    % resolution is 120 micron in both axes
Grid.BLCorner = [-Grid.Size(1)/2*Grid.Res(1);0]; % Bottom left corner is at y=0 and x is centered around the middle
Grid.dt       = 1/20e6;                          % 20 MHz sampling rate
Grid.tBounds  = [1; 1024];                       % use sampels #1 to #1024

% Acoustic properties
Medium.Vs    = 1540; % Speed of sound in tissue [m/s]
Medium.Grun  = 0.2;    % Gruneisen Coefficient [#]
Medium.Alpha = 0.5;    % Acoustic attenuation [dB/cm/Hz]

% Transducer elements parameters
Nd = 128; % number of detector elements
IRF = [-0.324892304172815;0.219788627689701;0.121053766285138;-0.129080581501210;-0.00696699183965304;0.0394727775122044;-0.00273608678533026;-0.0120324853855231;-0.0104804444675527;-0.0102585280483687;0.00151384809487678;0.00564257530091734;-0.00339179900930122;0.00341445484986071;0.00318142614718585;-0.00219336587049555;0.000708984054294864;0.000912933366712991;-0.000238776960496660;3.63310685201789e-05;0.000118795059741864;0.000124766004700987;0.000131774584150708;0.000151823433573222;0.000153371782553271;0.000139041642914682;0.000179637804931685;-0.000404356508806117;-0.00103091726483244;-0.00466864718117263;-0.00533985643173869;-0.00147315215891280;0.00549653998149071;0.0456064470093856;0.0840594195749200;0.149210371891805;-0.167687644783456];

% Noise model
A2DRangeUtil = 0.25;
A2DBits = 14;
RF_noise_amp = 0.05;
Img_noise_amp = 0.05;

Normalize_model = true;

%% first we generate an empty template with for all of the detector's elements
Elem = repmat(struct(...
    'Position',    [0,0,0], ... % Position w.r.t grid [m]
    'Direction',   [0; 0], ...  % Polar and azimuthal rotation angles [rad]
    'Size',        [0; 0], ...  % Width and height [m]
    'DirCut',      [0;0], ...   % Polar and azimuthal directivity cutoffs [Rad]
    'Impulse',     0, ...       % Impulse response function [#]
    'Sensitivity', 0, ...       % Relative sensitivity [#]
    'CentFreq',    0,...        % Central frequency [Hz]
    'DetNum', 0)...             % Element number
    ,1,Nd);                     % Replicate for each element

% Now we will fill the structure with content
for k = 1:Nd
    Elem(k).Position(1)    = 200e-6*(-Nd/2-0.5+(k-1)); % linear array with 200 micron element apart
    Elem(k).Position(2:3)  = [-0.5e-3,0];              % all elements are on y=-1/2 mm plane
    Elem(k).Direction(1:2) = pi/2*[1,1];               % all elements are pointing torward positive y (i.e. normal is [0 1 0])
    Elem(k).Size           = [200e-6, 4e-3];           % each element is 200 micro wide and 4 mm tall
    Elem(k).DirCut         = pi/4*[1,1];               % +- 45 degrees cutoff
    Elem(k).Impulse        = 1;                        % Use Ideal Impulse response function [#]
    Elem(k).Sensitivity    = 1;                        % All elements are 100% senstive [#]
    Elem(k).CentFreq       = 5e6;                      % 5 MHz central frequency [Hz]
    Elem(k).DetNum         = k;                        % No devision into sub elements
end

%% View geometry
Grid.TRCorner = Grid.BLCorner+Grid.Res.*(Grid.Size-1);             % Top-right corner coordinates [m]
figure(1),patch('XData',[Grid.BLCorner(1) Grid.TRCorner(1) Grid.TRCorner(1) Grid.BLCorner(1) Grid.BLCorner(1)],...
    'YData',[Grid.BLCorner(2) Grid.BLCorner(2) Grid.TRCorner(2) Grid.TRCorner(2) Grid.BLCorner(2)],...
    'ZData',zeros(5,1),'FaceColor',[0 0 1],'FaceAlpha',0.05);
hold all;
for el = 1:length(Elem)
    Elem(el).Normal    = [sin(Elem(el).Direction(1))*cos(Elem(el).Direction(2));...
        sin(Elem(el).Direction(1))*sin(Elem(el).Direction(2)); cos(Elem(el).Direction(1))]; % Element's normal
    [Xr,Yr,Zr] = Rect3D(Elem(el).Position,Elem(el).Size,Elem(el).Normal);
    patch('XData',Xr,'YData',Yr,'ZData',Zr,'FaceColor',[1 0 0],'FaceAlpha',0.5,'LineStyle','none');
    plot3(Elem(el).Position(1)+[0 ,3e-3*Elem(el).Normal(1)],Elem(el).Position(2)+[0 ,3e-3*Elem(el).Normal(2)],Elem(el).Position(3)+[0 ,3e-3*Elem(el).Normal(3)],'k');
end
axis equal, grid on; xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);

%% Get the forward model matrix - shouldn't take long with the current settings
[M,Grid] = Forward_Model(Grid,Elem,Medium,false,3);

%% As an example genrate a simple image and get the predicted RF
P0 = double(imread('coins.png'));   % Load one of Matlab's built in image
P0 = abs(P0 - median(P0(:))); % remove backgroud
P0 = P0./max(P0(:)); % normalize
P0_dB = 20*log10(P0); 
figure, imagesc(Grid.xAxis,Grid.yAxis,P0'); colormap hot; colorbar; title('Input Image of P0'); axis image; caxis([0 1]); % Show the input image
P0_noisy = P0 + Img_noise_amp*randn(Grid.Size'); % Add image noise

RF = reshape(M*P0_noisy(:),[Grid.NtD/Nd,Nd]); % Calculating the RF is as simple as a matrix-vector multiplication
RF_Filtered = zeros(size(RF));
for ch = 1:Nd
    RF_Filtered(:,ch) = conv(RF(:,ch),IRF,'same'); % Apply IRF
 end
RF_noisy = RF_Filtered + RF_noise_amp*randn(size(RF)); % Add RF noise
RF_noisy = Quantization(RF_noisy,A2DBits,A2DRangeUtil); % Quantize

figure, imagesc(Grid.tAxis,1:Nd,RF_noisy'); colormap gray; colorbar; title('Predicted RF'); % Show the predicted RF

%% Inverse using DAS and SPANNER
RF_Deconv = DeconvRF(RF_noisy,IRF,Grid.tBounds(2),0.1); % pre-processing step - deconvolve RF

% DAS recon
PA_DAS_Recon = DAS(RF_Deconv,Grid,Elem,Medium,0);
PA_DAS_Recon = PA_DAS_Recon./max(PA_DAS_Recon(:)); % normalize
PA_DAS_Recon_dB = 20*log10(PA_DAS_Recon);
figure, imagesc(Grid.xAxis,Grid.yAxis,PA_DAS_Recon'); colormap hot; colorbar; title('Predicted image using DAS'); axis image;  % Show the input image

% Normalize the model matrix to improve results. This step is optional
if Normalize_model
    W = sqrt(sum(M.^2,1));
    W(W==0) = 1;
    for k = 1:length(W)
        M(:,k) = M(:,k)./W(k);
    end
    W = reshape(W,Grid.Size');
    figure, imagesc(Grid.xAxis,Grid.yAxis,(W.^2)');  colormap jet; colorbar; title('Senstivity map'); axis image; % Show the senstivity image
end

% SPANNER recon
Init = zeros(Grid.Size');
PA_SPNR_Recon = SPANNER_reg(M,RF_Deconv(:),Init,1e-6,20,20,0.05,10);
if Normalize_model
	PA_SPNR_Recon = PA_SPNR_Recon./max(W,5);
end
PA_SPNR_Recon = PA_SPNR_Recon./max(PA_SPNR_Recon(:)); % normalize
PA_SPNR_Recon_dB = 20*log10(PA_SPNR_Recon);
figure, imagesc(Grid.xAxis,Grid.yAxis,PA_SPNR_Recon'); colormap hot; colorbar; title('Predicted image using SPANNER'); axis image;  % Show the input image

%% Utility function for drawing 3D rectangles
function [X,Y,Z] = Rect3D(Center,Size,Normal)
 
    if Normal(1) == 0 && Normal(2) ==0
        V1 = [1,0,0];
        V2 = [0,1,0];
    else
        V1 = cross(Normal,[0,0,1]);
        V1 = V1./norm(V1);
        V2 = cross(V1,Normal);
    end
 
    P1 = Center - V1*Size(1)/2 - V2*Size(2)/2;
    P2 = Center + V1*Size(1)/2 - V2*Size(2)/2;
    P3 = Center + V1*Size(1)/2 + V2*Size(2)/2;
    P4 = Center - V1*Size(1)/2 + V2*Size(2)/2;
 
    X = [P1(1), P2(1), P3(1), P4(1), P1(1)]; 
    Y = [P1(2), P2(2), P3(2), P4(2), P1(2)]; 
    Z = [P1(3), P2(3), P3(3), P4(3), P1(3)]; 
end
