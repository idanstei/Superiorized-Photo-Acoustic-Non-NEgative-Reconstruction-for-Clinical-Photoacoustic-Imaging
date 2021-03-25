function PA_DAS_Recon = DAS(RF_data,Grid,Det,Prop,DoFilt)

%% Time and Freq axis
[Grid.Nt,Nd] = size(RF_data);   % Number of Detectors
fs = 1/Grid.dt;
Nxy = prod(Grid.Size);

%% Calculate Delay Matrix
Delay = zeros(Nxy,Nd);
for d=1:Nd
    Delay(:,d) = fs/Prop.Vs*reshape(hypot(Grid.xMesh'-Det(d).Position(1),Grid.yMesh'-Det(d).Position(2)),Nxy,1);
end

%% Filtering
if DoFilt
    [M,Wn,beta,typ] = kaiserord(1e6*[1.5 2 5 6.5],[0 1 0],[10^(-60/20) 10^(0.1/20)-1 10^(-60/20)],fs); % get the filter order
    b = fir1(M,Wn,typ,kaiser(M+1,beta),'noscale'); % design filter
    [gd,~] = grpdelay(b,1); % Get the filter's group delay
    RF_data = filter(b,1,RF_data); %Apply filter
    Delay = Delay + mean(gd); % Compansate for the average group delay
end

%% Apodization
Apod(1,:) = ones(Nd,1); % chebwin(Num_Rx,60) / ones(Num_Rx,1) / gausswin(Num_Rx) / hann(Num_Rx) / tukeywin(Num_Rx);
RF_data = RF_data.*repmat(Apod(1,:),[Grid.Nt 1]); %Applay apodization
RF_data = hilbert(RF_data);% Hilbert Transform

%% Beamforming
PA_DAS_Recon = zeros(size(Grid.xMesh(:)));
t = 0.5:1:Grid.Nt-0.5;        % t axis in samples
for d = 1:Nd
        RF_Ch = RF_data(:,d); % Only pick the RF data for the current cahnnel
        PA_DAS_Recon = PA_DAS_Recon+interp1(t,RF_Ch,Delay(:,d),'linear',0); % Delay and Sum
end
PA_DAS_Recon = abs(reshape(PA_DAS_Recon,Grid.Size'))/Nd; % Envelope detection