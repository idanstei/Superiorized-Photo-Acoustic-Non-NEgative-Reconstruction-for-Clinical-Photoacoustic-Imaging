function RF_Deconv = DeconvRF(RF,IRF,Nt,Param)
%% Removes the effcts of the TGC on the RF signal

RF    = RF(1:Nt,:);                                           % Use only the requierd time window
% RFSum = sum(RF,2);
Nz = 0;%length(find(RFSum ==0));
L_IRF = length(IRF);
Nd =  size(RF,2);

%% for lsqr deconv
M_IRF = convmtx(IRF(:),Nt-Nz);
M_IRF = M_IRF(ceil(L_IRF/2):end-floor(L_IRF/2),:);
% DerivativeFilter = [13 -190 1305 -4680 9690 -12276 9690 -4680 1305 -190 13]/240;
% M_reg_LF = convmtx(gausswin(16),Nt-Nz);
% M_reg_LF = M_reg_LF (8:end-8,:);
% M_reg_HF = convmtx(DerivativeFilter',Nt-Nz);
% M_reg_HF = M_reg_HF (5:end-6,:);
Weights = full(sqrt(sum(M_IRF.^2,1)));
M_norm = bsxfun(@rdivide,M_IRF,Weights);
M_reg = Param*eye(Nt-Nz); %+M_reg_HF + M_reg_LF;

M = [M_norm; M_reg];

% For L1 sparse deconv
% M = Nt + L_IRF - 1;
% H = sparse(M,Nt);
% e = ones(Nt,1);
% for k = 0:L_IRF-1
%     H = H + spdiags(IRF(k+1)*e, -k, M, Nt);            % H : convolution matrix (sparse)
% end
% 
% phi_L1 = @(x) lam * abs(x);
% wfun_L1 = @(x) abs(x)/lam;
% dphi_L1 = @(x) lam * sign(x);

%% For Weiner deconv
% Sx = fft(RF(Nz+1:end,:));
% IRF = circshift([IRF(:); zeros(size(Sx,1)-L_IRF,1)],-0*floor(L_IRF/2));
% H = fft(IRF);
% Sn = Param;

%% Apply inverse for each ch
RF_Deconv = zeros(size(RF));
% M_InvReg = pinv(M_IRF'*M_IRF + Param^2*eye(Nt-Nz))*M_IRF'; % Calculate the regualrized inverse
% M_InvReg = TGSVDinv(M_IRF,Param);

for d = 1:Nd
%     RF_Deconv(:,d) = M_InvReg*RF(:,d);
%     RF_Deconv(Nz+1:Nt,d) = ifft(conj(H).*Sx(:,d)./max(eps,abs(H).^2.*Sx(:,d)+1),'symmetric');
    [x,~] = lsqr(M,[RF(Nz+1:Nt,d); zeros(Nt-Nz,1)],1e-6,1e5,[],[],zeros(Nt-Nz,1));
    RF_Deconv(Nz+1:Nt,d) = x./Weights(:); 
    RF_Deconv(:,d) = RF_Deconv(:,d)  - mean(RF_Deconv(:,d) ); % Remove the DC
    RF_Deconv(:,d) = RF_Deconv(:,d)*norm(RF(:,d))/norm(RF_Deconv(:,d)); % preserve norm
end

end


% [x1, cost1] = deconv_MM(y, phi_L1, wfun_L1, H, 50);
