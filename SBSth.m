function F=SBSth(Mc,spectrum1)

[Nw,No,No]=size(spectrum1);


%% Mode Dependent loss

%loss parameters
load('modeidx.mat','Ifull1')
AzimIdx=Ifull1(:,1);
RadIdx=Ifull1(:,2);

a=0; %loss slope azim idx, change to non zero for MDL
b=0; %loss slope rrad idx
n1=2; %order of polynomial loss
n2=1;

losscoeff= 0.0001+a*abs(AzimIdx).^n1+b*abs(RadIdx).^n2; %0.0001 to avoid division by zero
Modeloss=exp(-losscoeff);

%% Threshold calculation  
Mcloss= Mc.*Modeloss;
Tloss=sum(Mcloss);



Aj =(Mc.*(1-exp(-losscoeff))./losscoeff)./Tloss; %with MDL
% Aj=Mc/T; %without MDL


% Calculate the SBS for the launched modes
Pm=Aj;

Gm=zeros(Nw,No);

for i=1:No
    for j=1:No
        Gij=reshape(spectrum1(:,i,j),[Nw,1]);
        Gm(:,i)=Gm(:,i)+Gij*Pm(j);
    end
end

Gmpeak=max(Gm);

TI=max(Gmpeak);

SBSsupp=1./TI;

F=SBSsupp;
end





