
function F=phase2MC(phase0,A,X,Y,Nmp)


Nratio=10; %SLM macropixel in each direction
Nslm=Nratio*Nmp; %Number  of points in amplitude array


%SLM parameters

rslm=1;
dslm=2*rslm/Nslm;
xslm=linspace(-rslm,rslm,Nslm); %slm grid
yslm=linspace(-rslm,rslm,Nslm);

[Xslm,Yslm]=meshgrid(xslm,yslm);

sigmaSLM=0.1658*rslm/0.29; %SLM has gaussian profile in fourier plane
                           % Here largest NA is 0.29 and its itensity is
                           % that 0.1658 sigma. Change according to SLM
                           % beam profile
       

FFamplitude=exp(-(Xslm.^2+Yslm.^2)/(4*sigmaSLM^2));

Npad=20*Nslm;
farfieldpadded=zeros(Npad,Npad);

%% load vector fiber modes
% or use from argument of the function A

dx=abs(X(1,1)-X(1,2));
dy=dx;

%% generate far field
phase=kron(phase0,ones(Nratio));

farfield=FFamplitude.*exp(1i*phase);

farfieldpadded(Npad/2-Nslm/2+1:Npad/2+Nslm/2,Npad/2-Nslm/2+1:Npad/2+Nslm/2)=farfield;

% FFintensity=abs(farfieldpadded).^2;


%% fourier transform

nearfield= fftshift(fft2(farfieldpadded));


% NFintensity=abs(nearfield).^2;

nearfield1=nearfield(Npad/2-Npad/Nratio+1:Npad/2+Npad/Nratio,Npad/2-Npad/Nratio+1:Npad/2+Npad/Nratio);

% NFphase=angle(nearfield1);

%% keep only alternate elements

[Nxnf,Nynf]=size(nearfield1);

xcheck= mod(1:Nxnf,2)*2-1;
mat = (xcheck'*xcheck) > 0;

nearfield1=nearfield1.*mat;

for j=2:2:Nxnf
    nearfield1(:,j)=nearfield1(:,j-1)+nearfield1(:,j);
    nearfield1(:,j-1)=nearfield1(:,j);
end

%% normalize speckle pattern

a=10e-6; %core radius
xnf=linspace(-a,a,Nxnf);
ynf=linspace(-a,a,Nynf);

[Xnf,Ynf]=meshgrid(xnf,ynf);

B=interp2(Xnf,Ynf,nearfield1,X,Y);

B(isnan(B)) = 0;

IB=trapz(trapz(abs(B).^2))*dx*dy;

IB=abs(IB);
B=B/sqrt(IB);

%% generate mode content

[modes_num,Nx,Ny]=size(A);
Mc=zeros(1,modes_num);

    for ii=1:modes_num
        temp=reshape(A(ii,:,:),[Nx,Ny]);
        Mc(ii)=trapz(trapz(temp.*conj(B)))*dx*dy;
    end
    
Mc=abs(Mc).^2;

F=Mc';

end


