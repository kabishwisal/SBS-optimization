%% %% random phase (or phase loop)

tic
round_number=1; %start with 1, set i+1 after running round=i

load('normfibermodes.mat','A','X','Y') %run normfibermodes.m to modify file
load('GS-normalized-MDL.mat','spectrum1'); %load gain_spectrum for threshold calc


%%
MPnumber = 13 ; MPtot = MPnumber^2; %number of macropixel in the SLM
   StepNum = 1 ;  clear mnSetRand

modes_num=size(A,1);   
No=modes_num; %Number of optical modes
McStore=zeros(No,MPtot); %for storing mode content

   % If using mnSet, MP progresses from first column, second to...
   for m = 1:MPnumber
       for n = 1:MPnumber
           mnSet(StepNum,:) = [m,n] ;
           StepNum = StepNum + 1 ;
       end
   end 

  % mnSetRand is an mnSet with order shuffled
  StepNumSet = randperm(MPtot) ; 
   for StepNum = 1:MPtot 
       mnSetRand(StepNumSet(StepNum),:) = mnSet(StepNum,:) ;  
   end
  mnSet = mnSetRand ; 

 %round of optimizaation
 
 if round_number==1
  Phase0 = rand(MPnumber) ;% Starting phase pattern; or zeros(MPnumber) round 1
 else
  Phase0=PhaseOpt; %round i>1, store PhaseOpt from round i-1
 end


  Phase=Phase0;
  ScanRes = 0.1 ; ScanRange = 0:ScanRes:(1-ScanRes) ;

  % Optimization
  clear ThresholdScan
  MaxthMP=zeros(MPtot,1);
  ThresholdScanFull=zeros(length(ScanRange),MPtot);

  for StepNum = 1:MPtot
      
      m = mnSet(StepNum,1) ; n = mnSet(StepNum,2) ;
      
      ThresholdScan=zeros(length(ScanRange),1);
      
      for Scan = 1:length(ScanRange)
          
        disp([StepNum,Scan]) %for tracking
        
        Phase(m,n) =ScanRange(Scan)+Phase0(m,n);
        Ephase = exp(1i*Phase*2*pi) ;
        
        phase0=2*pi*Phase;
        Mc=phase2MC(phase0,A,X,Y,MPnumber); % obtain MC from phase
        
        ThresholdScan(Scan)=SBSth(Mc,spectrum1); %obtain threshold from MC
        
      end
      
      [maxth,pos] = max(ThresholdScan) ;
      ThresholdScanFull(:,StepNum)=ThresholdScan(:);
      MaxthMP(StepNum)=maxth;
      Phase(m,n) = ScanRange(pos)+Phase0(m,n);
      
      phase0=2*pi*Phase;
      Mc=phase2MC(phase0,A,X,Y,MPnumber); % obtain MC from phase
      McStore(:,StepNum)=Mc; %Store Mc for later analysis
  end
  

  %%% Phase is then optimized. 
  
BGS=BGSeff(Mc,spectrum1);
[pks,locs,w]=findpeaks(BGS);

PhaseOpt=Phase;  
phase0=2*pi*Phase; 
Mc=phase2MC(phase0,A,X,Y,MPnumber); % obtain MC from phase, [A,X,Y] are optical fields with coordinates
Opt_th=SBSth(Mc,spectrum1); %obtain threshold from MC

% Mc gives optimal mode content and Opt_th gives optimum SBS threshold

toc