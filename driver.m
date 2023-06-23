% Two species model example

clear all
clc

  close all
  %% Tissue Parameters
  T1pmean = [ 30 ]; % s
  T1pstdd = [ 10 ]; % s
  T1plb   = [ 5  ]; % s
  T1pub   = [ 45 ]; % s
  T1lmean = [ 25 ]; % s
  T1lstdd = [ 10 ]; % s
  T1llb   = [  5 ]; % s
  T1lub   = [ 45 ]; % s
  kplmean = [ .15 ];       % s
  kplstdd = [ .03 ];       % s
  kpllb   = [ .01 ];       % s
  kplub   = [ .35 ];       % s
  kvemean = [ 0.05 ];      % s
  kvestdd = [ .01  ];      % s
  kvelb   = [ 0.01 ];      % s
  kveub   = [ 0.20 ];      % s
  t0mean  = [ 4    ];      % s
  t0sttd  = [ 1.3  ];      % s
  t0lb    = [ 0    ];       % s
  t0ub    = [ 7    ];       % s
  alphamean  =  [2.5];
  alphasttd  =  [.3];
  betamean  =  [4.5];
  betasttd  =  [.3];
  tisinput=[T1pmean; T1pstdd; T1lmean; T1lstdd; kplmean; kplstdd; kvemean; kvestdd;t0mean;t0sttd;alphamean; alphasttd; betamean ; betasttd ];
  tisinputlbub=[T1plb; T1pub; T1llb;T1lub; kpllb; kplub; kvelb; kveub;t0lb;t0ub;alphamean-2*alphasttd;alphamean+2*alphasttd; betamean-2*betasttd ;betamean+2*betasttd ];
  
  %% Variable Setup
  Ntime = 30;
  currentTR = 3;
  TR_list = (0:(Ntime-1))*currentTR ;
  M0 = [0,0];
  ve = 0.95;
  %ve = 1.;
  VIF_scale_fact = [100;0];
  bb_flip_angle = 20;
  opts = optimset('lsqcurvefit');
  opts.TolFun = 1e-09;
  opts.TolX = 1e-09;
  opts.Display = 'off';
  params = struct('t0',[t0mean(1);0],'gammaPdfA',[alphamean(1)  ;1],'gammaPdfB',[betamean(1);1],...
      'scaleFactor',VIF_scale_fact,'T1s',[T1pmean(1),T1lmean(1)],'ExchangeTerms',[0,kplmean(1) ;0,0],...
      'TRList',TR_list,'PerfusionTerms',[kvemean(1),0],'volumeFractions',ve,...
      'fitOptions', opts)
  model = HPKinetics.NewMultiPoolTofftsGammaVIF_Edit2();
  
  
  %% Get true Mz
  %% Choose Excitation Angle
  FAType = {'Const'};
  %% HACK- @cmwalker code for initial conditions - https://github.com/fuentesdt/TumorHPMRI/blob/master/models/gPC/walker/ShowMxyPub.m
  for i = 1:numel(FAType)
      switch (FAType{i})
          case('Const') % Nagashima for lactate const 10 pyruvate
              tic
              E1(1) = exp(-currentTR *(1/T1pmean+kplmean));
              E1(2) = exp(-currentTR /T1lmean);
              for n = 1:Ntime
                  % 20deg for pyruvate 30deg for lactate - currently used in brain
                  flips(2,n) = 30*pi/180;
                  flips(1,n) = 20*pi/180;
              end
              params.FaList = flips ;
      end
  
      
      
  
      tic
      %% Fitting
      [t_axis,Mxy,Mz] = model.compile(M0.',params);
      toc
  end
  
  %% Plot initial guess
  plotinit = true;
  if plotinit
      % plot initial guess
      figure(1)
      plot(TR_list,Mxy(1,:),'b',TR_list,Mxy(2,:),'k')
      ylabel('Const Mxy')
      xlabel('sec')
      figure(2)
      plot(TR_list,params.FaList(1,:)*180/pi,'b',TR_list,params.FaList(2,:)*180/pi,'k')
      ylabel('Const FA (deg) ')
      xlabel('sec')
      figure(3)
      plot(TR_list,Mz(1,:),'b--',TR_list,Mz(2,:),'k--')
      hold
      plot(TR_list,Mz(1,:)./cos(params.FaList(1,:)),'b',TR_list,Mz(2,:)./cos(params.FaList(2,:)),'k')
      ylabel('Const Mz')
      xlabel('sec')
  
      % plot gamma
      jmA0    = VIF_scale_fact(1);
      jmalpha = alphamean(1);
      jmbeta  = betamean(1);
      jmt0    = t0mean(1);
      jmaif   = jmA0  * gampdf(TR_list - jmt0  , jmalpha , jmbeta);
      figure(4)
      plot(TR_list,jmaif ,'b')
      ylabel('aif')
      xlabel('sec')
  end
  
  
