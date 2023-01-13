function [para,paramax,paramin,betamap,alphamap]=initializepara_eakf(dailyincidence,num_ens,parafit)
%Z,D,mu,theta,alpha1,alpha2,...,alpha3142,beta1,...,beta3142
Zlow=2;Zup=5;%latency period
Dlow=2;Dup=5;%infectious period
mulow=0.2;muup=1;%relative transmissibility
thetalow=0.01;thetaup=0.3;%movement factor
alphalow=0.05;alphaup=1;%reporting rate
betalow=0.2;betaup=2.5;%transmission rate

num_loc=size(dailyincidence,1);

%define alpha
alphamap=4+(1:num_loc)';

%define beta
betamap=4+num_loc+(1:num_loc)';

load popd
PD=log10(popd);
PD_median=median(PD);

%Z,D,mu,theta,alpha1,alpha2,...,alpha3142,beta1,...,beta3142
paramin=[Zlow;Dlow;mulow;thetalow;ones(num_loc,1)*alphalow;ones(num_loc,1)*betalow];
paramax=[Zup;Dup;muup;thetaup;ones(num_loc,1)*alphaup;ones(num_loc,1)*betaup];

para=zeros(size(paramin,1),num_ens);

%parafit:beta;mu;Z;D;alpha;theta
%Z
para(1,:)=median(parafit(3,:))*ones(1,num_ens);
%D
para(2,:)=median(parafit(4,:))*ones(1,num_ens);
%mu
para(3,:)=median(parafit(2,:))*ones(1,num_ens);
%theta
para(4,:)=median(parafit(6,:))*ones(1,num_ens);

for l=1:num_loc
    %alpha
    para(alphamap(l),:)=parafit(5,ceil(rand(1,num_ens)*size(parafit,2)));
    %beta
    scale=0.8;
    factor=max(0.1,min(2,PD(l)/PD_median*scale));
    para(betamap(l),:)=parafit(1,ceil(rand(1,num_ens)*size(parafit,2)))*factor;
end

