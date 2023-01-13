function infer(lastsavefolder )
%% this is the code that runs the model inference.  This is used to fit the baseline scenario
if nargin==0
    lastsavefolder=''
end
 
savevaxstart=1;
if strcmp(lastsavefolder,'')
    startfromlast=0;
else
    startfromlast=1;
    startfromlastvax=1; % SET to 1 if you want to restart from the start of the vaccinations dec 14
    % optional reads in lastsavefolder
end

%save inference result Tsave days before the latest day
Tsave=14;%two weeks to deal with backfill
 
    load geomaps.mat
    load commutedata
    load population
    num_loc=size(part,1)-1;
    num_mp=size(nl,1);
    num_ens=100;%number of ensemble
    load countyfips
    load dailyincidence
    load dailydeaths
    statefips=readtable('statefips.csv');
    countyfips.statefips=str2double(countyfips.STATE);
    countyfips=outerjoin(countyfips,statefips,'LeftKeys','statefips','RightKeys','fips','RightVariables','state_name');
    
    num_times=size(dailyincidence,2);%total length of data
    T=num_times+6*7;%projection for 6 weeks
    
    Td=9;%average reporting delay
    a=1.85;%shape parameter of gamma distribution
    b=Td/a;%scale parameter of gamma distribution
    rnds=ceil(gamrnd(a,b,1e6,1));%pre-generate gamma random numbers
    
    Td_death=15.4+7;%average delay of death
    a=2.6;%shape parameter of gamma distribution
    b=Td_death/a;%scale parameter of gamma distribution
    rnds_d=ceil(gamrnd(a,b,1e6,1));%pre-generate gamma random numbers
    
    obs_case=zeros(size(dailyincidence));
    obs_death=zeros(size(dailydeaths));
    %% vaccine stuff

    sv=readtable('statevaxdata.csv');
    listoflocs=[statefips.state_abbr;{'BP2','DD2','IH2','LTC','VA2','US'}'];
    
    % trim Dec
    inddec=find(sv.Date<='2021-01-11');
    sv(inddec,:)=[]; % just 4 days of no data at US level. also remove through jan 11 everywhere
    
    for i=1:length(listoflocs);
        ind=find(strcmp(sv.Location,listoflocs(i)));
        temp=sv(ind,:);
        vax1(i,:)=temp.Administered_Dose1;
        vax2(i,:)=temp.Administered_Dose1_Recip;
    end
    
    mergevax=vax2; % to start, entire Administre_Dose1_Recip;
    mergevax(:,1:38)=vax1(:,1:38); % replace the first 46 days with Administed_Dose1;
    
    for i=1:length(listoflocs)
        fillin(i,:)=interp1([1 29],[0 mergevax(i,1)],1:29);
    end
    mergevax=[fillin mergevax];
    
    dailyvax=max(mergevax(:,2:end)-mergevax(:,1:end-1),0);
    
    for i=1:length(listoflocs)
        dailyvax7(i,:)=movmean(dailyvax(i,:),[6 0]);
    end
    
    tlagvax=24;
    interimvaxedS=zeros(num_mp,tlagvax,num_ens);%vaccinated susceptibles per timesetp
    totalvaxperloc=0;%cumulative doses per subpop
    total_vaccine_immune(num_mp,num_ens)=0;
    
    vaxbystate=[];
    %% smooth the data:  7 day moving average
 
    for l=1:num_loc
        obs_case(l,:)=movmean(dailyincidence(l,:),[6 0]);
        obs_death(l,:)=movmean(dailydeaths(l,:),[6 0]);
    end
    
    %%
    incidence=dailyincidence;%to initialize
    startday='02/21/2020'
    %set OEV
    OEV_case=zeros(size(dailyincidence));
    
    for l=1:num_loc
        for t=1:num_times
            obs_ave=mean(dailyincidence(l,max(1,t-6):t));
            OEV_case(l,t)=max(5,obs_ave^2/25);
            
        end
    end
    
    %adjusting inter-county movement
    load MI_inter
    %adjusting mobility starting from March 16, day 25
    MI_inter_relative=MI_inter(:,2:end);
    for t=25:size(MI_inter_relative,2)
        MI_inter_relative(:,t)=MI_inter_relative(:,t)./MI_inter_relative(:,t-1);
        MI_inter_relative(isnan(MI_inter_relative(:,t)),t)=0;
        MI_inter_relative(isinf(MI_inter_relative(:,t)),t)=0;
        MI_inter_relative(:,t)=min(MI_inter_relative(:,t),1);
    end
    MI_inter_relative(:,1:24)=1;
    
    C=C*ones(1,T);
    Cave=Cave*ones(1,T);
    for t=25:T%day 25: March 16
        C(:,t)=C(:,t-1);
        Cave(:,t)=Cave(:,t-1);
        for l=1:num_loc
            for j=part(l)+1:part(l+1)-1
                if t<=size(MI_inter_relative,2)
                    C(part(l),t)=C(part(l),t)+((1-MI_inter_relative(nl(j),t))*C(j,t));
                    Cave(part(l),t)=Cave(part(l),t)+((1-MI_inter_relative(nl(j),t))*Cave(j,t));
                    C(j,t)=(MI_inter_relative(nl(j),t)*C(j,t));
                    Cave(j,t)=(MI_inter_relative(nl(j),t)*Cave(j,t));
                end
            end
        end
    end
    
    
    infection_proj=zeros(num_loc,num_ens,T);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %run model unitl March 10, day 19
    [S,E,Ir,Iu,Seedc]=initialize(nl,part,C(:,1),num_ens,incidence);
    obs_temp=zeros(num_loc,num_ens,T);%records of reported cases
    death_temp=zeros(num_loc,num_ens,T);%records of death
    load('parafit1')
    %initialize parameters
    [para,paramax,paramin,betamap,alphamaps]=initializepara_eakf(dailyincidence,num_ens,parafit);
    
    paramax_ori=paramax;
    paramin_ori=paramin;
    para_ori=para;
    %%%%%%%%%%%%inflate variables and parameters
    lambda=2;
    para=mean(para,2)*ones(1,num_ens)+lambda*(para-mean(para,2)*ones(1,num_ens));
    para=checkbound_para(para,paramax,paramin);
    %     S=mean(S,2)*ones(1,num_ens)+lambda*(S-mean(S,2)*ones(1,num_ens));
    E=mean(E,2)*ones(1,num_ens)+lambda*(E-mean(E,2)*ones(1,num_ens));
    Ir=mean(Ir,2)*ones(1,num_ens)+lambda*(Ir-mean(Ir,2)*ones(1,num_ens));
    Iu=mean(Iu,2)*ones(1,num_ens)+lambda*(Iu-mean(Iu,2)*ones(1,num_ens));
    [S,E,Ir,Iu]=checkbound(S,E,Ir,Iu,C(:,1));
    %fix Z, D and theta
    %Z
    para(1,:)=parafit(3,1:num_ens);
    %D
    para(2,:)=parafit(4,1:num_ens);
    %mu
    para(3,:)=parafit(2,1:num_ens);
    %theta
    para(4,:)=parafit(6,1:num_ens);
    
    parastd=std(para,0,2);
    %auxiliary
    tempcase=zeros(1,1,T);
    tempdeath=zeros(1,1,T);
    if startfromlast==0%start from beginning
        for t=1:22%run unitl March 13
            t
            %%%%%%%%%%%seeding
            if t<=size(Seedc,2)
                [S,E,Ir,Iu]=seeding(S,E,Ir,Iu,nl,part,C(:,t),Seedc,t);
            end
            for k=1:num_ens
                [S(:,k),E(:,k),Ir(:,k),Iu(:,k),dailyIr_temp,dailyIu_temp]=model_eakf(nl,part,C(:,t),Cave(:,t),S(:,k),E(:,k),Ir(:,k),Iu(:,k),para(:,k),betamap,alphamaps);
                for l=1:num_loc
                    infection_proj(l,k,t)=min(sum(dailyIr_temp(part(l):part(l+1)-1))+sum(dailyIu_temp(part(l):part(l+1)-1)),500);
                end
                %reporting delay
                for l=1:num_loc
                    for j=part(l):part(l+1)-1
                        inci=min(dailyIr_temp(j),100);
                        if inci>0
                            rnd=datasample(rnds,inci);
                            for h=1:length(rnd)
                                if (t+rnd(h)<=T)
                                    obs_temp(l,k,t+rnd(h))=obs_temp(l,k,t+rnd(h))+1;
                                end
                            end
                        end
                    end
                end
                
            end
        end
    end
    
    cumu_dailyIr_post=zeros(num_mp,num_ens);
    cumu_dailyIu_post=zeros(num_mp,num_ens);
    
    %EAKF
    lambda=1.2;%inflation
    num_para=size(para,1);

    
    if startfromlast==0
        Tlast=22;
    else
        if   startfromlastvax==0
            load (['../../',lastsavefolder,'/output/temp/lastinference.mat']);
            para_post(:,:,1:size(para_post_last,3))=para_post_last;
        else
            load (['../../',lastsavefolder,'/origoutput/temp/lastinferenceVax.mat']);
        end
        S=S_last; E=E_last; Ir=Ir_last; Iu=Iu_last;
        obs_temp(:,:,1:size(obs_temp_last,3))=obs_temp_last;
        death_temp(:,:,1:size(death_temp_last,3))=death_temp_last;
        infection_proj(:,:,1:size(infection_proj_last,3))=infection_proj_last;
        para=para_last;
        if savevaxstart==1
            num_times=299;
        end
        
        clear *_last
    end
    
    for t=Tlast+1:num_times%start from March 14, and stop on the latest date
        
        tic
        paramin(4+(1:num_loc))=min(0.05+0.2*(t-22)/223,0.25);
        
        %fix Z, D and theta
        %Z
        para(1,:)=parafit(3,1:num_ens);
        %D
        para(2,:)=parafit(4,1:num_ens);
        %theta
        para(4,:)=parafit(6,1:num_ens);
        %seeding
        if t<=size(Seedc,2)
            [S,E,Ir,Iu]=seeding(S,E,Ir,Iu,nl,part,C(:,t),Seedc,t);
        end
        %re-initialize beta
        for l=1:num_loc
            if Seedc(l,t)>0
                para(betamap(l),:)=para_ori(betamap(l),:);
            end
        end
        %%%%%%%%%%%%%%%%%%
        S_temp=S; E_temp=E; Ir_temp=Ir; Iu_temp=Iu;
        %integrate forward one step
        dailyIr_prior=zeros(num_mp,num_ens);
        dailyIu_prior=zeros(num_mp,num_ens);
        for k=1:num_ens%run for each ensemble member
            [S(:,k),E(:,k),Ir(:,k),Iu(:,k)]=adjustmobility(S(:,k),E(:,k),Ir(:,k),Iu(:,k),nl,part,MI_inter_relative,t);
            [S(:,k),E(:,k),Ir(:,k),Iu(:,k),dailyIr_temp,dailyIu_temp]=model_eakf(nl,part,C(:,t),Cave(:,t),S(:,k),E(:,k),Ir(:,k),Iu(:,k),para(:,k),betamap,alphamaps);
            dailyIr_prior(:,k)=dailyIr_temp;
            dailyIu_prior(:,k)=dailyIu_temp;
        end
        
        %% vaccinate
        % vaccine starts Dec 14, or t=298
        
        
        dv=datenum(startday)+t -1;
        dayssincevax=t-297;
        txeff=0.90; %transmisison blocking efficacy
        
        if dayssincevax > 0
            
            for i=1:51
                ind=find(subpopstatemap==statefips.fips(i));
                if dayssincevax<=size(dailyvax7,2)
                    vaxperloc(ind,1)=dailyvax7(i,dayssincevax)*C(ind,t)/sum(C(ind,t));
                else
                    vaxperloc(ind,1)=dailyvax7(i,end)*C(ind,t)/sum(C(ind,t));
                end
            end
            vaxbystate(:,dayssincevax)=accumarray(subpopstatemap,vaxperloc);
            totalvaxperloc=totalvaxperloc+vaxperloc; %track cumulative doses for subpopulation
            unvaccinated=max(round((C(:,t)-totalvaxperloc)),0);
            
            
            for k=1:num_ens%run for each ensemble member
                
                % proportion of vaccines going to susecptible=
                % S/(N - already vaccinated)
                
                proportion_Seligible=(S(:,k)- squeeze(sum(interimvaxedS(:,:,k),2)))./unvaccinated;
                proportion_Seligible(unvaccinated==0)=0; % avoid divide by zero
                proportion_Seligible=max(proportion_Seligible,0); % avoid negtive
                proportion_Seligible=min(proportion_Seligible,1);
                
                interimvaxedS(:,tlagvax+1,k)=proportion_Seligible.*vaxperloc;
                
                % 3.5 weeks after vaccine, move successfully vaccinated
                % from S to R
                if dayssincevax>tlagvax
                    S(:,k)=max(S(:,k)-interimvaxedS(:,1,k)*txeff,0);
                    % store total effectively vaccinated to date
                    total_vaccine_immune(:,k)=total_vaccine_immune(:,k)+interimvaxedS(:,1,k)*txeff;
                end
                
            end
            interimvaxedS=interimvaxedS(:,2:end,:);% remove the first date, since those people have already mvoed out of S
            for goat=1:num_ens
                tempvs(:,goat)=accumarray(subpopstatemap,interimvaxedS(:,end,goat));
            end
            newlyvaxedSbystate(:,dayssincevax)=median(tempvs,2);
            trackSbystate(:,dayssincevax)= accumarray(subpopstatemap,median(S,2));
            
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%
        %integrate forward for six days, prepare for observation
        Tproj=6;
        obs_temp1=obs_temp;
        death_temp1=death_temp;
        for t1=t:t+Tproj-1
            for k=1:num_ens%run for each ensemble member
                [S_temp(:,k),E_temp(:,k),Ir_temp(:,k),Iu_temp(:,k)]=adjustmobility(S_temp(:,k),E_temp(:,k),Ir_temp(:,k),Iu_temp(:,k),nl,part,MI_inter_relative,t1);
                [S_temp(:,k),E_temp(:,k),Ir_temp(:,k),Iu_temp(:,k),dailyIr_temp,dailyIu_temp]=model_eakf(nl,part,C(:,t1),Cave(:,t1),S_temp(:,k),E_temp(:,k),Ir_temp(:,k),Iu_temp(:,k),para(:,k),betamap,alphamaps);
                %reporting delay
                for l=1:num_loc
                    for j=part(l):part(l+1)-1
                        inci=dailyIr_temp(j);
                        if inci>0
                            rnd=datasample(rnds,inci);
                            for h=1:length(rnd)
                                if (t1+rnd(h)<=T)
                                    obs_temp1(l,k,t1+rnd(h))=obs_temp1(l,k,t1+rnd(h))+1;
                                end
                            end
                        end
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for t1=min(t+Tproj,num_times)
            obs_ens=obs_temp1(:,:,t1);%observation at t1, prior
            
            %loop through local observations
            for l=1:num_loc
                
                %% %%%%%%%%%%%%%%%%%case
                
                if (sum(obs_case(l,1:t1))>0)
                    %Get the variance of the ensemble
                    obs_var = OEV_case(l,t1);
                    prior_var = var(obs_ens(l,:));
                    post_var = prior_var*obs_var/(prior_var+obs_var);
                    if prior_var==0%if degenerate
                        post_var=1e-3;
                        prior_var=1e-3;
                    end
                    prior_mean = mean(obs_ens(l,:));
                    post_mean = post_var*(prior_mean/prior_var + obs_case(l,t1)/obs_var);
                    %%%% Compute alpha and adjust distribution to conform to posterior moments
                    alpha = (obs_var/(obs_var+prior_var)).^0.5;
                    dy = post_mean + alpha*(obs_ens(l,:)-prior_mean)-obs_ens(l,:);
                    %Loop over each state variable (connected to location l)
                    %adjust related metapopulation
                    neighbors=part(l):part(l+1)-1;%metapopulation live in l
                    
                    for h=1:length(neighbors)
                        j=neighbors(h);
                        %S
                        %                         temp=S(j,:);
                        %                         A=cov(temp,obs_ens(l,:));
                        %                         rr=A(2,1)/prior_var;
                        %                         dx=rr*dy;
                        %                         S(j,:)=S(j,:)+dx;
                        %E
                        temp=E(j,:);
                        A=cov(temp,obs_ens(l,:));
                        rr=A(2,1)/prior_var;
                        dx=rr*dy;
                        E(j,:)=E(j,:)+dx;
                        %Ir
                        temp=Ir(j,:);
                        A=cov(temp,obs_ens(l,:));
                        rr=A(2,1)/prior_var;
                        dx=rr*dy;
                        Ir(j,:)=Ir(j,:)+dx;
                        %Iu
                        temp=Iu(j,:);
                        A=cov(temp,obs_ens(l,:));
                        rr=A(2,1)/prior_var;
                        dx=rr*dy;
                        Iu(j,:)=Iu(j,:)+dx;
                        %dailyIr
                        temp=dailyIr_prior(j,:);
                        A=cov(temp,obs_ens(l,:));
                        rr=A(2,1)/prior_var;
                        dx=rr*dy;
                        dailyIr_prior(j,:)=round(max(dailyIr_prior(j,:)+dx,0));
                        %dailyIu
                        temp=dailyIu_prior(j,:);
                        A=cov(temp,obs_ens(l,:));
                        rr=A(2,1)/prior_var;
                        dx=rr*dy;
                        dailyIu_prior(j,:)=round(max(dailyIu_prior(j,:)+dx,0));
                    end
                    
                    temp=para(alphamaps(l),:);
                    A=cov(temp,obs_ens(l,:));
                    rr=A(2,1)/prior_var;
                    dx=rr*dy;
                    para(alphamaps(l),:)=para(alphamaps(l),:)+dx;
                    %inflation
                    if std(para(alphamaps(l),:))<parastd(alphamaps(l))
                        para(alphamaps(l),:)=mean(para(alphamaps(l),:),2)*ones(1,num_ens)+lambda*(para(alphamaps(l),:)-mean(para(alphamaps(l),:),2)*ones(1,num_ens));
                    end
                    %adjust beta
                    temp=para(betamap(l),:);
                    A=cov(temp,obs_ens(l,:));
                    rr=A(2,1)/prior_var;
                    dx=rr*dy;
                    para(betamap(l),:)=para(betamap(l),:)+dx;
                    %inflation
                    if std(para(betamap(l),:))<parastd(betamap(l))
                        para(betamap(l),:)=mean(para(betamap(l),:),2)*ones(1,num_ens)+lambda*(para(betamap(l),:)-mean(para(betamap(l),:),2)*ones(1,num_ens));
                    end
                else
                    %adjust related metapopulation
                    neighbors=part(l):part(l+1)-1;%metapopulation live in l
                    for h=1:length(neighbors)
                        j=neighbors(h);
                        %S
                        S(j,:)=C(j,t);
                        %E
                        E(j,:)=0;
                        %Ir
                        Ir(j,:)=0;
                        %Iu
                        Iu(j,:)=0;
                        %dailyIr
                        dailyIr_prior(j,:)=0;
                        %dailyIu
                        dailyIu_prior(j,:)=0;
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            para=checkbound_para(para,paramax,paramin);
        end
        %% update posterior Ir and Iu
        dailyIr_post=dailyIr_prior;
        dailyIu_post=dailyIu_prior;
        
        cumu_dailyIr_post=cumu_dailyIr_post+dailyIr_post;
        cumu_dailyIu_post=cumu_dailyIu_post+dailyIu_post;
        
        if dayssincevax > 0
            for goat=1:num_ens
                tempnewIr(:,goat)=accumarray(subpopstatemap,dailyIr_post(:,goat));
            end
            newIrbystate(:,dayssincevax,1:3)=prctile(tempnewIr,[2.5 50 97.5],2);
            trackIrbystate(:,dayssincevax)= accumarray(subpopstatemap,median(dailyIr_post,2));
            trackIubystate(:,dayssincevax)= accumarray(subpopstatemap,median(dailyIu_post,2));
        end
        
        %% update obs_temp
        for k=1:num_ens
            for l=1:num_loc
                for j=part(l):part(l+1)-1
                    inci=dailyIr_post(j,k);
                    if inci>0
                        rnd=datasample(rnds,inci);
                        for h=1:length(rnd)
                            if (t+rnd(h)<=T)
                                obs_temp(l,k,t+rnd(h))=obs_temp(l,k,t+rnd(h))+1;
                            end
                        end
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%
        S=C(:,t)*ones(1,num_ens)-E-cumu_dailyIr_post-cumu_dailyIu_post-total_vaccine_immune;
        
        [S,E,Ir,Iu]=checkbound(S,E,Ir,Iu,C(:,t));
        if dayssincevax>0
            para_post(:,:,dayssincevax)=para;
        end
        %%%%prepare projection
        for k=1:num_ens
            %total infection
            for l=1:num_loc
                infection_proj(l,k,t)=sum(dailyIr_post(part(l):part(l+1)-1,k))+sum(dailyIu_post(part(l):part(l+1)-1,k));
            end
        end
        
        %%
        %if
        if dayssincevax==1 % save a last inference here
            
            Tlast=t;
            S_last=S; E_last=E; Ir_last=Ir; Iu_last=Iu;
            obs_temp_last=obs_temp;%records of reported cases
            death_temp_last=death_temp;%records of reported cases
            infection_proj_last=infection_proj;
            
            para_last=para;
            save('../output/temp/lastinferenceVax.mat','S_last','E_last','Ir_last','Iu_last',...
                'obs_temp_last','death_temp_last','infection_proj_last',...
                'para_last','Tlast','interimvaxedS','stateR0','vaxbystate',...
                'cumu_dailyIr_post','cumu_dailyIu_post',...
                'totalvaxperloc',   'total_vaccine_immune','newlyvaxedSbystate','newIrbystate');% add back in 'dec14stateS'
            clear *_last
        elseif  t== (num_times-Tsave)
            % save a last inference here
            
            Tlast=t;
            S_last=S; E_last=E; Ir_last=Ir; Iu_last=Iu;
            obs_temp_last=obs_temp;%records of reported cases
            death_temp_last=death_temp;%records of reported cases
            infection_proj_last=infection_proj;
            para_post_last=para_post;
            para_last=para;
            save('../output/temp/lastinference.mat','S_last','E_last','Ir_last','Iu_last',...
                'obs_temp_last','death_temp_last','infection_proj_last',...
                'para_post_last','para_last','Tlast','interimvaxedS','stateR0','vaxbystate',...
                'cumu_dailyIr_post','cumu_dailyIu_post',...
                'totalvaxperloc',   'total_vaccine_immune','newlyvaxedSbystate','newIrbystate');% add back in 'dec14stateS'
            clear *_last
        end
        %% save S on the first day of vaccine
        if dayssincevax==0
            % first day of vaccine save S
            for ee=1:100
                dec14stateS(:,ee)=accumarray(subpopcountymap,S(:,ee));
            end
            countyfips.dec14stateS=median(dec14stateS,2);
            
            writetable(countyfips,'../output/marta/dec14S.csv');
        end
        %% some state level outputs
        for j=1:100
   
            mbeta=para(betamap,:);
            malpha=para(alphamaps,:);
            mD=para(2,:);
            mmu=para(3,:);
            R0(:,j)=median(mbeta(:,j).*mD(j).*(malpha(:,j)+(1-malpha(:,j)).*mmu(j)),2);
        end
        medR0=median(R0,2);
        
        % get state level R0 by taking popuulatoin weighted avg
        for s=unique(countystatemap)';
            stateR0(s,t)=sum(medR0(countystatemap==s).*...
                population(countystatemap==s)/sum(population(countystatemap==s)));
            
        end
        %%
        if t==num_times
            % last day for JHU scenario fitting
            mkdir('../output/marta');
            csvwrite('../output/marta/stater0daily.csv',...
                [statefips.fips stateR0(statefips.fips,:)])

            for j=1:num_ens
                Sstate(:,j)=accumarray(subpopstatemap,S(:,j));
                Estate(:,j)=accumarray(subpopstatemap,E(:,j));
                Irstate(:,j)=accumarray(subpopstatemap,Ir(:,j));
                Iustate(:,j)=accumarray(subpopstatemap,Iu(:,j));
                total_vaccine_immune_state(:,j)=accumarray(subpopstatemap,total_vaccine_immune(:,j));
                %  statedec14s(:,j)=accumarray(countystatemap,dec14stateS(:,j));
            end
            
            statetable=[array2table([1:56]','VariableNames',{'fips'}),array2table(round(newlyvaxedSbystate))];
            statetable2=join(statefips,statetable,'Keys','fips','LeftVariables',{'state_name','fips'});
            writetable(statetable2,'../output/marta/newlyvaccinatedS.csv');
            
            newIrbystate=newIrbystate(statefips.fips,:,:);
            save('../output/marta/newIrbystate.mat','newIrbystate');
            
            statetable=[array2table([1:56]','VariableNames',{'fips'})...
                array2table(round(vaxbystate))];
            statetable3=join(statefips,statetable,'Keys','fips','LeftVariables',{'state_name','fips'});
            writetable(statetable3,'../output/marta/totalvaxbystate.csv');
            
            statetable=[array2table([1:56]','VariableNames',{'fips'})...
                array2table(round(median(total_vaccine_immune_state,2)))];
            statetable3=join(statefips,statetable,'Keys','fips','LeftVariables',{'state_name','fips'});
            writetable(statetable3,'../output/marta/total_vaccine_immunebystate.csv');
            
            writestate(statefips,Sstate,Estate,Irstate,Iustate);
            
        end
        
        toc
    end
    
    S_start=S; E_start=E; Ir_start=Ir; Iu_start=Iu;
    obs_temp_start=obs_temp;%records of reported cases
    death_temp_start=death_temp;%records of reported cases
    infection_proj_start=infection_proj;
    
    %clear interimvaxedS
    clear *_last
    clear *_ori
    clear *_prior
    clear *_temp
    clear *_temp1
    clear S E Ir Iu
    if savevaxstart==1
        save('../output/temp/saveeverythingVax.mat');
    else
        save('../output/baselinetemp/saveeverythingforpaper.mat');
    end     
 
return







function para = checkbound_paraini(para,paramax,paramin)
for i=1:size(para,1)
    temp=para(i,:);
    index=(temp<paramin(i))|(temp>paramax(i));
    index_out=find(index>0);
    index_in=find(index==0);
    %redistribute out bound ensemble members
    para(i,index_out)=datasample(para(i,index_in),length(index_out));
end


function para = checkbound_para(para,paramax,paramin)
for i=1:size(para,1)
    para(i,para(i,:)<paramin(i))=paramin(i)*(1+0.2*rand(sum(para(i,:)<paramin(i)),1));
    para(i,para(i,:)>paramax(i))=paramax(i)*(1-0.2*rand(sum(para(i,:)>paramax(i)),1));
end

function [S,E,Ir,Iu]=checkbound(S,E,Ir,Iu,C)
for k=1:size(S,2)
    S(S(:,k)<0,k)=0; E(E(:,k)<0,k)=0; Ir(Ir(:,k)<0,k)=0; Iu(Iu(:,k)<0,k)=0;
end


function [S,E,Ir,Iu]=adjustmobility(S,E,Ir,Iu,nl,part,MI_inter_relative,t)
num_loc=size(MI_inter_relative,1);
for l=1:num_loc
    for j=part(l)+1:part(l+1)-1
        if t<=size(MI_inter_relative,2)
            S(part(l))=S(part(l))+((1-MI_inter_relative(nl(j),t))*S(j));
            S(j)=(MI_inter_relative(nl(j),t)*S(j));
            E(part(l))=E(part(l))+((1-MI_inter_relative(nl(j),t))*E(j));
            E(j)=(MI_inter_relative(nl(j),t)*E(j));
            Ir(part(l))=Ir(part(l))+((1-MI_inter_relative(nl(j),t))*Ir(j));
            Ir(j)=(MI_inter_relative(nl(j),t)*Ir(j));
            Iu(part(l))=Iu(part(l))+((1-MI_inter_relative(nl(j),t))*Iu(j));
            Iu(j)=(MI_inter_relative(nl(j),t)*Iu(j));
        end
    end
end

function [obsdataout ]=postprocess_proj(obsdata,...
    num_loc,num_ens,obs_case, num_times )

temp=obsdata;
for l=1:num_loc
    total=squeeze(sum(temp(l,:,:),3))';
    total(:,2)=(1:num_ens)';
    total=sortrows(total,1);%from low to high
    obsdataout(l,:,:)=temp(l,total(:,2),:);
    
end

 




function [ obsdata]= delayedobservations( obsdata, num_loc,part,...
    dailyIr_temp,dailyIu_temp ,k,t,rnds,rnds_d,T)
 

%reporting delay
for l=1:num_loc
    for j=part(l):part(l+1)-1
        inci=dailyIr_temp(j);
        if inci>0
            rnd=datasample(rnds,inci);
            for h=1:length(rnd)
                if (t+rnd(h)<=T)
                    obsdata(l,k,t+rnd(h))=obsdata(l,k,t+rnd(h))+1;
                end
            end
        end
    end
end

 
obsdata=round(obsdata);
 
 