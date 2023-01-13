%% this is the code that runs the counterfactual projections. it uses as inputs files that are created in infer.m


clear;close all;
load('../output/temp/saveeverythingVax.mat');
futurepara=load('../output/baselinetemp/saveeverythingforpaper.mat','para_post');
startday='02/21/2020';
%% intervention scenarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% intervention

avertstarttime=t; 
T=size(dailyincidence,2); % set to num_times if comparing to observed. extend to size(obs_temp_start,3) to do projection as well.
load('countyfips.mat');
 dayssincevax=1;


%% Counterfactual: Scenario 1

increaserate=0;
reduce=0;
S=S_start; E=E_start; Iu=Iu_start; Ir=Ir_start;
obs_proj=obs_temp_start;

infection_proj=infection_proj_start;

wk=0;
for t=avertstarttime+1:T

    dayssincevax=dayssincevax+1;
    parapost=futurepara.para_post(:,:,dayssincevax);
    paratemp=parapost;
    beta=paratemp(betamap,:);
    alpha=paratemp(alphamaps,:);
    D=paratemp(2,:);
    tempmu=paratemp(3,:);


    tdiff=t-(num_times+1);
    if mod(tdiff,7)==0
        % save Reff
        wk=wk+1;
        for i=1:3142
            for j=1:num_ens
                Reff(i,j)=beta(i,j).*D(j).*(alpha(i,j)+(1-alpha(i,j)).*tempmu(j)).*sum(S(part(i):part(i+1)-1,j))...
                    ./population(i);
            end
        end
        Reff_nochange(:,:,wk)=prctile(Reff,[2.5 25 50 75 97.5],2);
    end



    %%
    for k=1:num_ens%run for each ensemble member

        %% continue with SEIR model
        [S(:,k),E(:,k),Ir(:,k),Iu(:,k)]=adjustmobility(S(:,k),E(:,k),Ir(:,k),Iu(:,k),nl,part,MI_inter_relative,t);
        [S(:,k),E(:,k),Ir(:,k),Iu(:,k),dailyIr_temp,dailyIu_temp]=model_eakf(nl,part,C(:,t),Cave(:,t),...
            S(:,k),E(:,k),Ir(:,k),Iu(:,k),paratemp(:,k),betamap,alphamaps);
        for l=1:num_loc
            infection_proj(l,k,t)=sum(dailyIr_temp(part(l):part(l+1)-1))+sum(dailyIu_temp(part(l):part(l+1)-1));
        end
        [obs_proj ]= ...
            delayedobservations( obs_proj,  num_loc,part,...
            dailyIr_temp,dailyIu_temp,k,t,rnds,rnds_d,T);

    end

end

%%

[obs_projnochange ]=postprocess_proj(obs_proj,...
    num_loc,num_ens,obs_case,num_times );

infection_projnochange=infection_proj;

save(['../output/temp/counterfact'],'obs_case',...
    'obs_projnochange','infection_projnochange' );


clear *projnochange 

%% Counterfactual: Scenario 3
dominus10=1
if dominus10==1

    S=S_start; E=E_start; Iu=Iu_start; Ir=Ir_start;
    obs_proj=obs_temp_start;
    
    infection_proj=infection_proj_start;
    
    dayssincevax=1;
    wk=0;
    for t=avertstarttime+1:T
        
        dayssincevax=dayssincevax+1;
        parapost=futurepara.para_post(:,:,dayssincevax);
        paratemp=parapost;
        paratemp(betamap,:)=paratemp(betamap,:)*0.9; % minus 10%;
        beta=paratemp(betamap,:);
        alpha=paratemp(alphamaps,:);
        D=paratemp(2,:);
        tempmu=paratemp(3,:);
        
        
        tdiff=t-(num_times+1);
        if mod(tdiff,7)==0
            % save Reff
            wk=wk+1;
            for i=1:3142
                for j=1:num_ens
                    Reff(i,j)=beta(i,j).*D(j).*(alpha(i,j)+(1-alpha(i,j)).*tempmu(j)).*sum(S(part(i):part(i+1)-1,j))...
                        ./population(i);
                end
            end
            Reff_minus10(:,:,wk)=prctile(Reff,[2.5 25 50 75 97.5],2);
        end
        
        
        
        %%
        for k=1:num_ens%run for each ensemble member
            
            %% continue with SEIR model
            [S(:,k),E(:,k),Ir(:,k),Iu(:,k)]=adjustmobility(S(:,k),E(:,k),Ir(:,k),Iu(:,k),nl,part,MI_inter_relative,t);
            [S(:,k),E(:,k),Ir(:,k),Iu(:,k),dailyIr_temp,dailyIu_temp]=model_eakf(nl,part,C(:,t),Cave(:,t),...
                S(:,k),E(:,k),Ir(:,k),Iu(:,k),paratemp(:,k),betamap,alphamaps);
            for l=1:num_loc
                infection_proj(l,k,t)=sum(dailyIr_temp(part(l):part(l+1)-1))+sum(dailyIu_temp(part(l):part(l+1)-1));
            end
            [obs_proj ]= ...
                delayedobservations( obs_proj,  num_loc,part,...
                dailyIr_temp,dailyIu_temp,k,t,rnds,rnds_d,T);
            
        end
    end
    
    %%
    
    [obs_projminus10 ]=postprocess_proj(obs_proj,...
        num_loc,num_ens,obs_case,num_times );
    
    infection_projminus10=infection_proj;
    

    save(['../output/temp/counterfact_minus10'],'obs_case' ,...
        'obs_projminus10','infection_projminus10' );
end
%% Counterfactual: Scenario 2
doplus10=1
if doplus10==1
    S=S_start; E=E_start; Iu=Iu_start; Ir=Ir_start;
    obs_proj=obs_temp_start;
    
    infection_proj=infection_proj_start;
    dayssincevax=1;
    wk=0;
    for t=avertstarttime+1:T
        
        dayssincevax=dayssincevax+1;
        parapost=futurepara.para_post(:,:,dayssincevax);
        paratemp=parapost;
        paratemp(betamap,:)=paratemp(betamap,:)*1.1; % plus 10%;
        beta=paratemp(betamap,:) ;
        alpha=paratemp(alphamaps,:);
        D=paratemp(2,:);
        tempmu=paratemp(3,:);
                
        tdiff=t-(num_times+1);
        if mod(tdiff,7)==0
            % save Reff
            wk=wk+1;
            for i=1:3142
                for j=1:num_ens
                    Reff(i,j)= beta(i,j).*D(j).*(alpha(i,j)+(1-alpha(i,j)).*tempmu(j)).*sum(S(part(i):part(i+1)-1,j))...
                        ./population(i);
                end
            end
            Reff_plus10(:,:,wk)=prctile(Reff,[2.5 25 50 75 97.5],2);
        end
                        
        %%
        for k=1:num_ens%run for each ensemble member            
            % continue with SEIR model
            [S(:,k),E(:,k),Ir(:,k),Iu(:,k)]=adjustmobility(S(:,k),E(:,k),Ir(:,k),Iu(:,k),nl,part,MI_inter_relative,t);
            [S(:,k),E(:,k),Ir(:,k),Iu(:,k),dailyIr_temp,dailyIu_temp]=model_eakf(nl,part,C(:,t),Cave(:,t),...
                S(:,k),E(:,k),Ir(:,k),Iu(:,k),paratemp(:,k),betamap,alphamaps);
            for l=1:num_loc
                infection_proj(l,k,t)=sum(dailyIr_temp(part(l):part(l+1)-1))+sum(dailyIu_temp(part(l):part(l+1)-1));
            end
            [obs_proj ]= ...
                delayedobservations( obs_proj,  num_loc,part,...
                dailyIr_temp,dailyIu_temp,k,t,rnds,rnds_d,T);            
        end        
    end
    

    [obs_projplus10 ]=postprocess_proj(obs_proj,...
        num_loc,num_ens,obs_case,num_times );    
    infection_projplus10=infection_proj;

    save(['../output/temp/counterfact_plus10'],'obs_case',...
        'obs_projplus10','infection_projplus10' );
end
 


%% functions
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

end
 