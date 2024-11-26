%SEIR-EKAF MODEL RUN WITH REAL DATA
function MODEL_RUN(nickname)

%Inter-county commuting is stored using neighbor list: nl, part and C.
%nl (neighbor list) and part (partition) both describe the network structure.
%For instance, the neighbors of location i are nl(part(i):part(i+1)-1).
%This neighbor set include i itself.
%C has the same size with nl.
%The commuters from location i to location j=nl(part(i)+x) is C(part(i)+x).

num_times=437; %from January 20, 2020 to  March 31, 2021

%%%%%%%%%%%%%%%%%loading data

load dailyincidence_real.csv
dailyincidence=dailyincidence_real(:,1:num_times);
clear dailyincidence_real
load commutedata.mat %load inter-county commuting data
load statecodes.mat
load population.mat
load parafit_rev.mat
load google_percw.csv

loc_of_interest = {'Baja California','California','Texas','Florida','New York','Pennsylvania','Estado de Mexico', ...
    'Distrito Federal','Jalisco','Veracruz de Ignacio de la Llave','Puebla', ...
    'Ontario','Quebec','British Columbia','Alberta','Manitoba'};


%reproductive number R
%Rt ≡ βD[α + (1 − α)μ]


%%%%%%%%%%%%%%%%%%%%%%%%%
num_ens=300; %number of ensemble members

%%%%%%%%%%%%%%%%%%%%% select today for figures names
mmdd=datestr(datetime('now'), 'mmdd');
all_file_name=strjoin([ "Model_Runs/" mmdd "_" nickname ".mat"],'');
lambda=1.01;%inflation in EAKF
lambda_beta=lambda;

%inflaction factor for the initial parameter values
lambda_initial_para=0.5;

%inflation factor for the values of out-of-bound para in checkbound_para.m
flact_checkpara=0;

%inflate observed variable, yes/no
inflate_obs = "yes";
%inflate all the other state variables, yes/no
inflate_sv = "no";
%%%%%%%%%%%%%%%%%%%%% reprobe parameters
doreprobe= "yes";   % yes/no if you're reprobinge -- I reinitialize alpha and beta, not reprobe
reprobe_percent=2;  % what percent of the ensemble to reprobe
reprobe_t=7;        % how often to reprobe (days)
reprobeS= "no";     % yes/no if you're reprobing State Variables


%reintialize a fraction of the alphas to the prior values +/- 10%
a_reinit_time=30000000; %no reeinitializiaton
a_reinit_perc=0; %percentage of ensamble to reinitialize

%increasing alphamin gradually
amin_increase=0.5; %percentence of daily increase for alpha min in checbound_para exponentioal
amin_increase_mex=0; %percentence of daily increase for alpha min in checbound_para linear for mx (no increase)
l_can = table2array(statecodes(strcmp(statecodes.Var3, 'CA'), 1));
l_usa = table2array(statecodes(strcmp(statecodes.Var3, 'US'), 1));
l_mex = table2array(statecodes(strcmp(statecodes.Var3, 'MX'), 1));


%%%%% OEV settings, OEV_case(l,t)=max(OEV_base,obs_ave^OEV_exp);
OEV_denom=100;
OEV_base=5;

num_loc=size(part,1)-1;
num_mp=size(nl,1);
num_times=size(dailyincidence,2);


%%%%%% GOOGLE ADJ %%%%%%

%the real_dailyincidence start from January 22nd while
%google adj starts form February 15th
%-there are 24 days between January 22nd and February 15th.-
google_percw=google_percw/100; %from 0 to 1
google_percw=[zeros(num_loc, 24), google_percw]; % adding 24 leading 0 to the dataset for each location
google_percw=google_percw+1; % easier moltiplation and avoids negatives
goog_lenght=size(google_percw,2);

%smooth the data: 7 day moving average
for l=1:num_loc
    for t=1:goog_lenght
        if (t+6)<=goog_lenght
            google_percw(l,t)=(mean(google_percw(l,max(1,t-6):min(t+6,goog_lenght))));
        else
            google_percw(l,t)=(mean(google_percw(l,max(1,goog_lenght-6):goog_lenght)));
        end
    end
end


C=C*ones(1,num_times);%Daily subpopulation size
Cave=Cave*ones(1,num_times);
for t=2:num_times
    for l=1:num_loc %state
        for j=part(l)+1:part(l+1)-1 % each metapop, state by state, exluding the first
            if t<=size(google_percw,2)
                C(part(l),t)=C(part(l),t)+(C(j,t)-ceil(google_percw(nl(j),t)*C(j,t))); %big metapop - jj
                Cave(part(l),t)=Cave(part(l),t)+((Cave(j,t)-ceil(google_percw(nl(j),t)*Cave(j,t))));
                C(j,t)=ceil(google_percw(nl(j),t)*C(j,t)); %small metapop - ij
                Cave(j,t)=ceil(google_percw(nl(j),t)*Cave(j,t));
            end
        end
    end
end



%%%%%% GOOGLE ADJ end %%%%%%

obs_case=zeros(size(dailyincidence));
%smooth the data: 7 day moving average
for l=1:num_loc
    for t=1:num_times
        if (t+3)<=num_times
            obs_case(l,t)=(mean(dailyincidence(l,max(1,t-3):min(t+3,num_times))));
        else
            obs_case(l,t)=(mean(dailyincidence(l,max(1,num_times-6):num_times)));
        end
    end
end
T=num_times;%number of days to keep track of reported cases
startday='01/20/2020';
Tstart=datetime(startday);
%set OEV, observation error variance
OEV_case=zeros(size(dailyincidence));
for l=1:num_loc
    for t=1:num_times
        obs_ave=mean(dailyincidence(l,max(1,t-6):t));
        %OEV_case(l,t)=max(5,obs_ave^2/100);
        OEV_case(l,t)=max(OEV_base,(obs_ave^2)/OEV_denom);
    end
end
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize model state variable

S=C(:,1)*ones(1,num_ens);

% seed a few cases in each location
rnd_matrix = randi([1, 9], size(S));
E=rnd_matrix;
rnd_matrix = randi([1, 9], size(S));
Ir=rnd_matrix;
rnd_matrix = randi([1, 9], size(S));
Iu=rnd_matrix;


%initialize parameter
[para,betamap,alphamaps]=initialize_para(num_loc,num_ens,parafit);


%initialize alpha
% alpha = weighted alpha +- random flactuation.
CAid = statecodes.Var1(strcmp(statecodes.Var3, 'CA'));
USid = statecodes.Var1(strcmp(statecodes.Var3, 'US'));
MXid = statecodes.Var1(strcmp(statecodes.Var3, 'MX'));

%US got the  aplha valuee from seroprevalnce studies
% Get the current variable
ids = USid;
fluctuation = lambda_initial_para * US_weighted_alpha;
randomFluctuations = -fluctuation + 2 * fluctuation * rand(size(para(alphamaps(ids), :)));
para(alphamaps(ids), :) = repmat(US_weighted_alpha, 1, num_ens) + randomFluctuations;
%out-of-bound alpha gets US_weighted_alpha
para(alphamaps(ids),para(alphamaps(ids),:)<alphalow)=US_weighted_alpha;
para(alphamaps(ids),para(alphamaps(ids),:)>alphaup)=US_weighted_alpha;

CTRYids = {CAid, MXid};
for i = 1:numel(CTRYids)    
    ids = CTRYids{i};
    mu_exp = 0.022; 
    randomFluctuations = exprnd(mu_exp, size(para(alphamaps(ids), :)));
    para(alphamaps(ids), :) = alphalow + randomFluctuations; % Shifted exponential distribution (min is alpha lower bound)
    para(alphamaps(ids),para(alphamaps(ids),:)>alphaup)=alphalow+0.0001;
end





for i = 1:num_loc
    %out-of-bound beta are resampled across the entire range
    para(betamap(i),para(betamap(i),:)<betalow)=random('uniform', betalow, betaup, size(para(betamap(i),para(betamap(i),:)<betalow)));
    para(betamap(i),para(betamap(i),:)>betaup)=random('uniform', betalow, betaup, size(para(betamap(i),para(betamap(i),:)>betaup)));
end
%out-of-bound Z,D,mu and theta are resampled across their ranges
for i=1:size(para,1)-length(alphamaps)-length(betamap)
    para(i,para(i,:)<paramin(i))=random('uniform', paramin(i), paramax(i), size(para(i,para(i,:)<paramin(i)))); % selecting random values between the bound for each ensamble out
    para(i,para(i,:)>paramax(i))=random('uniform', paramin(i), paramax(i), size(para(i,para(i,:)>paramax(i))));
end


para_ori=para; %used to re-initialize beta
para_amin_increase=para;  %used to re-initialize alpha, the range increase over time
%%%%%%%%%%%%%%%%%%%%% reprobe parameters

reprobe_whichpara=cat(1,alphamaps,betamap);  % reprobe only the pars I am estimating

[S,E,Ir,Iu]=checkbound(S,E,Ir,Iu,C(:,1));

parastd=std(para,0,2);%get ensemble spread of parameters

%create seedc to reinizialize beta before the start of a peak
Seedc=zeros(num_loc,num_times);

reported=sum(dailyincidence,2);%total incidence
reported(:,2)=(1:size(reported,1));
reported=sortrows(reported,-1);

seedlocs=reported(reported(:,1)>0,2);%cases >0

for l=1:size(seedlocs)
    seedloc=seedlocs(l);
    temp=dailyincidence(seedloc,:);
    %find peaks
    below=temp<=1;
    startpoints=[];
    cnt=0;
    last=1;
    for i=1:length(below)
        if below(i)==1
            cnt=cnt+1;
        else
            if cnt>=7
                startpoints=[startpoints,last];
                last=i;
            end
            cnt=0;
        end
    end
    for i=1:length(startpoints)
        index=find(temp>=5);
        index=index(index>=startpoints(i));
        if ~isempty(index)
            T0=index(1);%first reporting with at least 5 cases
            c=sum(temp(T0:min(T0+4,num_times)));
            Seedc(seedloc,max(1,T0-8))=c;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%% inizialize variables


dailyIr_prior_rec=zeros(num_loc,num_ens,T);%prior reported infection
dailyIu_prior_rec=zeros(num_loc,num_ens,T);%prior unreported infection

obs_var_rec=zeros(num_loc,T);
prior_var_rec=zeros(num_loc,T);
prior_mean_rec=zeros(num_loc,T);
post_var_rec=zeros(num_loc,T);
post_mean_rec=zeros(num_loc,T);
alpha_rec=zeros(num_loc,T);

paramin_rec=zeros(size(paramin,1),T);

dy_rec=zeros(num_loc,num_ens,T); %Kalman gain dy
dx_alpha_rec=zeros(num_loc,num_ens,T); %Kalman gain dx alpha
dx_beta_rec=zeros(num_loc,num_ens,T); %Kalman gain dx beta
dx_E_rec=zeros(num_loc,num_ens,T);
dx_Ir_rec=zeros(num_loc,num_ens,T);
dx_Iu_rec=zeros(num_loc,num_ens,T);
dx_dailyIr_prior_rec=zeros(num_loc,num_ens,T);
dx_dailyIu_prior_rec=zeros(num_loc,num_ens,T);




%initialize poseteriors (ends with rec=counts, ends with post=perc)
S_post=zeros(num_loc,num_times,num_ens);
S_rec=zeros(num_loc,num_times,num_ens);
E_post=zeros(num_loc,num_times,num_ens);
E_rec=zeros(num_loc,num_times,num_ens);
Ir_post=zeros(num_loc,num_times,num_ens);
Ir_rec=zeros(num_loc,num_times,num_ens);
Iu_post=zeros(num_loc,num_times,num_ens);
Iu_rec=zeros(num_loc,num_times,num_ens);
dailyIr_post_rec=zeros(num_loc,num_ens,T);
cumu_dailyIr_post_rec=zeros(num_loc,num_times,num_ens);
dailyIu_post_rec=zeros(num_loc,num_ens,T);
cumu_dailyIu_post_rec=zeros(num_loc,num_times,num_ens);


%initialize stavariables for yesterday checkbound, used when t>1
S_yesterday=zeros(num_mp,num_ens);
E_yesterday=zeros(num_mp,num_ens);
Ir_yesterday=zeros(num_mp,num_ens);
Iu_yesterday=zeros(num_mp,num_ens);


%initialize cumulative reported and unreported infections
cumu_dailyIr_post=zeros(num_mp,num_ens);
cumu_dailyIu_post=zeros(num_mp,num_ens);
%%%%%%%%%%%%%%%%%%%%%%%
num_para=size(para,1);%number of parameters
para_post=zeros(num_para,num_ens,num_times);%posterior parameters


para_amin_increase_rec=zeros(size(para_post));
for t=1:num_times
    %increading alphamin over time
    paramin(alphamaps(l_can))=paramin(alphamaps(l_can))*(1+amin_increase/100);
    paramin(alphamaps(l_usa))=paramin(alphamaps(l_usa))*(1+amin_increase/100);
    paramin(alphamaps(l_mex))=paramin(alphamaps(l_mex))+(amin_increase_mex/100);

    paramin_rec(alphamaps,t)=paramin(alphamaps);
    delta_paramin=paramin_rec(alphamaps,t)-paramin_rec(alphamaps,1);
    para_amin_increase(alphamaps,:)=para_ori(alphamaps,:)+delta_paramin; %%used as para_ori to reinitialize alpha in checkbound_para
    para_amin_increase_rec(:,:,t)=para_amin_increase;
    %re-initialize n% of alpha every re_init_time days
    if mod(t,a_reinit_time) == 0
        %select n% of ensamble to reintialize
        n_sel_ens=round(a_reinit_perc/100 * num_ens); %how many ens to select
        sel_ens = randperm(num_ens, n_sel_ens); %which ens are selected
        para(alphamaps, sel_ens) = para_ori(alphamaps, sel_ens);
    end


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
    checkbound_para(para,paramax,paramin,para_amin_increase,alphamaps,betamap,flact_checkpara);



    %integrate forward one step
    dailyIr_prior=zeros(num_mp,num_ens);
    dailyIu_prior=zeros(num_mp,num_ens);

    for k=1:num_ens %run for each ensemble member
        [S(:,k),E(:,k),Ir(:,k),Iu(:,k)]=adjustmobility(S(:,k),E(:,k),Ir(:,k),Iu(:,k),nl,part,google_percw,t);
        [S(:,k),E(:,k),Ir(:,k),Iu(:,k),dailyIr_temp,dailyIu_temp]=integrate_model(nl,part,C(:,t),Cave(:,t),S(:,k),E(:,k),Ir(:,k),Iu(:,k),para(:,k),betamap,alphamaps);
        dailyIr_prior(:,k)=dailyIr_temp;
        dailyIu_prior(:,k)=dailyIu_temp;
    end


    %%%%%% inflate observed
    if inflate_obs == "yes"
        dailyIr_prior=mean(dailyIr_prior,2)*ones(1,num_ens)+lambda*(dailyIr_prior-mean(dailyIr_prior,2)*ones(1,num_ens));
    end
    if inflate_sv == "yes"
        dailyIu_prior=mean(dailyIu_prior,2)*ones(1,num_ens)+lambda*(dailyIu_prior-mean(dailyIu_prior,2)*ones(1,num_ens));
        S=mean(S,2)*ones(1,num_ens)+lambda*(S-mean(S,2)*ones(1,num_ens));
        E=mean(E,2)*ones(1,num_ens)+lambda*(E-mean(E,2)*ones(1,num_ens));
        Ir=mean(Ir,2)*ones(1,num_ens)+lambda*(Ir-mean(Ir,2)*ones(1,num_ens));
        Iu=mean(Iu,2)*ones(1,num_ens)+lambda*(Iu-mean(Iu,2)*ones(1,num_ens));
    end

    for i=1:num_loc
        for j=1:num_ens
            dailyIr_prior_rec(i,j,t)=sum(dailyIr_prior(part(i):part(i+1)-1,j));
            dailyIu_prior_rec(i,j,t)=sum(dailyIu_prior(part(i):part(i+1)-1,j));
        end
    end

    %%%  EAKF  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    obs_ens=dailyIr_prior_rec(:,:,t);
    %loop through local observations
    for l=1:num_loc
        %%%%%%%%%%%%%%%%%%%case
        %Get the variance of the ensemble
        obs_var = OEV_case(l,t);
        obs_var_rec(l,t)=obs_var;

        prior_var = var(obs_ens(l,:));
        post_var = prior_var*obs_var/(prior_var+obs_var);

        if prior_var==0 %if degenerate
            post_var=1e-3;
            prior_var=1e-3;
        end

        prior_mean = mean(obs_ens(l,:));
        post_mean = post_var*(prior_mean/prior_var + obs_case(l,t)/obs_var);


        prior_var_rec(l,t)=prior_var;
        prior_mean_rec(l,t)=prior_mean;
        post_var_rec(l,t)=post_var;
        post_mean_rec(l,t)=post_mean;


        %%%% Compute alpha and adjust distribution to conform to posterior moments
        alpha = (obs_var/(obs_var+prior_var)).^0.5;
        alpha_rec(l,t)=alpha;
        dy = post_mean + alpha*(obs_ens(l,:)-prior_mean)-obs_ens(l,:);
        dy_rec(l,:,t)=dy'; %save Kalman gain for plotting

        %Loop over each state variable (connected to location l)
        %adjust related metapopulation
        neighbors=part(l):part(l+1)-1;%metapopulation live in l
        for h=1:length(neighbors)
            j=neighbors(h);
            %E
            temp=E(j,:);
            A=cov(temp,obs_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            E(j,:)=E(j,:)+dx;
            dx_E_rec(l,:,t)=dx_E_rec(l,:,t)+dx;
            %Ir
            temp=Ir(j,:);
            A=cov(temp,obs_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            Ir(j,:)=Ir(j,:)+dx;
            dx_Ir_rec(l,:,t)=dx_Ir_rec(l,:,t)+dx;
            %Iu
            temp=Iu(j,:);
            A=cov(temp,obs_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            Iu(j,:)=Iu(j,:)+dx;
            dx_Iu_rec(l,:,t)=dx_Iu_rec(l,:,t)+dx;
            %dailyIr
            temp=dailyIr_prior(j,:);
            A=cov(temp,obs_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            dailyIr_prior(j,:)=round(max(dailyIr_prior(j,:)+dx,0));
            dx_dailyIr_prior_rec(l,:,t)=dx_dailyIr_prior_rec(l,:,t)+dx;
            %dailyIu
            temp=dailyIu_prior(j,:);
            A=cov(temp,obs_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            dailyIu_prior(j,:)=round(max(dailyIu_prior(j,:)+dx,0));
            dx_dailyIu_prior_rec(l,:,t)=dx_dailyIu_prior_rec(l,:,t)+dx;
        end
        %adjust alpha
        temp=para(alphamaps(l),:);
        A=cov(temp,obs_ens(l,:));
        rr=A(2,1)/prior_var;
        dx=rr*dy;
        para(alphamaps(l),:)=para(alphamaps(l),:)+dx;

        dx_alpha_rec(l,:,t)=dx'; %save Kalman gain for plotting

        %inflation alpha
        if std(para(alphamaps(l),:))<parastd(alphamaps(l))
            para(alphamaps(l),:)=mean(para(alphamaps(l),:),2)*ones(1,num_ens)+lambda*(para(alphamaps(l),:)-mean(para(alphamaps(l),:),2)*ones(1,num_ens));
        end

        %adjust beta
        temp=para(betamap(l),:);
        A=cov(temp,obs_ens(l,:));
        rr=A(2,1)/prior_var;
        dx=rr*dy;
        para(betamap(l),:)=para(betamap(l),:)+dx;

        dx_beta_rec(l,:,t)=dx'; %save Kalman gain for plotting

        %inflation beta
        if std(para(betamap(l),:))<parastd(betamap(l))
            para(betamap(l),:)=mean(para(betamap(l),:),2)*ones(1,num_ens)+lambda_beta*(para(betamap(l),:)-mean(para(betamap(l),:),2)*ones(1,num_ens));
        end

    end

    para = checkbound_para(para,paramax,paramin,para_amin_increase,alphamaps,betamap,flact_checkpara);

    %%% add reprobing of select parameters

    if doreprobe=="yes"
        if mod(t,reprobe_t)==0 % for every reprobe_t timesteps

            % randomly pick reprobe_percent ensemble members to resample
            num_reprobe=round(num_ens*reprobe_percent/100);
            reprobeind=randi([1 num_ens],num_reprobe,1);

            % for r=reprobe_whichpara % loop through the parameters
            for nn=reprobeind' % loop through the sampled members
                para(reprobe_whichpara,nn)= para_ori(reprobe_whichpara,nn);
                if reprobeS=="yes"
                    S(:,nn)=  ceil(randi([5 95],num_mp,1) /100.*C(:,t));
                    E(:,nn)=  ceil(randi([0 20],num_mp,1) /100.*C(:,t));
                    Ir(:,nn)= ceil(randi([0 20],num_mp,1) /100.*C(:,t));
                    Iu(:,nn)= ceil(randi([0 20],num_mp,1) /100.*C(:,t));

                end
                [S,E,Ir,Iu]=checkbound_yesterday(S,E,Ir,Iu,C(:,t),S_yesterday,E_yesterday,Ir_yesterday,Iu_yesterday,t);
                para = checkbound_para(para,paramax,paramin,para_amin_increase,alphamaps,betamap,flact_checkpara);
            end
        end
    end

    %update posterior Ir and Iu
    dailyIr_post=dailyIr_prior;
    dailyIu_post=dailyIu_prior;

    cumu_dailyIr_post=cumu_dailyIr_post+dailyIr_post;
    cumu_dailyIu_post=cumu_dailyIu_post+dailyIu_post;

    %%%%%%%%%%%%%%%%update S
    S=C(:,t)*ones(1,num_ens)-E-cumu_dailyIr_post-cumu_dailyIu_post;
    %%%%%%%%%%%%%%%%
    %[S,E,Ir,Iu]=checkbound(S,E,Ir,Iu,C(:,t));
    [S,E,Ir,Iu]=checkbound_yesterday(S,E,Ir,Iu,C(:,t),S_yesterday,E_yesterday,Ir_yesterday,Iu_yesterday,t);

    %%%%%%%%save stavariables for yesterday checkbound
    S_yesterday=S;
    E_yesterday=E;
    Ir_yesterday=Ir;
    Iu_yesterday=Iu;


    %%%%%%%%save posterior statevariables
    for i=1:num_loc
        for j=1:num_ens
            S_post(i,t,j)=sum(S(part(i):part(i+1)-1,j))./population(i);
            S_rec(i,t,j)=sum(S(part(i):part(i+1)-1,j)); %real numbers, not percentage
            E_post(i,t,j)=sum(E(part(i):part(i+1)-1,j))./population(i);
            E_rec(i,t,j)=sum(E(part(i):part(i+1)-1,j));   %real numbers, not percentage
            Ir_post(i,t,j)=sum(Ir(part(i):part(i+1)-1,j))./population(i);
            Ir_rec(i,t,j)=sum(Ir(part(i):part(i+1)-1,j));   %real numbers, not percentage
            Iu_post(i,t,j)= sum(Iu(part(i):part(i+1)-1,j))./population(i);
            Iu_rec(i,t,j)= sum(Iu(part(i):part(i+1)-1,j));   %real numbers, not percentage
            dailyIr_post_rec(i,j,t)=sum(dailyIr_post(part(i):part(i+1)-1,j)); %real numbers, not percentage
            dailyIu_post_rec(i,j,t)=sum(dailyIu_post(part(i):part(i+1)-1,j));  %real numbers, not percentage
            cumu_dailyIr_post_rec(i,t,j)=sum(cumu_dailyIr_post(part(i):part(i+1)-1,j)); %real numbers, not percentage
            cumu_dailyIu_post_rec(i,t,j)=sum(cumu_dailyIu_post(part(i):part(i+1)-1,j)); %real numbers, not percentage
        end
    end
    para_post(:,:,t)=para;

    % Update and display the progress bar
    fprintf('%s [%s%s] %d/%d %.2f%%\r', nickname, repmat('#', 1, round(t/num_times*20)), repmat('-', 1, 20-round(t/num_times*20)),t, num_times, t/num_times*100);

end
fprintf('\n'); % to move to the next line after the loop completes
%calculate means
dailyIr_prior_rec_mean = squeeze(mean(dailyIr_prior_rec,2));
dailyIu_prior_rec_mean = squeeze(mean(dailyIu_prior_rec,2));
dailyIr_post_rec_mean = squeeze(mean(dailyIr_post_rec,2));
dailyIu_post_rec_mean = squeeze(mean(dailyIu_post_rec,2));
cumu_dailyIr_post_rec_mean=squeeze(mean(cumu_dailyIr_post_rec,3));
cumu_dailyIu_post_rec_mean=squeeze(mean(cumu_dailyIu_post_rec,3));
E_rec_mean=mean(E_rec,3);
S_rec_mean=mean(S_rec,3);
Ir_rec_mean=mean(Ir_rec,3);
Iu_rec_mean=mean(Iu_rec,3);

cumu_dailyIr_post_mean=cumu_dailyIr_post_rec_mean./population;
cumu_dailyIu_post_mean=cumu_dailyIu_post_rec_mean./population;

save(all_file_name)
Plotting_paper(all_file_name)


