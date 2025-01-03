%para: Z,D,mu,theta,alpha1,alpha2,...,alpha3142,beta1,...,beta3142

load population.mat
load dailyincidence_real.csv
load statecodes.mat

num_sims=1;

parafit=zeros(6,100,num_sims);



num_loc=size(population,1);
Zlow=2;Zup=5;%latency period
Dlow=2;Dup=5;%infectious period
mulow=0.2;muup=0.45;%relative transmissibility
thetalow=0;thetaup=0.02;%movement factor
alphalow=0.025;alphaup=0.60;%reporting rate
betalow=0.2;betaup=4;%transmission rate


paramin=[Zlow;Dlow;mulow;thetalow;ones(num_loc,1)*alphalow;ones(num_loc,1)*betalow];
paramax=[Zup;Dup;muup;thetaup;ones(num_loc,1)*alphaup;ones(num_loc,1)*betaup];

ca_infection_induced_seroprevalence = 0.019 ;  % July https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2789086
days_to_consider=6*30; 
weighted_alphaca=sum(dailyincidence_real(:,1:days_to_consider),2)./(population*ca_infection_induced_seroprevalence); %percentage of people that got reported
statecodes.Var4 = weighted_alphaca;
statecodeca = statecodes.Var4(strcmp(statecodes.Var3, 'CA'));
median(statecodeca);
%ans =
%    0.0269
statecodes(:,4) = [];

us_infection_induced_seroprevalence = 0.035; % July 2020  https://jamanetwork.com/journals/jama/fullarticle/2784013#joi210100f3
days_to_consider=6*30; 
weighted_alphaus=sum(dailyincidence_real(:,1:days_to_consider),2)./(population*us_infection_induced_seroprevalence); %percentage of people that got reported
statecodes.Var4 = weighted_alphaus;
statecodesus = statecodes.Var4(strcmp(statecodes.Var3, 'US'));
median(statecodesus);  %median of the state alphas in us based on the infection induce seroprevalence
% ans =
%     0.2563
statecodes(:,4) = [];

mx_infection_induced_seroprevalence = 0.035; % february 2020  10.3390/microorganisms9040850
days_to_consider=38; 
weighted_alphamx=sum(dailyincidence_real(:,1:days_to_consider),2)./(population*mx_infection_induced_seroprevalence); %percentage of people that got reported
statecodes.Var4 = weighted_alphamx;
statecodemx = statecodes.Var4(strcmp(statecodes.Var3, 'MX'));
median(statecodemx);
% ans =
%     0
statecodes(:,4) = [];
% CA and MX alpha calcualted from infection induced seroprevalence is too low
% the prior distribution will be assign in the initialization step in the model

US_weighted_alpha=median(statecodesus);

for n_sim=1:num_sims
    parafit(1,:,n_sim)= squeeze(Zlow + (Zup-Zlow)*rand(100,1));
    parafit(2,:,n_sim)= squeeze(Dlow + (Dup-Dlow)*rand(100,1));
    parafit(3,:,n_sim)= squeeze(mulow + (muup-mulow)*rand(100,1));
    parafit(4,:,n_sim)= squeeze(thetalow + (thetaup-thetalow)*rand(100,1));
    parafit(5,:,n_sim)= squeeze(alphalow + (alphaup-alphalow)*rand(100,1)); % only for alpha US

    parafit(6,:,n_sim)=squeeze(normrnd((betaup-betalow)/2, 0.75, [1, 100])); %normally dist
    parafit(6,(parafit(6,:,n_sim) < betalow),n_sim) = betalow+(betaup-betalow)*rand(); %out of bound values back inside randomly
    parafit(6,(parafit(6,:,n_sim) > betaup),n_sim) = betalow+(betaup-betalow)*rand(); %out of bound values back inside randomly

end

parafit=squeeze(parafit);


clear statecodes dailyincidence_real days_to_consider rowsToUpdate population n_sim statecodesus statecodesmx statecodeca ...
    ca_infection_induced_seroprevalence mx_infection_induced_seroprevalence us_infection_induced_seroprevalence ...
    weighted_alphamx weighted_alphaus weighted_alphaca



save parafit_rev.mat