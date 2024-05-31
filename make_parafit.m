% parafit is the matrix used to initialized the parameters:
%para: Z,D,mu,theta,alpha1,alpha2,...,alpha3142,beta1,...,beta3142

load population.mat
load dailyincidence_real.csv
load statecodes.mat

num_sims=1; %number of simulation/dimentions for parafit

parafit=zeros(6,100,num_sims);

num_loc=size(population,1);
%parameter lower and upper bounds
Zlow=2;Zup=5; %latency period
Dlow=2;Dup=5; %infectious period
mulow=0.2;muup=0.45; %relative transmissibility
thetalow=0;thetaup=0.02; %movement factor
alphalow=0.025;alphaup=0.60; %reporting rate
betalow=0.2;betaup=4; %transmission rate


paramin=[Zlow;Dlow;mulow;thetalow;ones(num_loc,1)*alphalow;ones(num_loc,1)*betalow];
paramax=[Zup;Dup;muup;thetaup;ones(num_loc,1)*alphaup;ones(num_loc,1)*betaup];

us_infection_induced_seroprevalence = 0.035; %July 2020  https://jamanetwork.com/journals/jama/fullarticle/2784013#joi210100f3
days_to_consider=6*30; %first six months of pandemic
weighted_alphaus=sum(dailyincidence_real(:,1:days_to_consider),2)./(population*us_infection_induced_seroprevalence); %percentage of people that got reported
statecodes.Var4 = weighted_alphaus;
statecodesus = statecodes.Var4(strcmp(statecodes.Var3, 'US'));
median(statecodesus)  %median of the state alphas in us based on the infection induce seroprevalence


% Manually changing the content of the 4th column to 0.08 for the selected
% rows (wighted alpha)
rowsToUpdate = strcmp(statecodes.Var3, 'MX');
statecodes.Var4(rowsToUpdate) = 0.08;
rowsToUpdate = strcmp(statecodes.Var3, 'CA');
statecodes.Var4(rowsToUpdate) = 0.12;
rowsToUpdate = strcmp(statecodes.Var3, 'US');
statecodes.Var4(rowsToUpdate) = 0.25;

weighted_alpha=statecodes.Var4;

for n_sim=1:num_sims
    parafit(1,:,n_sim)= squeeze(Zlow + (Zup-Zlow)*rand(100,1));
    parafit(2,:,n_sim)= squeeze(Dlow + (Dup-Dlow)*rand(100,1));
    parafit(3,:,n_sim)= squeeze(mulow + (muup-mulow)*rand(100,1));
    parafit(4,:,n_sim)= squeeze(thetalow + (thetaup-thetalow)*rand(100,1));
    parafit(5,:,n_sim)= squeeze(alphalow + (alphaup-alphalow)*rand(100,1));

    parafit(6,:,n_sim)=squeeze(normrnd((betaup-betalow)/2, 0.75, [1, 100])); %normally dist
    parafit(6,(parafit(6,:,n_sim) < betalow),n_sim) = betalow+(betaup-betalow)*rand(); %out of bound values back inside randomly
    parafit(6,(parafit(6,:,n_sim) > betaup),n_sim) = betalow+(betaup-betalow)*rand(); %out of bound values back inside randomly

end

parafit=squeeze(parafit);
clear statecodes dailyincidence_real days_to_consider rowsToUpdate population n_sim
save parafit_b-rnd_higher-alow.mat