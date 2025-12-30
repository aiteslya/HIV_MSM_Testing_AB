function [] = Main_Traj_Week_create_State(varargin)
% input parameters
% par_counter: number of the parameter set to use (in the respective
% configuration file)
% batch_n: number of the batch for the current set of the parameters

% this function is the main module for calibration procedure, used to 1.
% identify valid parameter values 2. generate calibrated parameter state

% Default values
par_counter = 1;
batch_n = 1;

% Check number of inputs
if nargin == 2
    par_counter = varargin{1};
    batch_n = varargin{2};
    % Ensure inputs are integers
    if ~isscalar(par_counter) || ~isnumeric(par_counter) || par_counter ~= floor(par_counter)
        error([mfilename, ': par_counter must be an integer.']);
    end
    if ~isscalar(batch_n) || ~isnumeric(batch_n) || batch_n ~= floor(batch_n)
        error([mfilename, ': batch_n must be an integer.']);
    end
elseif nargin ~= 0
    error([mfilename, ': Function expects either 0 or 2 input arguments.']);
end

% time step: 1 week

% define handles
% population registries
Population = LargeData();
Pop_hiv = LargeData();

% snapshots of the population
Casual_Snap_time = LargeData();

Rels_steady = LargeData();
Rels_steady_hist = LargeData();

alive_ind = LargeData();
infect_alive = LargeData();
steady_dur = LargeData();
casual_dur = LargeData();

infect_diagn_time = LargeData();

total_infect_undiagn_ind = LargeData();

% initialize variables
myname = mfilename; % get own name of this script file

min_age = 15;
max_age = 75;
num_year = max_age - min_age;

band_width=10;
age_edges=min_age:band_width:max_age; % define edges of age bands, 10 years long
num_age=numel(age_edges)-1; % number of age bins

year_week = 52; % duration of a year in weeks
year_day = 365; % duration of a year in days
week = 7; % duration of a week in days
half_year = floor(year_week/2);
N_real = 200000;

% risk score deliniation
risk_edge = [0, 9, 27, Inf];

% control seed
rng(par_counter*100+batch_n); 
%%
% reading of the distributions and parameters
% folder names 
DataFolderStr='Data';
% location of configuration files
ConfFolderStr='Configuration_calibr';
ConfFileStr=['config_pars_',num2str(par_counter),'.txt']; % most of the parameter definitions in this
% file will be replaced by seperate distributions
ConfInitInfectFileStr='HIV_spec.txt';
% configuration of additional screening
AddScreenConfig='Screen00.txt';
%% Added March 3 - update the cluster
% testing rates configuration
test_file_str = ['test_rate_pars_', num2str(par_counter),'.txt'];

% file names

% Birth age propensities: based on age of sexual debut
AgeBirthPropStr='age_sexual_debut_prop.csv';

% mean probability of dying in the 
% individual age band (needed for the simulation), exponential distribution
ProbFileStr='ProbDeathAge.csv';

out_folder_str='Output';
% destination of the outputs settings
OutFolderStr=[out_folder_str,'_pars_',num2str(par_counter),'_batch_',num2str(batch_n)];
if exist(OutFolderStr,'dir')~=7
    mkdir(OutFolderStr);
end

ACSDistrFolderStr='Distributions/ACS';

% read the distribution of propensity to acquire non-steady partners
for counter_age=1:num_age
    % load distributions for individuals in a steady partnership
    file_name=['Non_SteadyACSNonStYesStPartAge',num2str(age_edges(counter_age)),'_',num2str(age_edges(counter_age+1)),'.csv'];
    file_path = fullfile(ACSDistrFolderStr, file_name);
    % Read the CSV file into a table
    cas_prop_Yes_st{counter_age} = readmatrix(file_path);
    % load distributions for individuals in a steady partnership
    file_name=['Non_SteadyACSNonStNoStPartAge',num2str(age_edges(counter_age)),'_',num2str(age_edges(counter_age+1)),'.csv'];
    file_path = fullfile(ACSDistrFolderStr, file_name);

    cas_prop_No_st{counter_age} = readmatrix(file_path);
end

% read the distribution of propensity to acquire steady partners
folder_str = 'Distributions/EMIS2017';
file_name = 'steady_parts_prop.csv';
steady_prop_age = readmatrix(fullfile(folder_str, file_name));

% Read the distribution of propensity to be 'born' at different ages
AgeBirthProp = readmatrix(fullfile(folder_str, AgeBirthPropStr));
AgeBirthProp = AgeBirthProp(:, 2);

% read death probability of dying while within an age bracket (so per whole
% 10 year bracket)
DeathProb=readmatrix(fullfile(DataFolderStr,ProbFileStr)); % needed for simulation of death, exponential distribution is hardcoded

% loading configuration files
% read the configuration file for the overall processes
cfg = cfgRead(fullfile(ConfFolderStr,ConfFileStr));
% read configuration file for initialization of infection stages
cfg_infect = cfgRead(fullfile(ConfFolderStr,ConfInitInfectFileStr));
% read configuration file for additional screening
cfg_screen = cfgRead(fullfile(ConfFolderStr,AddScreenConfig));
% read configuration file for testing
cfg_test = cfgRead(fullfile(ConfFolderStr, test_file_str));
test_rates = [ cfg_test.test_rate_1 cfg_test.test_rate_2 cfg_test.test_rate_3 cfg_test.test_rate_4];
% load distribution for testing rates, dependent on the frequency of
% changing partners
test_distr = readmatrix(fullfile('Configuration_calibr', ['test_distr_pars_',num2str(par_counter),'.csv']));

c_10 = 1/cfg.AgeBr;% ageing rate assuming 10-year age bracket, years^{-1}

T_burn_d = year_week*cfg.T_burn; % convert burn-in time of demographics and sexual networks from days to years

% read PrEP configuration file
PrEPConfFileStr = 'PrEP_cfg.txt';
cfg_PrEP = cfgRead(fullfile(ConfFolderStr,PrEPConfFileStr));

% create distribution of durations of partnerships
% steady partnerships, the distribution is in days
steady_dur_distr = create_st_dur_distr(cfg.rho/year_day);

% casual partnerships, the distribution is in days
%casual_dur_distr=create_cas_dur_distr(cfg.rho_cas/year_day);
casual_dur_distr = readmatrix('Distributions/casual_dur_distr_NS.csv');
% load age assortativity and serosorting matrices
FolderStrDistr='Distributions';

% age
% steady
file_name=fullfile(FolderStrDistr,'Steady_age.csv');
age_attr_steady=readmatrix(file_name);
% casual
file_name=fullfile(FolderStrDistr,'Casual_age.csv');
age_attr_casual=readmatrix(file_name);

% serosorting
% steady
file_name=fullfile(FolderStrDistr,'Steady_sero_mix.csv');
sero_attr_steady=readmatrix(file_name);
% casual
file_name=fullfile(FolderStrDistr,'Casual_sero_mix.csv');
sero_attr_casual=readmatrix(file_name);

% % PrEP
% file_name=fullfile(FolderStrDistr,'PrEP.csv');
% % read PrEP uptake (in years^-1)
% prep_uptake = readmatrix(file_name);
% % convert to weekly rates
% prep_uptake = prep_uptake/year_week;

% condom use
file_name = fullfile(FolderStrDistr,'condom.csv');
% read condom use (in years^-1)
% first two rows - ego during steady partnership
% last two rows - ego casual partnership
% odd rows - no PrEP
% even rows with Prep
condom  =readmatrix(file_name);

% load steady state distribution of steady partners
file_name = 'SteadyEquilGeneral.csv';
file_path = fullfile(FolderStrDistr,'EMIS2017', file_name);
steady_equil_distr = readmatrix(file_path); % no steady partner, one steady partner

% load the data for the yearly number of immigrated and diagnosed
diagn_immigr_data = readmatrix('Data/diagn_immigr_annual.csv');

% compare the running time to the available data, whether we need to
% extrapolate extra data points
final_year = (cfg.T-cfg.T_burn-0.5) + 2015+1; % we generate one extra year
if final_year>diagn_immigr_data(end,1)
    extra_years = (diagn_immigr_data(end,1)+1):1:final_year;

    % calculate the linear approximation to be used in the simulation
    diagn_immigr_LM = fitlm(diagn_immigr_data(:,1), diagn_immigr_data(:,2));
    
    % Visualize the fitted model
    temp = feval(diagn_immigr_LM, extra_years);
    diagn_immigr_FV = [diagn_immigr_data; [extra_years; temp]'];
else
    diagn_immigr_FV = diagn_immigr_data;
end

% rescale the data on diagn_immigr_FV

mult_fact = cfg.N/N_real;
diagn_immigr_FV(:,2) = mult_fact*diagn_immigr_FV(:,2);


%%

% index to access entries in population table
pop_index.id = 1;% number, unique
pop_index.birth = 2; % number (0 being the start of the simulation, can be negative)
pop_index.death = 3; % number
pop_index.age = 4; % number, positive
pop_index.sp1_id = 5; % number, id steady partner
pop_index.sp1_stdate = 6; % start of the current partnership
pop_index.nsteady = 7; % number of current steady partners
pop_index.casual_prop = 8; % propensity to acquire casual partners
pop_index.age_bin = 9; % auxiliary entry, age bin, base - 10 year age-bands
pop_index.steady_prop = 10; % propensity to acquire steady partners

% index to access entries in HIV table
hiv_index.id=1; % id, structure corresponding to the id in population table
hiv_index.status=2; % HIV status
hiv_index.date_inf=3; % date of infection
hiv_index.date_diagn=4; % date when diagnosis was received
hiv_index.date_art=5; % date when ART was started
hiv_index.date_sup=6; % date when viral suppression has been achieved
% date of expected cross to the new infection stage
hiv_index.date_new_stage = 7;
% testing rate per year
hiv_index.test_rate = 8;
hiv_index.ever_tested = 9;
hiv_index.testing_bin = 10;

% HIV status coding
HIV_st.susc=0;
HIV_st.A1=1;
HIV_st.A23=2;
HIV_st.A45=3;
HIV_st.C=4;
HIV_st.EA=5;
HIV_st.LA=6;

HIV_st.suscPrep=7;

HIV_st.A1D=8;
HIV_st.A23D=9;
HIV_st.A45D=10;
HIV_st.CD=11;
HIV_st.EAD=12;
HIV_st.LAD=13;

HIV_st.A1T=14;
HIV_st.A23T=15;
HIV_st.A45T=16;
HIV_st.CT=17;
HIV_st.EAT=18;
HIV_st.LAT=19;

HIV_st.A1S=20;
HIV_st.A23S=21;
HIV_st.A45S=22;
HIV_st.CS=23;
HIV_st.EAS=24;
HIV_st.LAS=25;

% allocate containers of the outputs
% time variable for the time series summary
% here T is the total running time in years, includes all the burning-in
% stages
% n_points: extrapolation points, ultimately, more is better, especially
% for 
tdistr=linspace(0,cfg.T*year_week,cfg.n_points);
% age distribution in the end of the run
Age_distr_ER=zeros(1, num_age);
% total population size time series
PS_TS=zeros(1 ,cfg.n_points);
% steady partnerships time series
% number of individuals without a steady partner
Steady0_TS=zeros(1,cfg.n_points);
% number of individuals with one steady partner
Steady1_TS=zeros(1,cfg.n_points); 
% duration of steady partnerships, used for calibration and de-bug
Steady_Rels_dur_Stats=zeros(1,1);

% time series which tracks number of individuals with XX casual partners in
% the last 12 months, XX=0..maximum number of partners +10

% set maximum number of casual partners
max_num_casual = max(cellfun(@(x) find(x>0, 1, 'last'), [cas_prop_No_st, cas_prop_Yes_st]));

% padding with extra number of partners to account 
% Action: 1500 was done for the debug, i think we should ultimately have +10 not
% +1500, let us remember to circle back
max_num_casual=max_num_casual+1500;

% distribution of the number of steady partners within last 12 months taken
% at the end of the simulation run, used for calibration and debug, 11
% since this is the number of partners EMIS-2017 tabulated, with 11-th
% entry denoting >=10 partners
num_steady_distr_12mon=zeros(1,11);
% number of casual partners within last 6 months, distribution 
num_casual_distr_6mon=zeros(1,max_num_casual+1);

Casual_dur_Stats=zeros(1,1);% means

% infection time series
S_TS=zeros(1,cfg.n_points);
IA1_TS=zeros(1,cfg.n_points);
IA23_TS=zeros(1,cfg.n_points);
IA45_TS=zeros(1,cfg.n_points);
IC_TS=zeros(1,cfg.n_points);
IEA_TS=zeros(1,cfg.n_points);
ILA_TS=zeros(1,cfg.n_points);
SPrep_TS=zeros(1,cfg.n_points);

IA1D_TS=zeros(1,cfg.n_points);
IA23D_TS=zeros(1,cfg.n_points);
IA45D_TS=zeros(1,cfg.n_points);
ICD_TS=zeros(1,cfg.n_points);
IEAD_TS=zeros(1,cfg.n_points);
ILAD_TS=zeros(1,cfg.n_points);

IA1T_TS=zeros(1,cfg.n_points);
IA23T_TS=zeros(1,cfg.n_points);
IA45T_TS=zeros(1,cfg.n_points);
ICT_TS=zeros(1,cfg.n_points);
IEAT_TS=zeros(1,cfg.n_points);
ILAT_TS=zeros(1,cfg.n_points);

IA1S_TS=zeros(1,cfg.n_points);
IA23S_TS=zeros(1,cfg.n_points);
IA45S_TS=zeros(1,cfg.n_points);
ICS_TS=zeros(1,cfg.n_points);
IEAS_TS=zeros(1,cfg.n_points);
ILAS_TS=zeros(1,cfg.n_points);

Dead_TS=zeros(1,cfg.n_points);
HIVDead_TS=zeros(1,cfg.n_points);

CSDead_TS=zeros(1,cfg.n_points);
CSHIVDead_TS=zeros(1,cfg.n_points);

% new infections (incidence) time series 
new_infect_TS=zeros(1,cfg.n_points);
cs_new_infect_TS=zeros(1,cfg.n_points);

% new diagnoses (incidence) time series 
new_diagn_TS=zeros(1,cfg.n_points);
cs_new_diagn_TS=zeros(1,cfg.n_points);

% proportion of diagnosed through the strategy out of the total diagnosed

prop_diagn_AHI=zeros(1,1); 

% information about sources and targets of HIV infection
% the infector
inf_source_age_distr=zeros(1,num_age);
inf_source_HIVSt_distr=zeros(1,25); % 6 stages for undiagnosed, diagnosed, treated 
% propensity to enter casual partnerships
inf_source_casual_prop_distr=zeros(1, max_num_casual + 50); % 50 is padding

% the infectee
inf_target_age_distr=zeros(1,num_age);
inf_target_HIVSt_distr=zeros(1,25);
inf_target_casual_prop_distr=zeros(1, max_num_casual + 50);

% initialize arrays that describe the total number of newly diagnosed per
% trajectory (from the start of the infection dynamics
new_total_diagn_ensemble=zeros(1,1);
new_AHI_diagn_ensemble=zeros(1,1);

%%
% parameter calculation
% calculate birth rate (year) to achieve the supplied population size
% Action: ultimately remove the code between *****

% ***************************************************

DeathRates = -c_10*log(1 - DeathProb);

% re-code AgeBirthProp into p
delta_p = AgeBirthProp(17)/44;
p(1:16) = AgeBirthProp(1:16);
p(16:num_year) = delta_p;
c = ones(1, num_year);
mu = nan(1, num_year);

for age_bin_counter = 1:6
    mu(1,(1:10)+(age_bin_counter-1)*10) = zeros(1,10)+DeathRates(age_bin_counter);
end

[lambda, N_solutions] = solve_lambda(cfg.N, mu, c, p);

N_prop = N_solutions./sum(N_solutions);

% convert rates from yearly to daily
% demographic
lambda_w=lambda/year_week;

% ***************************************************
% initally the probability was probability of dying between ages x and 
% x+n, we need this probability to be per day (so the duration of age 
% bracket in years times duration of the bracket in weeks;

DeathProbs_w=DeathProb/(cfg.AgeBr*year_week);

% details about the individuals ddiagnosed through AHI program

diagn_det_ahi=[];

% details about the people detected in a regular program
diagn_det = [];


% We set up an internal array, designated to mark time points where the data
% will twrite step is being reached, otherwise only the current state of the
% system is retained for the purposes of dynamics propagatopn

% definition of write step to record HIV time series
write_step = 2;

N_steps=ceil(1+(cfg.T*year_week-mod(cfg.T*year_week,write_step))/write_step);
twrite=0:write_step:((N_steps*write_step)-1);

% disp('The whole trajectory');
% tic

% hash maps for casual partnerships
all_casual_parts = containers.Map('KeyType', 'double', 'ValueType', 'any');
rec_casual_parts = containers.Map('KeyType', 'double', 'ValueType', 'any');

% Main simulation
% Allocating the arrays and setting the population to its initial state
% ids
% the size of the population table is allocated to be twice of the size
% predicted by the analog deterministic model

Population.Data=zeros(numel(fieldnames(pop_index)), 2*cfg.N);

% set ids of the initial population, the rest of ids will remain 0 -
% indication that respective columns contain no valid entries

Population.Data(pop_index.id,1:1:cfg.N)=1:1:cfg.N;
Id_counter=cfg.N+1; % the id of the next person that will enter the population

% date of birth, years
Population.Data(pop_index.birth,1:1:2*cfg.N)=-300;

% date of death, years, obsolete
Population.Data(pop_index.death,1:1:2*cfg.N)=-300;

% age, years
Population.Data(pop_index.age,1:1:2*cfg.N)=-1;

% snapshot of acquisition patterns, taken every 4 weeks (13 of these per
% year) throughout social network burn in phase
Casual_Snap_time.Data = nan(2*cfg.N, cfg.T_burn*13);

four_week_counter = 1;

% allocate time series containers
% population age distribution in time, size: number of time sampling points x age cathegories
yAge=zeros(N_steps,num_age);

% time series for the number of individuals with and without a steady partners
ySteadyPars=zeros(N_steps,2);

% time series for number of individuals in various HIV stages
yHIV=zeros(ceil((cfg.T*year_week-cfg.T_burn*year_week)/write_step),numel(fieldnames(HIV_st)));

% time series for total dead individuals
yDead=zeros(N_steps,1);
% time series death due to HIV infection
yHIVDead=zeros(N_steps,1);

% set ages of the individuals in the population at the start of the
% simulation
% use age distribution at the equilibrium, obtained during calibration of
% entrance rate, lambda

% Pop_age_bins=randsample(age_bins, cfg.N, true, N_prop);
Pop_age = randsample(min_age:1:(max_age-1), cfg.N, true, N_prop);

% set ages
Population.Data(pop_index.age,1:cfg.N) = Pop_age;
% assign age bands according to the age
Population.Data(pop_index.age_bin,1:cfg.N) = discretize(Pop_age,min_age:band_width:max_age);
Pop_age_bins = Population.Data(pop_index.age_bin,1:cfg.N);
% record age-distribution of the population at the start of the
% simulation
for bin_counter=1:6
    yAge(1,bin_counter)=sum(Population.Data(pop_index.age_bin,:)==bin_counter);
end

% set birth date
Population.Data(pop_index.birth,1:1:cfg.N)=-Population.Data(pop_index.age,1:1:cfg.N);

% initialize steady partnerships descriptions, prior to initiation
% everyone is single
Population.Data(pop_index.sp1_id,1:1:2*cfg.N)=-1;
Population.Data(pop_index.sp1_stdate,1:1:2*cfg.N)=-40000;

% initialize the propensity to form casual relationships
Population.Data(pop_index.casual_prop, 1:1:2*cfg.N) = -1;
% initialize the propensity to form steady partnerships
Population.Data(pop_index.steady_prop, 1:1:2*cfg.N) = -1;


max_steady_parts = 10;

for bin_counter=1:num_age
     
    ind = find(Pop_age_bins==bin_counter);
    % assign the propensities to form casual partnerships for the existing individuals
    Population.Data(pop_index.casual_prop, ind) = randsample(0:(numel(cas_prop_No_st{bin_counter})-1),numel(ind),true,cas_prop_No_st{bin_counter});

    % assign the propensities to form steady partnerships for the
    % existing individuals
    Population.Data(pop_index.steady_prop, ind) = randsample(0:max_steady_parts, numel(ind), true, steady_prop_age(bin_counter,:));
    
end    

% read transition rates from diagnosis to ART

art_rate = readmatrix('Data/art_rate.csv');
  
% array keeping track of indices of live individuals
% dynamically changing list
alive_ind.Data = 1:1:cfg.N;

% array keeping track of relationship duration 
steady_dur.Data = [];
% array keeping track of currently active steady relationships
Rels_steady.Data = [];
% array keeping track of steady relationships that were taking place within the
% last 12 months
Rels_steady_hist.Data = [];

% duration of casual partnerships
casual_dur.Data = [];

% age of infection at the time of diagnosis
infect_diagn_time.Data = [];

% initial creation of steady partnerships

tcounter=1;
% non-dimensionalize the distribution of relationships durations
steady_dur_distr=steady_dur_distr/sum(steady_dur_distr);

% initialize HIV dynamics state of the population
Pop_hiv.Data=zeros(length(fieldnames(hiv_index)), 2*cfg.N);
% initialize ids to the same values as they appear in the array
% Population
Pop_hiv.Data(hiv_index.id,:)=Population.Data(pop_index.id,:);
% Prior to seeding of infection everyone is susceptible
Pop_hiv.Data(hiv_index.status,:)=HIV_st.susc;
Pop_hiv.Data(hiv_index.date_inf,:)=-40000;
Pop_hiv.Data(hiv_index.date_diagn,:)=-40000;
Pop_hiv.Data(hiv_index.date_art,:)=-40000;
Pop_hiv.Data(hiv_index.date_sup,:)=-40000;
Pop_hiv.Data(hiv_index.date_new_stage,:)=-40000;
Pop_hiv.Data(hiv_index.test_rate,:) = 0;
Pop_hiv.Data(hiv_index.ever_tested,:) = -40000;
Pop_hiv.Data(hiv_index.testing_bin,:) = 0;


N_pairs = floor(cfg.N*(1-steady_equil_distr(1))/2);
[Population, Rels_steady, ~] = create_N_steady_pairs(Population,pop_index,tcounter,Rels_steady,alive_ind,age_attr_steady,sero_attr_steady,Pop_hiv,hiv_index,all_casual_parts, HIV_st, N_pairs);

% re-assign propensity to form casual partnerships to these who are in
% the steady partnerhsip

% find the individuals
ids_steady = unique(sort([Rels_steady.Data(:,1); Rels_steady.Data(:,2)]));

%  assign the propensity
if ~isempty(ids_steady)
    num_ids_steady = numel(ids_steady);
    inds_steady = find(ismember(Population.Data(pop_index.id,:),ids_steady'));
    % Retrieve bin_ages
    bin_ages = Population.Data(pop_index.age_bin,inds_steady);
    % Generate arrays of random numbers for each group
    rands = rand(1, num_ids_steady);

    Population.Data(pop_index.casual_prop,inds_steady) = arrayfun(@(bin_age,r) set_casual_prop(cas_prop_Yes_st{bin_age},r), bin_ages,rands);        

end

% Calculate number of pairs formed
n_single_initial = cfg.N;
n_single = sum(Population.Data(pop_index.nsteady,alive_ind.Data)==0);
num_pairs = (n_single_initial - n_single) / 2;

% Generate durations of steady partnerships that were created
j_values = randsample(1:numel(steady_dur_distr),num_pairs,true,steady_dur_distr)/week;

% Assign durations
Rels_steady.Data((end-num_pairs+1):end,3) = tcounter;
Rels_steady.Data((end-num_pairs+1):end,4) = tcounter + j_values;
Rels_steady_hist.Data=Rels_steady.Data;

% assign testing rates based on the propensity to acquire casual
% partners, this is the only instance when it will be done - everywhere
% else it will depend on the number of partners

% testing edges for the number of partners
test_edge = [0, 3, 5, 11, Inf];
parts_num = discretize(Population.Data(pop_index.casual_prop, alive_ind.Data),test_edge);

for parts_bin_counter = 1:(numel(test_edge) - 1)
    inds = alive_ind.Data(parts_num ==parts_bin_counter);
    Pop_hiv.Data(hiv_index.test_rate, inds) = randsample(test_rates, numel(inds), true, test_distr(:,parts_bin_counter));
    Pop_hiv.Data(hiv_index.testing_bin, inds) = parts_bin_counter;
end

% record number of people who have 0 and 1 steady partners at the
% start of the simulation

ySteadyPars(1,1)=n_single;
ySteadyPars(1,2)=2*num_pairs;

% non-dimensionalize the distribution of relationships durations, days
casual_dur_distr=casual_dur_distr/sum(casual_dur_distr);

burn_fl=0;
infect_alive.Data=[];
fl_collect=0;

% arrays that will contain the data about who infected whom
% source of infection schema: [id age casual_prop HIV_st_source]
infect_source_data=[];
% target of infection schema: [id age casual_prop HIV_st_source Susc_st_target]
infect_target_data=[];

% definitive log of events

max_events = 1.25*cfg.N;

infectionLog = InfectionLog(max_events);
testingLog = TestingLog(max_events);

% number of new infections
% schema: [t num]
num_new_infect=[];
num_new_infec_sum = zeros(1,3);

% number of new diagnoses
% schema: [t num]
num_new_diagn=[];

% number of tests
% schema: [time number]
total_tested = [];
total_diagn = [];

% set HIV containers prior to HIV seeding to be empty, this needs to
% happen here since Remove_List function needs them

total_infect_undiagn_ind.Data = [];% container for these who are infected but not diagnosed yet
total_diagn_ind=[];% container for these who are diagnosed but not on treatment
total_treat_ind=[];% container for these who are on treatment but not on treatment yet
total_sup_ind=[];% suppressed

% containers for infections that entered through birth/immigration
num_immigr_inf = [];
num_immigr_diagn = [];
num_immigr_suppr = [];

twrite_counter=2;
t_write_HIV_counter = 1;

% indices of individuals enrolled in PrEP programme
% to be updated in PrEP uptake, PrEP drop out, and testing
PrEP_inds.Data  = [];

% set up PrEP flag, at the start of the calibration the PrEP programme has
% not started yet
prep_fl = 0;

disp('Full trajectory run time is')
tic

for tcounter=2:1:(year_week*cfg.T)+1 % convert time from years to weeks
     % disp('Time step');
     % tic   
    % disp('Ageing');   
    % tic

    % ageing module - shift the age
    if mod(tcounter,year_week)==0
        [Population, ind_new_bracket] = Age(Population,pop_index,alive_ind);
        if numel(ind_new_bracket)>0 % update the data for individuals who moved into the new bracket

            % Define bin_ages and nsteady
            bin_ages = Population.Data(pop_index.age_bin, ind_new_bracket);
            nsteady = Population.Data(pop_index.nsteady, ind_new_bracket);
            
            % Divide ind_new_bracket into groups by nsteady
            ind_new_brackets_0 = ind_new_bracket(nsteady==0);
            ind_new_brackets_1 = ind_new_bracket(nsteady==1);
            
            % bin_ages corresponding to each group
            bin_ages_0 = bin_ages(nsteady==0);
            bin_ages_1 = bin_ages(nsteady==1);
            
            % Generate arrays of random numbers for each group
            rands_0 = rand(1, numel(ind_new_brackets_0));
            rands_1 = rand(1, numel(ind_new_brackets_1));

            old_prop = Population.Data(pop_index.casual_prop, ind_new_bracket);
            
            % Call set_casual_prop on each group with corresponding random numbers
            Population.Data(pop_index.casual_prop, ind_new_brackets_0) = arrayfun(@(bin_age, r) set_casual_prop(cas_prop_No_st{bin_age}, r), bin_ages_0, rands_0);
            Population.Data(pop_index.casual_prop, ind_new_brackets_1) = arrayfun(@(bin_age, r) set_casual_prop(cas_prop_Yes_st{bin_age}, r), bin_ages_1, rands_1);

            % re-set propensity to enter steady partnerships
            for bin_counter = 1:num_age
                inds = ind_new_bracket(Population.Data(pop_index.age_bin, ind_new_bracket) == bin_counter);
                Population.Data(pop_index.steady_prop, inds) = randsample(0:max_steady_parts, numel(inds), true, steady_prop_age(bin_counter,:));
            end

            if prep_fl & sum(old_prop>Population.Data(pop_index.casual_prop,ind_new_bracket))
                %cond = (old_prop > Population.Data(pop_index.casual_prop,ind_new_bracket)) & Pop_hiv.Data(hiv_index.status, ind_new_bracket) == HIV_st.suscPrep; 
                % new condition, checking if a person is a prep user
                % (susceptible or infected), rather than just a susceptible
                % Prep user
                cond = (old_prop > Population.Data(pop_index.casual_prop,ind_new_bracket)) & (ismember(ind_new_bracket, PrEP_inds.Data));
                if sum(cond) > 0
                    ind_change = ind_new_bracket(cond);
                    old_prop_change = old_prop(cond);
                    [Pop_hiv, inds_changed_PreP, PrEP_inds] = update_PrEP_thresh(Pop_hiv, Population, alive_ind, infect_alive, pop_index, hiv_index, HIV_st, ind_change, old_prop_change, cfg_PrEP, test_edge, PrEP_inds, test_distr, rec_casual_parts, test_rates, cfg, current_year);
                    ind_new_bracket = setdiff(ind_new_bracket, inds_changed_PreP);
                end
            end 
            
        end
    end
%        toc
% 
% %        death module 
    % disp('Death');
    % tic 
    
    [Population, Pop_hiv, fl_collect,num_died,alive_ind,steady_dur,casual_dur,Rels_steady,dead_inds,infect_alive,reset_casual_prop,total_infect_undiagn_ind,total_diagn_ind,total_treat_ind,total_sup_ind, PrEP_inds, all_casual_parts] = Death(Population,Pop_hiv,fl_collect,pop_index,hiv_index,DeathProbs_w,tcounter,alive_ind,steady_dur,casual_dur,Rels_steady,infect_alive,total_infect_undiagn_ind,total_diagn_ind,total_treat_ind,total_sup_ind,all_casual_parts, PrEP_inds); 

    if ~exist('reset_casual_prop')
        reset_casual_prop=[];
    end

    if numel(reset_casual_prop)>0
        % Define bin_ages
        bin_ages = Population.Data(pop_index.age_bin,reset_casual_prop);
        % Generate arrays of random numbers for each group
        rands = rand(1, numel(reset_casual_prop));

        old_prop = Population.Data(pop_index.casual_prop, reset_casual_prop);

        Population.Data(pop_index.casual_prop,reset_casual_prop) = arrayfun(@(bin_age,r) set_casual_prop(cas_prop_No_st{bin_age},r), bin_ages, rands);

        if prep_fl & sum(old_prop>Population.Data(pop_index.casual_prop,reset_casual_prop))
            %cond = (old_prop > Population.Data(pop_index.casual_prop,reset_casual_prop)) & Pop_hiv.Data(hiv_index.status, reset_casual_prop) == HIV_st.suscPrep; 

            % new condition, checking if a person is a prep user
            % (susceptible or infected), rather than just a susceptible
            % Prep user
            cond = (old_prop > Population.Data(pop_index.casual_prop,reset_casual_prop)) & (ismember(reset_casual_prop, PrEP_inds.Data));
            if sum(cond) > 0
                ind_change = reset_casual_prop(cond);
                old_prop_change = old_prop(cond);
                [Pop_hiv, inds_changed_PreP, PrEP_inds] = update_PrEP_thresh(Pop_hiv, Population, alive_ind, infect_alive, pop_index, hiv_index, HIV_st, ind_change, old_prop_change, cfg_PrEP, test_edge, PrEP_inds, test_distr, rec_casual_parts, test_rates, cfg, current_year);
                reset_casual_prop = setdiff(reset_casual_prop, inds_changed_PreP);
            end
        end 

    end

    assert(numel(dead_inds)==num_died);
    yDead(tcounter)=num_died;

   % disp('Birth');
   % tic
    % birth module 

    if tcounter>=2.5*year_week
        current_year = 2014 + floor((tcounter-year_week/2)/year_week);
    else
        current_year = 2014 + floor(tcounter/year_week);
    end

    % current_year
    lambda_diagn_w = 1-exp(-diagn_immigr_FV(current_year-2013, 2)/year_week);
    if current_year>=2016 && current_year<=2020
        %current_year
        theta = art_rate(2, current_year-2015);
    elseif current_year>2020
        theta = art_rate(2,end);
    else
        theta = art_rate(2, 1);
    end

    % convert rate to probability of occuring per unit time
    % switch ModelConstants.unit
    %     case 'week'
    %         theta_prob = 1 - exp(-theta/52);
    %     case 'day'
    %         theta_prob = 1 - exp(-theta/365);
    % end
    [Population,Pop_hiv,Id_counter,alive_ind,ind_new_bracket, newly_infect3, new_diagn, new_suppr, infect_alive] = Birth_Immigr_create(Population, Pop_hiv, pop_index, hiv_index,lambda_w,tcounter,Id_counter,alive_ind, cfg, lambda_diagn_w, infect_alive, HIV_st, AgeBirthProp);
    
    if ~exist('ind_new_bracket')
        ind_new_bracket=[];
    end

    num_new_born = numel(ind_new_bracket);

    if num_new_born > 0 
        % Grab age bins for new born individuals
        bin_ages = Population.Data(pop_index.age_bin,ind_new_bracket);
        % Generate random numbers for each group
        rands = rand(1, num_new_born);

        % Propensity for casual partnerships
        Population.Data(pop_index.casual_prop,ind_new_bracket) = arrayfun(@(bin_age,r) set_casual_prop(cas_prop_No_st{bin_age},r), bin_ages,rands);

        % Generate random numbers for each group
        rands = rand(1, num_new_born);
        
        % Propenstity for steady partnerships
        Population.Data(pop_index.steady_prop,ind_new_bracket) = arrayfun(@(bin_age,r) set_steady_prop(steady_prop_age(bin_age,:),r), bin_ages,rands);

        % assign testing rates, using propensities to acquire casual
        % partnerships as proxy
        parts_num = discretize(Population.Data(pop_index.casual_prop, ind_new_bracket),test_edge);

        for parts_bin_counter = 1:(numel(test_edge) - 1)
            inds = ind_new_bracket(parts_num ==parts_bin_counter);
            if numel(inds)>0
                Pop_hiv.Data(hiv_index.test_rate, inds) = randsample(test_rates, numel(inds), true, test_distr(:,parts_bin_counter));
                Pop_hiv.Data(hiv_index.testing_bin, inds) = parts_bin_counter;
            end
        end

        % write up the immigration
        if numel(newly_infect3)>0
            num_immigr_inf = [num_immigr_inf; tcounter numel(newly_infect3)];
        end

        if numel(new_diagn)>0
            num_immigr_diagn = [num_immigr_diagn; tcounter numel(new_diagn)];
        end

        if numel(new_suppr)>0
            num_immigr_suppr = [num_immigr_suppr; tcounter numel(new_suppr)];
        end
        
    end

    % toc
    % 
    % disp('Steady partnership creation');
    % tic
    num_pairs_old=size(Rels_steady.Data,1);

    % create steady partnerships
    [Population,Rels_steady]=create_steady(Population,pop_index, tcounter,Rels_steady,alive_ind, age_attr_steady,sero_attr_steady,Pop_hiv,hiv_index,all_casual_parts, HIV_st, cfg, steady_dur_distr); 
    
    num_pairs_new=size(Rels_steady.Data,1);
    if num_pairs_new>num_pairs_old % new pairs were created
        % append new pairs to the history catalogue
        Rels_steady_hist.Data=[Rels_steady_hist.Data; Rels_steady.Data((num_pairs_old+1):end,:)];

        % resample propensity to form casual partnerships
        ids_steady=(sort(unique([Rels_steady.Data((num_pairs_old+1):(num_pairs_new),1); Rels_steady.Data((num_pairs_old+1):(num_pairs_new),2)])))';
        
        if ~isempty(ids_steady)
            
            num_ids_steady = numel(ids_steady);
            inds_steady = find(ismember(Population.Data(pop_index.id,:),ids_steady'));
            % Define bin_ages
            bin_ages = Population.Data(pop_index.age_bin,inds_steady);
            % Generate arrays of random numbers for each group
            rands = rand(1, num_ids_steady);

            % handle change in propensities

            old_prop = Population.Data(pop_index.casual_prop, inds_steady);

            Population.Data(pop_index.casual_prop,inds_steady) = arrayfun(@(bin_age,r) set_casual_prop(cas_prop_Yes_st{bin_age},r), bin_ages,rands);        
    
            if prep_fl & sum(old_prop>Population.Data(pop_index.casual_prop, inds_steady))
                %cond = (old_prop > Population.Data(pop_index.casual_prop, inds_steady)) & Pop_hiv.Data(hiv_index.status, inds_steady) == HIV_st.suscPrep; 
                
                % new condition, checking if a person is a prep user
                % (susceptible or infected), rather than just a susceptible
                % Prep user
                cond = (old_prop > Population.Data(pop_index.casual_prop,inds_steady)) & (ismember(inds_steady, PrEP_inds.Data));
                
                if sum(cond) > 0
                    ind_change = inds_steady(cond);
                    old_prop_change = old_prop(cond);
                    [Pop_hiv, inds_changed_PreP, PrEP_inds] = update_PrEP_thresh(Pop_hiv, Population, alive_ind, infect_alive, pop_index, hiv_index, HIV_st, ind_change, old_prop_change, cfg_PrEP, test_edge, PrEP_inds, test_distr, rec_casual_parts, test_rates, cfg, current_year);
                    inds_steady = setdiff(inds_steady, inds_changed_PreP);
                end
            end

        end

    end
    
    % toc
    % 
    % disp('Dissolution of steady partnership');
    % tic

    % break up steady pairs
    pairs_break=find(Rels_steady.Data(:,4)<=tcounter);
    
    if numel(pairs_break)>0
        [Population, Rels_steady,steady_dur, ind_new_bracket]=break_up_steady(Population,pop_index,Rels_steady,pairs_break,steady_dur);
    end

%       resample propensity to form casual partnerships

    if numel(ind_new_bracket)>0
        bin_ages=Population.Data(pop_index.age_bin,ind_new_bracket);
        rands = rand(1, numel(ind_new_bracket));

        % handle change in propensities
        old_prop = Population.Data(pop_index.casual_prop, ind_new_bracket);

        Population.Data(pop_index.casual_prop,ind_new_bracket) = arrayfun(@(bin_age,r) set_casual_prop(cas_prop_No_st{bin_age},r), bin_ages,rands);

        if prep_fl & sum(old_prop>Population.Data(pop_index.casual_prop,ind_new_bracket))
            %cond = (old_prop > Population.Data(pop_index.casual_prop,ind_new_bracket)) & Pop_hiv.Data(hiv_index.status, ind_new_bracket) == HIV_st.suscPrep; 
            
            % new condition, checking if a person is a prep user
            % (susceptible or infected), rather than just a susceptible
            % Prep user
            cond = (old_prop > Population.Data(pop_index.casual_prop,ind_new_bracket')) & (ismember(ind_new_bracket', PrEP_inds.Data));
                
            
            if sum(cond) > 0
                ind_change = (ind_new_bracket(cond))';
                old_prop_change = old_prop(cond);
                [Pop_hiv, inds_changed_PreP, PrEP_inds] = update_PrEP_thresh(Pop_hiv, Population, alive_ind, infect_alive, pop_index, hiv_index, HIV_st, ind_change, old_prop_change, cfg_PrEP, test_edge, PrEP_inds, test_distr, rec_casual_parts, test_rates, cfg, current_year);
                
                ind_new_bracket = setdiff(ind_new_bracket, inds_changed_PreP);
            end
        end 

    end

%         toc
% % 
%         disp('Creation of casual partnerships');
%         tic

    % create casual partnerships
    % at the start of each partnership there is a sexual intercourse,
    % therefore, potentially an HIV transmission can take place
    % keep track of previously existing infections
    N_infect_old =  numel(infect_alive.Data);
            
    [Pop_hiv,infect_alive,infect_source_data,infect_target_data, all_casual_parts, rec_casual_parts, new_rels_casual, newly_infect1, infectionLog] = create_casual(Population,Pop_hiv,pop_index,tcounter,alive_ind,casual_dur_distr, age_attr_casual, sero_attr_casual, condom, burn_fl, infect_source_data, infect_target_data, infect_alive, hiv_index, cfg, HIV_st, all_casual_parts, rec_casual_parts, infectionLog);
    % Note: new_rels_casual are the partnerships for whom infection
    % event was simulated
%        toc
% % 

    N_infect_new=numel(infect_alive.Data);
    % February 27 debug point
    % if N_infect_new>N_infect_old
    %     num_new_infect=[num_new_infect; tcounter N_infect_new-N_infect_old];
    % end

%  % 
    % disp('Dissolution of casual partnerships');
    % tic
    % dissolve casual partnerships
    % the wrapper arounf break_up_casual_v3 will need to be ultimately
    % removed
    [casual_dur, all_casual_parts] = break_up_casual(casual_dur, all_casual_parts, tcounter);
%         
%      toc  
    

    if mod(tcounter,4)==0 & tcounter>half_year

        %disp('Maintanence of lists and testing rates');
        % tic
        
        % ids of these who had casual partners within the last 6
        % months

        keysList = cell2mat(keys(rec_casual_parts));

        % indices of these who had casual partners within the last 6
        % months
        [~, keysList_inds] = ismember(keysList, Population.Data(pop_index.id,alive_ind.Data));

        % end of temp

        % disregard deceased individuals
        
        validInds = keysList_inds ~= 0;
        keysList = keysList(validInds);
        keysList_inds = keysList_inds(validInds);
        
        % indices of these who did not have a casual partner within the last 6 months
        ind_no_cas_parts_6_mon = setdiff(alive_ind.Data, keysList_inds);

        % retrieve ids 
        id_no_cas_parts_6_mon = Population.Data(pop_index.id,ind_no_cas_parts_6_mon);

        % record the partner bin into the population snapshot
        Casual_Snap_time.Data(id_no_cas_parts_6_mon, four_week_counter) = 1;
 

        % find these who can be tested
        if numel(ind_no_cas_parts_6_mon)>0
            ind_can_test_no_cas_parts = ind_no_cas_parts_6_mon(Pop_hiv.Data(hiv_index.status,ind_no_cas_parts_6_mon)<HIV_st.suscPrep);
        end

        % find these whose testing bin was larger than 1 and re-set
        % their testing rate

        ind_can_test_no_cas_parts = ind_can_test_no_cas_parts(Pop_hiv.Data(hiv_index.testing_bin, ind_can_test_no_cas_parts)>1);
        % re-set testing bin
        Pop_hiv.Data(hiv_index.testing_bin, ind_can_test_no_cas_parts) = 1;

        Pop_hiv.Data(hiv_index.test_rate, ind_can_test_no_cas_parts) = randsample(test_rates, numel(ind_can_test_no_cas_parts), true, test_distr(:,1));

        num_rec_parts = nan(length(keysList), 1);

        for k = 1:length(keysList)
            currentKey = keysList(k);

            % old code
            % updated_partnerships = [];
            % 
            % for p = rec_casual_parts(currentKey)
            %     % Check if the partnership's end date is more than 182 days ago
            %     if p.end_date > (tcounter - half_year)
            %         updated_partnerships = [updated_partnerships, p];
            %     end
            % end
            % 
            % end of old code

            temp = rec_casual_parts(currentKey); 
            updated_partnerships = temp(arrayfun(@(x) x.end_date, temp)+half_year>tcounter); 
            
            rec_casual_parts(currentKey) = updated_partnerships;
            num_rec_parts(k,1) = numel(updated_partnerships);
            %  If no partnerships remain for an individual, remove them from the map
            if isempty(updated_partnerships)
                remove(rec_casual_parts, currentKey);
            end
        end
        
        % retrieve ids and convert them to indices
        %ids = keysList;
        % % ids = ids(:);
        % % [~,inds] =  ismember(ids, Population.Data(pop_index.id,:));
        % % 
        % recs_not_found = find(inds==0);
        % 
        % ids(recs_not_found,:) = [];
        % num_rec_parts(recs_not_found,:) = [];
        % inds(recs_not_found,:) = [];
%
        Casual_Snap_time.Data(keysList,four_week_counter) = discretize(num_rec_parts, risk_edge);

        % retain only these records where individuals are undiagnosed
        % and are not on PrEP
        HIV_st_remove = Pop_hiv.Data(hiv_index.status,keysList_inds)>=7;
        keysList_inds(:, HIV_st_remove) = [];
        num_rec_parts(HIV_st_remove,:) = [];
        keysList(:, HIV_st_remove) = [];
        
        % remove these whose testing bins have not changed
        % bins individuals according to the number of partners
        bin_num_parts = discretize(num_rec_parts, test_edge);

        % compare the current testing bin with previous, retain only
        % these individuals whose testing bin has changed
        same_test_bin = (bin_num_parts' - Pop_hiv.Data(hiv_index.testing_bin,keysList_inds))==0;

        keysList_inds(:,same_test_bin) = [];
        bin_num_parts(same_test_bin) = [];
        num_rec_parts(same_test_bin) = [];
        Pop_hiv.Data(hiv_index.testing_bin, keysList_inds) = bin_num_parts';
        
        % adjust testing rates for people who did have partners
      
        for bin_counter = 1:(numel(test_edge) - 1)
            % determine which 
            log_ind = bin_num_parts ==bin_counter;

            % sample
            Pop_hiv.Data(hiv_index.test_rate,keysList_inds(log_ind)) = randsample(test_rates, sum(log_ind), true, test_distr(:,bin_counter));
        end

        four_week_counter = four_week_counter + 1;

    end
% 
   

    % HIV dynamics

    % initialization of HIV dynamics, executed once
    if ~burn_fl & tcounter>=T_burn_d+1
        burn_fl=1;
        % initialize the infection
        % seed infection
        [Pop_hiv,inds_not_diagn,inds_not_treat,inds_not_sup,inds_sup] = seed_infection(Population,Pop_hiv,pop_index,hiv_index,alive_ind,cfg_infect, HIV_st, tcounter,cfg, rec_casual_parts, theta);
        infect_alive.Data=[infect_alive.Data, inds_not_diagn,inds_not_treat,inds_not_sup,inds_sup];


        total_infect_undiagn_ind.Data=inds_not_diagn;% container for these who are infected but not diagnosed yet
        total_diagn_ind=inds_not_treat;% container for these who are diagnosed but not on treatment
        total_treat_ind=inds_not_sup;% container for these who are on treatment but not on treatment yet
        total_sup_ind=inds_sup;% container for these who are suppressed

        % initialize containers for diagnosed during the main
        % simulation
        new_total_diagn_ind=[];
        new_AHI_diagn_ind=[];

        % write up as a time series the initialization point
        for HIV_St_counter=0:1:numel(fieldnames(HIV_st))-1
            yHIV(t_write_HIV_counter,HIV_St_counter+1)=sum(Pop_hiv.Data(hiv_index.status,alive_ind.Data)==HIV_St_counter);
        end

        t_write_HIV_counter = t_write_HIV_counter + 1; 
    end

    if burn_fl % burn in is over, run the HIV dynamics

        % death due to AIDS
        % find these who are in the HIV stage
        ind_AIDS1=alive_ind.Data(find(Pop_hiv.Data(hiv_index.status,alive_ind.Data)==HIV_st.EA));
        ind_AIDS2=alive_ind.Data(find(Pop_hiv.Data(hiv_index.status,alive_ind.Data)==HIV_st.LA));
        ind_AIDS3=alive_ind.Data(find(Pop_hiv.Data(hiv_index.status,alive_ind.Data)==HIV_st.EAD));
        ind_AIDS4=alive_ind.Data(find(Pop_hiv.Data(hiv_index.status,alive_ind.Data)==HIV_st.LAD));
        ind_AIDS5=alive_ind.Data(find(Pop_hiv.Data(hiv_index.status,alive_ind.Data)==HIV_st.EAT));
        ind_AIDS6=alive_ind.Data(find(Pop_hiv.Data(hiv_index.status,alive_ind.Data)==HIV_st.LAT));

        ind_AIDS=sort([ind_AIDS1,ind_AIDS2,ind_AIDS3,ind_AIDS4,ind_AIDS5,ind_AIDS6]);

        % disp('HIV death');
        % tic 

        ind_dead=[];
        if numel(ind_AIDS)>0
            [Population,fl_collect,Pop_hiv,Rels_steady,steady_dur,casual_dur,alive_ind,num_died_HIV,ind_dead,infect_alive,total_infect_undiagn_ind,total_diagn_ind,total_treat_ind,total_sup_ind,reset_casual_prop, all_casual_parts, PrEP_inds] = HIV_death(Population,fl_collect,Pop_hiv,Rels_steady,steady_dur,casual_dur,alive_ind,pop_index,hiv_index,cfg,tcounter,infect_alive,HIV_st,total_infect_undiagn_ind,total_diagn_ind,total_treat_ind,total_sup_ind, all_casual_parts, PrEP_inds);
            assert(numel(ind_dead)==num_died_HIV);
        end

        if ~exist('reset_casual_prop')
            reset_casual_prop=[];
        end

        if numel(reset_casual_prop)>0
            % Define bin_ages
            bin_ages = Population.Data(pop_index.age_bin,reset_casual_prop);
            % Generate arrays of random numbers for each group
            rands = rand(1, numel(reset_casual_prop));

            % handle change in propensities

            old_prop = Population.Data(pop_index.casual_prop, reset_casual_prop);

            Population.Data(pop_index.casual_prop,reset_casual_prop) = arrayfun(@(bin_age,r) set_casual_prop(cas_prop_No_st{bin_age},r), bin_ages, rands);

            if prep_fl & sum(old_prop>Population.Data(pop_index.casual_prop, reset_casual_prop))
                %cond = (old_prop > Population.Data(pop_index.casual_prop, reset_casual_prop)) & Pop_hiv.Data(hiv_index.status, reset_casual_prop) == HIV_st.suscPrep; 

                % new condition, checking if a person is a prep user
                % (susceptible or infected), rather than just a susceptible
                % Prep user
                cond = (old_prop > Population.Data(pop_index.casual_prop,reset_casual_prop)) & (ismember(reset_casual_prop, PrEP_inds.Data));

                if sum(cond) > 0
                    ind_change = reset_casual_prop(cond);
                    old_prop_change = old_prop(cond);
                    [Pop_hiv, inds_changed_PreP, PrEP_inds] = update_PrEP_thresh(Pop_hiv, Population, alive_ind, infect_alive, pop_index, hiv_index, HIV_st, ind_change, old_prop_change, cfg_PrEP, test_edge, PrEP_inds, test_distr, rec_casual_parts, test_rates, cfg, current_year);
                    reset_casual_prop = setdiff(reset_casual_prop, inds_changed_PreP);
                end
            end 

        end
        

        yHIVDead(tcounter)=num_died_HIV;

        if num_died_HIV>0
            % update the ids for individuals at various stages of care
            % cascade
            total_infect_undiagn_ind.Data=setdiff(total_infect_undiagn_ind.Data,ind_dead);% container for these who are infected but not diagnosed yet
            total_diagn_ind=setdiff(total_diagn_ind,ind_dead);% container for these who are diagnosed but not on treatment
            total_treat_ind=setdiff(total_treat_ind,ind_dead);% container for these who are on treatment but not on treatment yet
            total_sup_ind=setdiff(total_sup_ind,ind_dead);
            
            % added december 4 2024
            newly_infect1 = setdiff(newly_infect1, ind_dead);
            newly_infect3 = setdiff(newly_infect3, ind_dead);
            new_diagn = setdiff(new_diagn, ind_dead);
            new_suppr = setdiff(new_suppr, ind_dead);
            % end of additions done on december 4 2024
        end

        % toc 

        % disp('HIV advance');
        % tic 
        % advancement to the next stage of infection
        [Pop_hiv] = Infection_Advance(Pop_hiv,hiv_index,cfg,infect_alive,HIV_st,tcounter);
        %toc
        % 
        % disp('Infection');
        % tic 

        % infection
        infect_undiagn_old = total_infect_undiagn_ind.Data;

        N_infect_old=numel(infect_alive.Data);
        if numel(infect_alive.Data)>0
            [Pop_hiv,infect_alive,newly_infect2,infect_source_data,infect_target_data, infectionLog] = infection(Population,Pop_hiv,pop_index,hiv_index,cfg,infect_alive,tcounter,HIV_st,condom,new_rels_casual,infect_source_data,infect_target_data, all_casual_parts, Rels_steady, infectionLog,rec_casual_parts);
        end

        infect_alive.Data = [infect_alive.Data, newly_infect1, newly_infect2, newly_infect3];

        % NOTE newly infected are not added to the list of individuals
        % who can potentially be diagnosed in the same step

        N_infect_new=numel(infect_alive.Data);
        if (numel(newly_infect1)+numel(newly_infect2)+numel(newly_infect3))>0
            num_new_infect = [num_new_infect; tcounter (numel(newly_infect1)+numel(newly_infect2)+numel(newly_infect3))];
            num_new_infec_sum(1,1) = num_new_infec_sum(1,1) + numel(newly_infect1);
            num_new_infec_sum(1,2) = num_new_infec_sum(1,2) + numel(newly_infect2);
            num_new_infec_sum(1,3) = num_new_infec_sum(1,3) + numel(newly_infect3);
        end

      %toc
      if tcounter>=(cfg_infect.prep_prop_init_time)*year_week

%           % PrEP processes
%             
            if prep_fl == 0
                prep_fl = 1;
            end
            
            if numel(alive_ind.Data)>numel(infect_alive.Data) % there are individuals who are not infected
                 [Pop_hiv, PrEP_inds] = PrEPUptake_prop_thresh(Population, Pop_hiv, hiv_index, HIV_st, alive_ind, infect_alive, cfg, pop_index, PrEP_inds, cfg_PrEP, current_year);
            end
    
%            toc
%  
%           PrEP dropout
%      
           % disp('PrEP drop');
           % tic 
            
            [Pop_hiv, PrEP_inds, inds_drop] = PrEPDrop(Pop_hiv,hiv_index,HIV_st,cfg, PrEP_inds);

            if numel(inds_drop)>0 % re-set testing rates based on the number of casual partners
                ids = Population.Data(pop_index.id, inds_drop);
                id_counter = 1;
                for id = ids
                    % retrieve number of partners
                   
                    if isKey(rec_casual_parts, id)
                        num_parts = numel(rec_casual_parts(id));
                        num_parts_bin = discretize(num_parts,test_edge);
                    else
                        num_parts_bin = 1;
                    end
                    Pop_hiv.Data(hiv_index.test_rate, inds_drop(id_counter)) = randsample(test_rates,1, true,test_distr(:, num_parts_bin));
                    Pop_hiv.Data(hiv_index.testing_bin, inds_drop(id_counter)) = num_parts_bin;

                    id_counter = id_counter + 1;
                end
                                    
            end
      end
        % cART dynamics

        % disp('diagnose');
        % tic 

        % diagnose
        num_AHI_diagn_ind_old=numel(new_AHI_diagn_ind);
        if numel(infect_undiagn_old)>0 % if number of infected and undiagnosed prior to this stage is above zero
            %
            [Pop_hiv, newly_diagn,prop_diagn_AHI,new_total_diagn_ind,new_AHI_diagn_ind, PrEP_inds, newly_tested, infect_diagn_time, testingLog] = HIV_test(Pop_hiv,hiv_index,HIV_st,infect_undiagn_old,prop_diagn_AHI,tcounter,new_total_diagn_ind,new_AHI_diagn_ind,PrEP_inds, alive_ind, infect_diagn_time, testingLog, Population, rec_casual_parts, pop_index);
            
            if numel(newly_diagn)>0
                total_tested = [total_tested; tcounter numel(newly_diagn)];
            end

            if numel(newly_tested)>0
                total_diagn = [total_diagn; tcounter numel(newly_tested)];
            end

            if numel(newly_diagn)>0
                num_new_diagn=[num_new_diagn; tcounter numel(newly_diagn)];
                if num_AHI_diagn_ind_old<numel(new_AHI_diagn_ind) % record who are the individuals detected in AHI trajectory
                    % schema [date age, HIV_st before detection, recent
                    % number of casual partnerships tcounter]
                    for ind_counter=(num_AHI_diagn_ind_old+1):(numel(new_AHI_diagn_ind))
                        id=Population.Data(pop_index.id,new_AHI_diagn_ind(ind_counter));
                        age=Population.Data(pop_index.age,new_AHI_diagn_ind(ind_counter));
                        HIV_stat=Pop_hiv.Data(hiv_index.status,new_AHI_diagn_ind(ind_counter))-13;
                        
                        num_parts=numel(retrievePartners(rec_casual_parts,id));
                        
                        diagn_det_ahi=[diagn_det_ahi; id age HIV_stat num_parts tcounter];
                    end
                end
                for ind_counter=1:numel(newly_diagn)
                    ind = newly_diagn(ind_counter);
                    id=Population.Data(pop_index.id, ind);
                    age=Population.Data(pop_index.age,ind);
                    % the calculation of epidemiological status during diagnosis will
                    % depend on whether a person was detected as a part
                    % of AHI or regular screening
                    h_st = Pop_hiv.Data(hiv_index.status, ind);
                    if h_st>=HIV_st.A1T
                        HIV_stat=h_st-13;
                    else
                        HIV_stat=h_st-7;
                    end

                    num_parts=numel(retrievePartners(rec_casual_parts,id));
                    % collect the time between infection and diagnosis
                    % retrieve time of infection
                    time_infect = Pop_hiv.Data(hiv_index.date_inf,ind);
                    infect_age = tcounter-time_infect;
                    
                    diagn_det=[diagn_det; id age HIV_stat num_parts tcounter infect_age];                        
                    
                end
            end
        else
            newly_diagn=[];
        end  

        % toc 
% % % % % % 
        % disp('Start treatment');
        % tic 
        % update the indices for the infected individuals who were not
        % diagnosed yet by adding infected in this turn to the
        % difference of old undiagnosed uninfected with individuals who
        % are newly diagnosed
        total_infect_undiagn_ind.Data=[setdiff(total_infect_undiagn_ind.Data,newly_diagn), newly_infect1, newly_infect2, newly_infect3];
        
        % start treatment
        if numel(total_diagn_ind)>0
            [Pop_hiv,newly_treat,total_diagn_ind] = start_treatment(Pop_hiv,hiv_index,HIV_st, total_diagn_ind,tcounter, theta, cfg);
        else
            newly_treat=[];
        end

        % all diagnosed (before and during simulation) who are alive
        % only add newly diagnosed after starting of the treatment - it
        % is impossible that most of them will get the treatment on the
        % same day
        if ~exist("new_AHI_diagn_ind")
            new_AHI_diagn_ind=[];
        end
        if numel(new_AHI_diagn_ind)>num_AHI_diagn_ind_old % we want to remove AHI individuals as upon diagnosis they will proceed to treatment immediately
            total_diagn_ind = [total_diagn_ind,setdiff(newly_diagn,new_AHI_diagn_ind((num_AHI_diagn_ind_old+1):(numel(new_AHI_diagn_ind))))];
        else
            total_diagn_ind = [total_diagn_ind,setdiff(newly_diagn,new_AHI_diagn_ind)];
        end   

        % add newly immigrated diagnosed individuals who are not on
        % treatment
        total_diagn_ind = [total_diagn_ind, new_diagn];
%        toc
% % 
% %             % get suppressed
% %             
        % disp('Get suppressed');
        % tic 

        if numel(total_treat_ind)>0
            [Pop_hiv,newly_sup,total_treat_ind ]=get_suppressed(Pop_hiv,hiv_index,HIV_st,cfg,total_treat_ind, tcounter);
        else
            newly_sup=[];
        end


        total_treat_ind=[total_treat_ind,newly_treat];
        % add individuals who were diagnosed in this step as a part of
        % AHI strategy to the treated compartment
        if numel(new_AHI_diagn_ind)>num_AHI_diagn_ind_old
            total_treat_ind=[total_treat_ind,new_AHI_diagn_ind((num_AHI_diagn_ind_old+1):(numel(new_AHI_diagn_ind)))];
        end

       % toc 
       % 
       % disp('Drop CART');
       % tic 
       % drop out of cART

        if numel(total_sup_ind)>0
            [Pop_hiv,total_diagn_ind,total_sup_ind]=cART_drop(Pop_hiv,hiv_index,HIV_st,cfg,total_diagn_ind,total_sup_ind,tcounter);
        end

        total_sup_ind=[total_sup_ind,newly_sup];


        % add people who entered the population infected
        infect_alive.Data = [infect_alive.Data new_diagn, new_suppr];
    %toc
    end

 
%       tabulate the population, later rewrite to do it only at a handful of
%       points

    if mod(tcounter,write_step)==0
        
        % disp('Upkeep at the end of the timestep');
        % tic 

        for counter=1:1:num_age
            yAge(twrite_counter,counter)=sum(Population.Data(pop_index.age_bin,alive_ind.Data)==counter);
        end

        % tabulate relationship statistics
        % steady partnerships
        % record number of people who have 0 and 1steady partners
        for n_steady_counter=0:1
            % ySteadyPars(tcounter,n_steady_counter+1)=sum(Population.Data(pop_index.nsteady,alive_ind.Data)==n_steady_counter);
            ySteadyPars(twrite_counter,n_steady_counter+1)=sum(Population.Data(pop_index.nsteady,alive_ind.Data)==n_steady_counter);
        end

        if burn_fl
            % tabulate HIV dynamics
            for HIV_St_counter=0:1:numel(fieldnames(HIV_st))-1
                % yHIV(tcounter,HIV_St_counter+1)=sum(Pop_hiv.Data(hiv_index.status,alive_ind.Data)==HIV_St_counter);
                yHIV(t_write_HIV_counter,HIV_St_counter+1)=sum(Pop_hiv.Data(hiv_index.status,alive_ind.Data)==HIV_St_counter);
            end

            t_write_HIV_counter = t_write_HIV_counter + 1; 
        end
        %toc
        twrite_counter = twrite_counter + 1;

        
    end

    %elapsedTime(tcounter-1)=toc;
    %mean(elapsedTime)*numel(elapsedTime)
   % if mod(tcounter,365)==0
   %      toc;
   %      disp('year');
   %      tic;
   % end
end

toc

% write summary of the run
% age distribution at the end of the run
Age_distr_ER(1, :) = yAge(end,:)./sum(yAge(end,:));
% allocate the containers for the outputs collection
% total population size
PS_TS(1, :)=interp1(twrite,sum(yAge,2),tdistr);
% number of steady partnerships
Steady0_TS(1,:) = interp1(twrite,ySteadyPars(:,1),tdistr);
Steady1_TS(1,:) = interp1(twrite,ySteadyPars(:,2),tdistr);
Steady_Rels_dur_Stats(1,1) = mean(steady_dur.Data);

Steady_dur_Stats = histcounts(steady_dur.Data,'normalization','probability');
Casual_dur_Stats = histcounts(casual_dur.Data,'normalization','probability');


if ~isempty(num_new_infect)   
    if size(num_new_infect,1)>1
        % Summing infections for each unique time step
        [unique_steps, ~, idx] = unique(num_new_infect(:,1)); % Get unique time steps
        summed_infections = accumarray(idx, num_new_infect(:,2)); % Sum infections per step
        
        % Create the cleaned array
        num_new_infect_clean = [unique_steps, summed_infections];

        new_infect_TS(1,:) = interp1(num_new_infect_clean(:,1),num_new_infect_clean(:,2),tdistr);
        cs_new_infect_TS(1,:) = interp1(num_new_infect_clean(:,1),cumsum(num_new_infect_clean(:,2)),tdistr);
    else
         ind = min(find(tdistr>num_new_infect(:,1),1),numel(tdistr));
         new_infect_TS(1, ind) = num_new_infect(1,2);
         cs_new_infect_TS(1,ind) = cumsum(num_new_infect(1,2));
    end

    if tcounter>(cfg_screen.T_extra_s*year_week+10)
        inf_source_age_distr(1,:)=histcounts(infect_source_data(:,3),age_edges);
        inf_source_HIVSt_distr(1,:)=histcounts(infect_source_data(:,5),0.5:1:25.5); % 6 stages for undiagnosed, diagnosed, treated 
        inf_source_casual_prop_distr(1,:)=histcounts(infect_source_data(:,4),0:250);
    
        inf_target_age_distr(1,:)=histcounts(infect_target_data(:,3),age_edges);
        inf_target_HIVSt_distr(1,:)=histcounts(infect_target_data(:,6),0.5:1:25.5); % 6 stages for undiagnosed, diagnosed, treated 
        inf_target_casual_prop_distr(1,:)=histcounts(infect_target_data(:,4),0:250);
    end 
else
   disp('No one was infected');
end

if ~isempty(total_tested)   
    if size(total_tested,1)>1
        % Summing infections for each unique time step
        [unique_steps, ~, idx] = unique(total_tested(:,1)); % Get unique time steps
        summed_tested = accumarray(idx, total_tested(:,2)); % Sum infections per step
        
        % Create the cleaned array
        num_new_tested_clean = [unique_steps, summed_tested];

        new_tested_TS(1,:) = interp1(num_new_tested_clean(:,1),num_new_tested_clean(:,2),tdistr);
        cs_new_tested_TS(1,:) = interp1(num_new_tested_clean(:,1),cumsum(num_new_tested_clean(:,2)),tdistr);
    else
         ind = min(find(tdistr>total_tested(:,1),1),numel(tdistr));
         new_tested_TS(1, ind) = total_tested(1,2);
         cs_new_tested_TS(1,ind) = cumsum(total_tested(1,2));
    end
else
   disp('No one was tested');
end

if ~isempty(total_diagn)   
    if size(total_diagn,1)>1
        % Summing infections for each unique time step
        [unique_steps, ~, idx] = unique(total_diagn(:,1)); % Get unique time steps
        summed_diagn = accumarray(idx, total_diagn(:,2)); % Sum infections per step
        
        % Create the cleaned array
        num_new_diagn_clean = [unique_steps, summed_diagn];

        new_diagn_TS(1,:) = interp1(num_new_diagn_clean(:,1),num_new_diagn_clean(:,2),tdistr);
        cs_new_diagn_TS(1,:) = interp1(num_new_diagn_clean(:,1),cumsum(num_new_diagn_clean(:,2)),tdistr);
    else
         ind = min(find(tdistr>total_diagn(:,1),1),numel(tdistr));
         new_diagn_TS(1, ind) = total_diagn(1,2);
         cs_new_diagn_TS(1,ind) = cumsum(total_diagn(1,2));
    end
else
   disp('No one was tested');
end

% steady partnerships
ind_rels=find((tcounter-Rels_steady_hist.Data(:,4))>year_week);
Rels_steady_hist.Data(ind_rels,:)=[];

for ind=alive_ind.Data % we need to iterate through all individuals since there may be people without partnerships
    id=Population.Data(pop_index.id,ind);
    num_rels=sum(Rels_steady_hist.Data(:,1)==id)+sum(Rels_steady_hist.Data(:,2)==id);
    if num_rels<10
        num_steady_distr_12mon(1,num_rels+1)=num_steady_distr_12mon(1,num_rels+1)+1;
    else
        num_steady_distr_12mon(1,end)=num_steady_distr_12mon(1,end)+1;
    end
end
% non-dimensionalize
num_steady_distr_12mon(1,:)=num_steady_distr_12mon(1,:)./sum(num_steady_distr_12mon(1,:));

% summary of the rate of acquisition of casual partners
% casual partnerships

% clean up before the write up

keysList = keys(rec_casual_parts);

for k = 1:length(keysList)
    currentKey = keysList{k};
    updated_partnerships = [];
    
    for p = rec_casual_parts(currentKey)
        % Check if the partnership's end date is more than 182 days ago
        if p.end_date > (tcounter - year_week/2)
            updated_partnerships = [updated_partnerships, p];
        end
    end
    
    rec_casual_parts(currentKey) = updated_partnerships;
    
    %  If no partnerships remain for an individual, remove them from the map
    if isempty(updated_partnerships)
        remove(rec_casual_parts, currentKey);
    end
end

for ind=alive_ind.Data
    id=Population.Data(pop_index.id,ind);
    num_rels=numel(retrievePartners(rec_casual_parts,id));

    if num_rels+1 <= size(num_casual_distr_6mon,2)
        num_casual_distr_6mon(1, num_rels+1) = num_casual_distr_6mon(1, num_rels+1) + 1;
    else
        disp(['In simulation ',num2str(batch_n),' number of casual partners exceeds allocated size: ',  num2str(num_rels+1)]);
    end
end

% non-dimensionalize
num_casual_distr_6mon(1,:)=num_casual_distr_6mon(1,:)./sum(num_casual_distr_6mon(1,:));

% debugging figures
Tot_pop = sum(yHIV,2);

Tot_S = yHIV(:,1) + yHIV(:,8);
Tot_PrEP = yHIV(:,8);

Tot_Infect = Tot_pop - Tot_S;
Tot_Prev = Tot_Infect./cfg.N;

Tot_Infect_Undiag = sum(yHIV(:,2:7),2);
Tot_Infect_A1 = yHIV(:,2) + yHIV(:,9) + yHIV(:,15) + yHIV(:,21);
Tot_Infect_A23 = yHIV(:,3) + yHIV(:,10) + yHIV(:,16) + yHIV(:,22);
Tot_Infect_A45 = yHIV(:,4) + yHIV(:,11) + yHIV(:,17) + yHIV(:,23);
Tot_Infect_C = yHIV(:,5) + yHIV(:,12) + yHIV(:,18) + yHIV(:,24);
Tot_Infect_EA = yHIV(:,6) + yHIV(:,13) + yHIV(:,19) + yHIV(:,25);
Tot_Infect_LA = yHIV(:,7) + yHIV(:,14) + yHIV(:,20) + yHIV(:,26);

% tabulate cascade of care
Tot_Diagn = Tot_Infect - Tot_Infect_Undiag;
Tot_Diagn_per = Tot_Diagn./Tot_Infect;
Tot_Treat = Tot_Diagn - sum(yHIV(:, 9:14),2);
Tot_Treat_per = Tot_Treat./Tot_Diagn;
Tot_Sup = Tot_Treat - sum(yHIV(:, 15:20), 2);
Tot_Sup_per = Tot_Sup./Tot_Treat;

%% Debug point for steady partnership acquisition by age

age_matr = readmatrix('Distributions/EMIS2017/SteadyCurrAgeProportion.csv');
age_bin_parts = zeros(2,6);
age_bin_prop = zeros(2, 6);
data_prop = zeros(2, 6);

for age_bin = 1:6
    age_bin_parts(1,age_bin) = sum(Population.Data(pop_index.age_bin,alive_ind.Data)==age_bin & Population.Data(pop_index.nsteady,alive_ind.Data)==0);
    age_bin_parts(2,age_bin) = sum(Population.Data(pop_index.age_bin,alive_ind.Data)==age_bin & Population.Data(pop_index.nsteady,alive_ind.Data)==1);

    age_bin_prop(1,age_bin) = age_bin_parts(1,age_bin)./sum(age_bin_parts(:,age_bin));
    age_bin_prop(2,age_bin) = 1 - age_bin_prop(1,age_bin);
    
    data_prop(1,age_bin) = age_matr(1,age_bin)./sum(age_matr(:,age_bin));
    data_prop(2,age_bin) = 1 - data_prop(1,age_bin);
end


% figure(1);subplot(1,2,1);bar(age_bin_parts');
% figure(1);subplot(1,2,2);bar(age_matr');
% 
% figure(2);subplot(1,2,1);bar(age_bin_prop');
% ylim([0,0.91])
% figure(2);subplot(1,2,2);bar(data_prop');
% ylim([0,0.91]);

% serosorting analysis
sero_rels = zeros(3,1);
for row_counter = 1:length(Rels_steady.Data)
    ids = Rels_steady.Data(row_counter,1:2);
    id1 = ids(1,1);
    id2 = ids(1,2);
    ind1 = alive_ind.Data(find(Population.Data(pop_index.id,alive_ind.Data)==id1));
    ind2 = alive_ind.Data(find(Population.Data(pop_index.id,alive_ind.Data)==id2));
    st1 = Pop_hiv.Data(hiv_index.status,ind1);
    st2 = Pop_hiv.Data(hiv_index.status,ind2);
    if st1<=HIV_st.suscPrep & st2<=HIV_st.suscPrep
        sero_rels(1,1) = sero_rels(1,1) + 1;
    elseif st1>HIV_st.suscPrep & st2>HIV_st.suscPrep
        sero_rels(3,1) = sero_rels(3,1) + 1;
    else
        sero_rels(2,1) = sero_rels(2,1) + 1;
    end
end

% assortativity in casual partnerships

% Initialize a set (map) to store unique partnerships
unique_partnerships = containers.Map('KeyType', 'char', 'ValueType', 'any');

% Get all individual IDs (keys in the map)
all_ids = keys(all_casual_parts);

for i = 1:length(all_ids)
    id1 = all_ids{i};  % Individual ID
    partnerships = all_casual_parts(id1);  % Struct array of partnerships

    for j = 1:length(partnerships)
        id2 = partnerships(j).id;  % Extract partner ID
        
        % Ensure id1 < id2 for consistency
        if id1 < id2
            key = sprintf('%d_%d', id1, id2);
        else
            key = sprintf('%d_%d', id2, id1);
        end
        
        % Store in the map to ensure uniqueness
        if ~isKey(unique_partnerships, key)
            unique_partnerships(key) = sort([id1, id2]);
        end
    end
end

partnerships_cell = values(unique_partnerships);  % Get cell array of pairs
unique_relationships = vertcat(partnerships_cell{:});

% serosorting analysis
sero_cas_rels = zeros(3,1);
for row_counter = 1:length(unique_relationships)
    ids = unique_relationships(row_counter,:);
    id1 = ids(1,1);
    id2 = ids(1,2);
    ind1 = alive_ind.Data(find(Population.Data(pop_index.id,alive_ind.Data)==id1));
    ind2 = alive_ind.Data(find(Population.Data(pop_index.id,alive_ind.Data)==id2));
    st1 = Pop_hiv.Data(hiv_index.status,ind1);
    st2 = Pop_hiv.Data(hiv_index.status,ind2);
    if st1<=HIV_st.suscPrep & st2<=HIV_st.suscPrep
        sero_cas_rels(1,1) = sero_cas_rels(1,1) + 1;
    elseif st1>HIV_st.suscPrep & st2>HIV_st.suscPrep
        sero_cas_rels(3,1) = sero_cas_rels(3,1) + 1;
    else
        sero_cas_rels(2,1) = sero_cas_rels(2,1) + 1;
    end
end

% Extract partnership data
steady_partnerships = Rels_steady_hist.Data;  % [id1, id2, start_time, end_time]

% Filter partnerships that started at time 52 or later
valid_partnerships = steady_partnerships(steady_partnerships(:,4) >= 52, :);

% Get all alive individuals at current time (105)
alive_ids = Population.Data(pop_index.id, alive_ind.Data);  % Extract alive individual IDs

% Initialize a dictionary where keys are alive individual IDs, values are partner counts
num_partners = containers.Map('KeyType', 'double', 'ValueType', 'double');

% Initialize all alive individuals with 0 partners
for i = 1:length(alive_ids)
    num_partners(alive_ids(i)) = 0;
end

% Count the number of partners per individual
for i = 1:size(valid_partnerships,1)
    id1 = valid_partnerships(i,1);
    id2 = valid_partnerships(i,2);

    % Increment partner count only for individuals who are alive at time 105
    if isKey(num_partners, id1)
        num_partners(id1) = num_partners(id1) + 1;
    end
    if isKey(num_partners, id2)
        num_partners(id2) = num_partners(id2) + 1;
    end
    
end

% Extract the number of partners per alive individual
partner_counts = zeros(length(alive_ids),1);
age_bins = zeros(length(alive_ids),1);

for i = 1:length(alive_ids)
    id = alive_ids(i);
    
    % Retrieve age bin for the individual
    age_bin = Population.Data(pop_index.age_bin, Population.Data(pop_index.id,:) == id);
    age_bins(i) = age_bin;

    % Retrieve partner count (safe since all alive individuals are already initialized)
    partner_counts(i) = num_partners(id);
end

% Compute distributions per age bin
max_partners = max(partner_counts); % Find max number of partners
age_bin_counts = zeros(6, max_partners + 1);  % Rows: age bins, Columns: partner counts

for i = 1:length(alive_ids)
    age_bin = age_bins(i);
    partners = partner_counts(i);
    age_bin_counts(age_bin, partners + 1) = age_bin_counts(age_bin, partners + 1) + 1;
end

% Normalize the distribution for each age bin
for bin = 1:6
    total_individuals = sum(age_bin_counts(bin, :));
    if total_individuals > 0
        age_bin_counts(bin, :) = age_bin_counts(bin, :) / total_individuals;
    end
end

% Display the normalized distribution
disp('Normalized distribution of the number of partners per age bin:');
disp(age_bin_counts);


%%


% Record outputs into CSV files
% population
%equilibrium distribution at the end of each run per each trajectory
writematrix(Age_distr_ER,[OutFolderStr,'/Age_distr_ER.csv']);
%population size time series, one per trajectory
writematrix(PS_TS,[OutFolderStr,'/PS_TS.csv']);
%steady partnerships
%time series for proportion of single individuals, one per trajectory
writematrix(Steady0_TS,[OutFolderStr,'/Steady0_TS.csv']);
%time series for proportion of individuals in a partnership, one per trajectory
writematrix(Steady1_TS,[OutFolderStr,'/Steady1_TS.csv']);
%average duration, one per trajectory
writematrix(Steady_Rels_dur_Stats,[OutFolderStr,'/Rels_dur_Stats.csv']);
%distribution of the number of steady partners within 12 months
writematrix(num_steady_distr_12mon,[OutFolderStr,'/num_steady_distr_6mon.csv']);

writematrix(Casual_dur_Stats,[OutFolderStr,'/Casual_dur_Stats.csv']);
%distribution of the number of steady partners within 6 months
writematrix(num_casual_distr_6mon,[OutFolderStr,'/num_casual_distr_6mon.csv']);


% epidemiological process
writematrix(S_TS,[OutFolderStr,'/S_TS.csv']);
writematrix(IA1_TS,[OutFolderStr,'/IA1_TS.csv']);
writematrix(IA23_TS,[OutFolderStr,'/IA23_TS.csv']);
writematrix(IA45_TS,[OutFolderStr,'/IA45_TS.csv']);
writematrix(IC_TS,[OutFolderStr,'/IC_TS.csv']);
writematrix(IEA_TS,[OutFolderStr,'/IEA_TS.csv']);
writematrix(ILA_TS,[OutFolderStr,'/ILA_TS.csv']);
writematrix(SPrep_TS,[OutFolderStr,'/SPrep_TS.csv']);

writematrix(IA1D_TS,[OutFolderStr,'/IA1D_TS.csv']);
writematrix(IA23D_TS,[OutFolderStr,'/IA23D_TS.csv']);
writematrix(IA45D_TS,[OutFolderStr,'/IA45D_TS.csv']);
writematrix(ICD_TS,[OutFolderStr,'/ICD_TS.csv']);
writematrix(IEAD_TS,[OutFolderStr,'/IEAD_TS.csv']);
writematrix(ILAD_TS,[OutFolderStr,'/ILAD_TS.csv']);

writematrix(IA1T_TS,[OutFolderStr,'/IA1T_TS.csv']);
writematrix(IA23T_TS,[OutFolderStr,'/IA23T_TS.csv']);
writematrix(IA45T_TS,[OutFolderStr,'/IA45T_TS.csv']);
writematrix(ICT_TS,[OutFolderStr,'/ICT_TS.csv']);
writematrix(IEAT_TS,[OutFolderStr,'/IEAT_TS.csv']);
writematrix(ILAT_TS,[OutFolderStr,'/ILAT_TS.csv']);

writematrix(IA1S_TS,[OutFolderStr,'/IA1S_TS.csv']);
writematrix(IA23S_TS,[OutFolderStr,'/IA23S_TS.csv']);
writematrix(IA45S_TS,[OutFolderStr,'/IA45S_TS.csv']);
writematrix(ICS_TS,[OutFolderStr,'/ICS_TS.csv']);
writematrix(IEAS_TS,[OutFolderStr,'/IEAS_TS.csv']);
writematrix(ILAS_TS,[OutFolderStr,'/ILAS_TS.csv']);

writematrix(Dead_TS,[OutFolderStr,'/Dead_TS.csv']);
writematrix(HIVDead_TS,[OutFolderStr,'/HIVDead_TS.csv']);
writematrix(CSDead_TS,[OutFolderStr,'/CSDead_TS.csv']);
writematrix(CSHIVDead_TS,[OutFolderStr,'/CSHIVDead_TS.csv']);

writematrix(new_infect_TS,[OutFolderStr,'/new_infect_TS.csv']);
writematrix(cs_new_infect_TS,[OutFolderStr,'/cs_new_infect_TS.csv']);

writematrix(new_diagn_TS,[OutFolderStr,'/new_diagn_TS.csv']);
writematrix(cs_new_diagn_TS,[OutFolderStr,'/cs_new_diagn_TS.csv']);

writematrix(new_tested_TS,[OutFolderStr,'/new_tested_TS.csv']);
writematrix(cs_new_tested_TS,[OutFolderStr,'/cs_new_tested_TS.csv']);

% who infects who
if tcounter>(cfg_screen.T_extra_s*year_week+10)
    writematrix(inf_source_age_distr,[OutFolderStr,'/inf_source_age_distr.csv']);
    writematrix(inf_source_HIVSt_distr,[OutFolderStr,'/inf_source_HIVSt_distr.csv']);
    writematrix(inf_source_casual_prop_distr,[OutFolderStr,'/inf_source_casual_prop_distr.csv']);

    writematrix(inf_target_age_distr,[OutFolderStr,'/inf_target_age_distr.csv']);
    writematrix(inf_target_HIVSt_distr,[OutFolderStr,'/inf_target_HIVSt_distr.csv']);
    writematrix(inf_target_casual_prop_distr,[OutFolderStr,'/inf_target_casual_prop_distr.csv']);
end

%write durations of steady and casual partnerships, from a single run
writematrix(Steady_dur_Stats,[OutFolderStr,'/Steady_dur_Stats.csv']);
writematrix(Casual_dur_Stats,[OutFolderStr,'/Casual_dur_Stats.csv']);

% number of new total diagnosed
writematrix(new_total_diagn_ensemble,[OutFolderStr,'/new_total_diagn_ensemble.csv']);
writematrix(new_AHI_diagn_ensemble,[OutFolderStr,'/new_AHI_diagn_ensemble.csv']);

%number of AHI within the trajectory

% proportion of AHI diagnosed out of all diagnosis
writematrix(prop_diagn_AHI,[OutFolderStr,'/prop_diagn_AHI.csv']);

if ~isempty(diagn_det_ahi)
    % details of individuals detected through AHI trajectory
    % schema [id age HIV at the time of detection recent number of casual parts]
    writematrix(diagn_det_ahi(:,2),[OutFolderStr,'/Det_diagn_AHI_Age.csv']);
    writematrix(diagn_det_ahi(:,3),[OutFolderStr,'/Det_diagn_AHI_HIV_st.csv']);
    writematrix(diagn_det_ahi(:,4),[OutFolderStr,'/Det_diagn_AHI_Num_Cas.csv']);
    writematrix(diagn_det_ahi(:,5),[OutFolderStr,'/Det_diagn_AHI_Time_det.csv']);
end

if ~isempty(diagn_det)
    % details of individuals detected through AHI trajectory
    % schema [id age HIV at the time of detection recent number of casual parts]
    writematrix(diagn_det(:,2),[OutFolderStr,'/Det_diagn_Age.csv']);
    writematrix(diagn_det(:,3),[OutFolderStr,'/Det_diagn_HIV_st.csv']);
    writematrix(diagn_det(:,4),[OutFolderStr,'/Det_diagn_Num_Cas.csv']);
    writematrix(diagn_det(:,5),[OutFolderStr,'/Det_diagn_Time_det.csv']);
    writematrix(diagn_det(:,6),[OutFolderStr,'/infect_age.csv']);
end

writematrix(infect_diagn_time.Data, fullfile(OutFolderStr, 'age_infect_at_diagn.csv'));

% write HIV dynamics for calibration with the time series
% record outputs and time array corresponding to them (convenience)

% % time array
% writematrix(t_HIV_write, fullfile(OutFolderStr,'t_HIV_write.csv'));

% total population
writematrix(Tot_pop, fullfile(OutFolderStr,'Tot_pop.csv'));

% susceptible
writematrix(Tot_S, fullfile(OutFolderStr,'Tot_PrEP.csv'));
writematrix(Tot_PrEP, fullfile(OutFolderStr,'Tot_PrEP.csv'));

% infected
writematrix(Tot_Infect, fullfile(OutFolderStr,'Tot_Infect.csv'));
writematrix(Tot_Prev, fullfile(OutFolderStr,'Tot_Prev.csv'));
writematrix(Tot_Infect_Undiag, fullfile(OutFolderStr,'Tot_Infect_Undiag.csv'));

writematrix(Tot_Infect_A1, fullfile(OutFolderStr,'Tot_Infect_A1.csv'));
writematrix(Tot_Infect_A23, fullfile(OutFolderStr,'Tot_Infect_A23.csv'));
writematrix(Tot_Infect_A45, fullfile(OutFolderStr,'Tot_Infect_A45.csv'));    
writematrix(Tot_Infect_C, fullfile(OutFolderStr,'Tot_Infect_C.csv'));        
writematrix(Tot_Infect_EA, fullfile(OutFolderStr,'Tot_Infect_EA.csv'));            
writematrix(Tot_Infect_LA, fullfile(OutFolderStr,'Tot_Infect_LA.csv'));   

writematrix(Tot_Diagn, fullfile(OutFolderStr,'Tot_Diagn.csv'));
writematrix(Tot_Diagn_per, fullfile(OutFolderStr,'Tot_Diagn_per.csv'));

writematrix(Tot_Treat, fullfile(OutFolderStr,'Tot_Treat.csv'));
writematrix(Tot_Treat_per, fullfile(OutFolderStr,'Tot_Treat_per.csv'));

writematrix(Tot_Sup, fullfile(OutFolderStr,'Tot_Sup.csv'));
writematrix(Tot_Sup_per, fullfile(OutFolderStr,'Tot_Sup_per.csv'));

% % write age of infection
% if ~isempty(num_new_infect)
%     writematrix(infect_source_data(:,1), fullfile(OutFolderStr,'infect_source_time.csv'));
%     writematrix(infect_source_data(:,5), fullfile(OutFolderStr,'infect_source_source.csv'));
%     writematrix(infect_source_data(:,6), fullfile(OutFolderStr,'infect_source_age_infect.csv'));
% end

infectionLog.trimData();
infectionLog.saveToFile(fullfile(OutFolderStr,'infection_log.mat'));

testingLog.trimData();
testingLog.saveToFile(fullfile(OutFolderStr,'testing_log.mat'));

% preserve state
if ~exist('State', 'dir')
    mkdir('State');
end

State_str = ['State/State_par_',num2str(par_counter),'_batch_',num2str(batch_n)];
if ~exist(State_str, 'dir')
    mkdir(State_str);
end

% Path and file name for saving the object
filename = fullfile(State_str, 'Population.mat');

% Save the Population object to file
save(filename, 'Population', '-v7.3'); % Using '-v7.3' to handle larger data efficiently

% Path and file name for saving the object
filename = fullfile(State_str, 'Pop_hiv.mat');

% Save the Pop_hiv object to file
save(filename, 'Pop_hiv', '-v7.3'); % Using '-v7.3' to handle larger data efficiently

% Path and file name for saving the object
filename = fullfile(State_str, 'Rels_steady.mat');

% Save the Rels_steady object to file
save(filename, 'Rels_steady', '-v7.3'); % Using '-v7.3' to handle larger data efficiently

% Path and file name for saving the object
filename = fullfile(State_str, 'Rels_steady_hist.mat');

% Save the Rels_steady_hist object to file
save(filename, 'Rels_steady_hist', '-v7.3'); % Using '-v7.3' to handle larger data efficiently

% Path and file name for saving the object
filename = fullfile(State_str, 'alive_ind.mat');

% Save the alive_ind object to file
save(filename, 'alive_ind', '-v7.3'); % Using '-v7.3' to handle larger data efficiently

% Path and file name for saving the object
filename = fullfile(State_str, 'infect_alive.mat');

% Save the infect_alive object to file
save(filename, 'infect_alive', '-v7.3'); % Using '-v7.3' to handle larger data efficiently

%
filename = fullfile(State_str, 'all_casual_parts.mat');
save(filename, 'all_casual_parts');
filename = fullfile(State_str, 'rec_casual_parts.mat');
save(filename, 'rec_casual_parts');

% Path and file name for saving the object
filename = fullfile(State_str, 'inds_not_diagn.csv');

% Save the inds_not_diagn object to file
writematrix(inds_not_diagn, filename);

% Path and file name for saving the object
filename = fullfile(State_str, 'inds_not_sup.csv');

% Save the inds_not_sup object to file
writematrix(inds_not_sup, filename);

% Path and file name for saving the object
filename = fullfile(State_str, 'inds_not_treat.csv');

% Save the inds_not_treat object to file
writematrix(inds_not_treat, filename);

% Path and file name for saving the object
filename = fullfile(State_str, 'inds_sup.csv');

% Save the inds_not_treat object to file
writematrix(inds_sup, filename);

% Path and file name for saving the object
filename = fullfile(State_str, 'total_infect_undiagn_ind.mat');

% Save the total_infect_undiagn_ind object to file
save(filename, 'total_infect_undiagn_ind', '-v7.3');
%writematrix(total_infect_undiagn_ind, filename);

% Path and file name for saving the object
filename = fullfile(State_str, 'total_sup_ind.csv');

% Save the total_sup_ind object to file
writematrix(total_sup_ind, filename);

% Path and file name for saving the object
filename = fullfile(State_str, 'total_diagn_ind.csv');

% Save the total_diagn_ind object to file
writematrix(total_diagn_ind, filename);

% Path and file name for saving the object
filename = fullfile(State_str, 'total_treat_ind.csv');

% Save the total_treat_ind object to file
writematrix(total_treat_ind, filename);

% Save the current state of the random number generator
rngState = rng;

% Define the file path for saving the RNG state
rngStateFilePath = fullfile(State_str, 'rng_state.mat');

% Save the state
save(rngStateFilePath, 'rngState');

% Path and file name for saving the object
filename = fullfile(State_str, 'PrEP_inds.mat');

% Save the alive_ind object to file
save(filename, 'PrEP_inds', '-v7.3'); % Using '-v7.3' to handle larger data efficiently


end