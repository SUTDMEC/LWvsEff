%%% Knapsack demonstration version 2: no dependencies between weight and
%%% efficiency

% Ïnteger weights of items
% N = 12;
% weights = randint(1,N,[1 1000])
% %%
% % Values of the items (don't have to be integers)
% values = randint(1,N,[1 100])
clear
clc
close all

tic
%Basic
% capacity=7;
% weights = [1 1 1 1 2 2 3];
% values  = -[1 1 2 3 1 3 5];
m_baseline=1711; %kg baseline vehicle weight from EPA Lotus Engineering report page 141
PP=23935; %USD MSRP 2014 Ford Fusion
FP=0.948; %USD/L
km=21561.6; %km
life=11.4; %years
BE=9; %L/100km baseline fusion fuel consumption

%maximum that the consumer is willing to pay for fuel consumption reduction
max_upfront_cost=1000;

eff=[0 0; .6 567.86; .3 210.71 ; 4 2389.29; 5.4 6417.86; 6.8 8025]; % with 40 % markup %column 1: marginal reduction in L/100km for measure, column 2: marginal cost of eff measure
%EPA study (from summary page in LW2_data_2)
% MARGINAL cost of lw=[0 0; 18.9000000000000,-6.26000000000000;68.3200000000000,-3.51000000000000;6.16000000000000,-2.48000000000000;0.530000000000000,-1.01000000000000;1.50000000000000,-0.360000000000000;16.3400000000000,-0.320000000000000;0,0;7.52000000000000,0.330000000000000;12.7000000000000,0.380000000000000;1.82000000000000,6.48000000000000;30.2500000000000,1.22000000000000;0.890000000000000,1.58000000000000;66.8300000000000,2.10000000000000;0.0800000000000000,2.45000000000000;1.07000000000000,2.79000000000000;42,3.06000000000000;2.37000000000000,3.17000000000000;2.44000000000000,3.92000000000000;32.7500000000000,5.15000000000000];
lw=[0 0;18.9000000000000,-118.314000000000;68.3200000000000,-239.803200000000;6.16000000000000,-15.2768000000000;0.530000000000000,-0.535300000000000;1.50000000000000,-0.540000000000000;16.3400000000000,-5.22880000000000;0,0;7.52000000000000,2.48160000000000;12.7000000000000,4.82600000000000;1.82000000000000,11.7936000000000;30.2500000000000,36.9050000000000;0.890000000000000,1.40620000000000;66.8300000000000,140.343000000000;0.0800000000000000,0.196000000000000;1.07000000000000,2.98530000000000;42,128.520000000000;2.37000000000000,7.51290000000000;2.44000000000000,9.56480000000000;32.7500000000000,168.662500000000];

%calculate the weights (upfront purchase costs) and values (impact of
%technology on total lifetime vehicle costs) of the lightweighting and
%efficiency technology
eff_measures={'NoEff'	'downsize/turbo'	'start-stop'	'full hybrid'	'plug-in hybrid'	'electric veh.'};
lw_measures={'NoLW','Transmission System','Body System (Group -A-) BIW & Closures','Body System (Group -D-) Glazing and Body Mechatronics','Lighting System','Driveline System','Frame and Mounting System','Fluid & Misc','Exhaust System','Fuel System','Steering System','Engine System','Electrical Dis. And Electronic Control System','Suspension System','Info, Gage, and Warning System','In-Vehicle Entertainment System','Body System (Group -B-) Interior','Body System (Group -C-) Exterior','Climate Control System','Brake System'};


dE=eff(1,1)+(0.000606*(BE-eff(1,1))+0.000708 )*0; % change in efficiency based on the change in mass selected, using the sensitivity calculated between hybrids, ICE vehicles, and EV's to 1 kg loss in weight
baseline_TCO=PP+(BE-dE)*FP*km*life/100+0+0; %the standard vehicle's total lifetime costs, without additional technology costs
wtcount=0;
%run through all of the efficiency options available
for i = 1:length(eff)

    weight(i)=floor(eff(i,2));
    if weight(i)<=1
        weight(i)=1;
        wtcount=wtcount+1;
    end
    dE_choice=eff(i,1);
    reduced_TCO(i)=PP+(BE-dE_choice)*FP*km*life/100+weight(i); 
    value(i)=baseline_TCO-reduced_TCO(i);

end
ind=i+1;
%run through the lightweighting options available
for j=1:length(lw)

    weight(ind)=floor(lw(j,2));%make the weights integer values
    %handle cases where weight is negative
    if weight(ind)<=1
        weight(ind)=1;
        wtcount=wtcount+1;
    end
    dE_choice=(0.000606*BE+0.000708 )*lw(j,1); % add the direct efficiency reduction from the powertrain to the sensitivity-adusted efficiency reduction from the weight reduction
    reduced_TCO(ind)=PP+(BE-dE_choice)*FP*km*life/100+weight(ind); 
    value(ind)=baseline_TCO-reduced_TCO(ind);     
    ind=ind+1;
end

[best, amount] = LWcase_knapsack(weight,value, max_upfront_cost);
items = find(amount);

TCO_reduction=sum(value(items));
upfront_cost=sum(weight(items));
%find the maximum reduction run in the problem
[maxval ,maxind]=max(TCO_reduction);

%compile the lw measures string
lw_str=[];
combined_measures=[eff_measures lw_measures];
for t=1:length(items)
    lw_str=[lw_str ', ' combined_measures{items(t)}  ];
end
fprintf(['\tTotal lifetime cost of ' num2str(baseline_TCO-maxval) ' USD at this upfront cost: ' num2str(upfront_cost(maxind)) '  found for with this hybrid:'   eff_measures{maxind} ' and this set of measures: '  lw_str   '\n'])
toc
keyboard

%values for different scenario plots
figure
hold
cc=cool(length(value));

plot(baseline_TCO-value);

ylabel('Total cost of ownership USD');
% title('Reduction values');
set(gca,'XTick',1:length(value))
set( gca(), 'XTickLabel',[eff_measures lw_measures ])
rotateXLabels( gca(), 45 )

%values for different scenario plots
figure
hold

plot(weight);

% title('Up-front cost values values');
ylabel('Up-front scenario cost (USD)');
set(gca,'XTick',1:length(weight))
set( gca(), 'XTickLabel', [eff_measures lw_measures ] )
rotateXLabels( gca(), 45 )

% figure
% subplot(2,1,1)
% plot(baseline_TCO-best);
% ylabel('Best total cost of ownership USD');
% axis([1 length(items) 0 max(best)+10000])
% set(gca,'XTick',1:length(items))
% set( gca(), 'XTickLabel', eff_measures )
% rotateXLabels( gca(), 45 )
% subplot(2,1,2)
% hold
% 
% stairs(amount)
% 
% axis([0 j+1 -.1 1.1]) %set the axis to 
% 
% ylabel('Technology Selected?');
% set(gca,'XTick',1:length(weights{m}))
% set(gca,'YTick',[0 1])
% set( gca(), 'XTickLabel', lw_measures )
% rotateXLabels( gca(), 45 )

