%%% Knapsack demonstration
% Ïnteger weights of items
% N = 12;
% weights = randint(1,N,[1 1000])
% %%
% % Values of the items (don't have to be integers)
% values = randint(1,N,[1 100])
try
% clear
% clc
close all

tic
%Basic
% capacity=7;
% weights = [1 1 1 1 2 2 3];
% values  = -[1 1 2 3 1 3 5];
m_baseline=1711; %kg baseline vehicle weight from EPA Lotus Engineering report page 141
PP=23935; %USD MSRP 2014 Ford Fusion
FP=0.95; %USD/L
km=21562; %km
life=11.4; %years
BE=9; %L/100km baseline fusion fuel consumption

max_upfront_vector=[10000 5000 2000 1000 750 500 400 300 200 100 50 1]; %maximum up-front cost

plot_scen=0; %flag to plot scenarios

eff=[0 0; .6 567.86; .3 210.71 ; 4 2389.29; 5.4 6417.86; 6.8 8025]; % with 40 % markup %column 1: marginal reduction in L/100km for measure, column 2: marginal cost of eff measure
%EPA study (from summary page in LW2_data_2)
% MARGINAL cost of lw=[0 0; 18.9000000000000,-6.26000000000000;68.3200000000000,-3.51000000000000;6.16000000000000,-2.48000000000000;0.530000000000000,-1.01000000000000;1.50000000000000,-0.360000000000000;16.3400000000000,-0.320000000000000;0,0;7.52000000000000,0.330000000000000;12.7000000000000,0.380000000000000;1.82000000000000,6.48000000000000;30.2500000000000,1.22000000000000;0.890000000000000,1.58000000000000;66.8300000000000,2.10000000000000;0.0800000000000000,2.45000000000000;1.07000000000000,2.79000000000000;42,3.06000000000000;2.37000000000000,3.17000000000000;2.44000000000000,3.92000000000000;32.7500000000000,5.15000000000000];
lw=[0 0;18.9000000000000,-118.314000000000;68.3200000000000,-239.803200000000;6.16000000000000,-15.2768000000000;0.530000000000000,-0.535300000000000;1.50000000000000,-0.540000000000000;16.3400000000000,-5.22880000000000;0,0;7.52000000000000,2.48160000000000;12.7000000000000,4.82600000000000;1.82000000000000,11.7936000000000;30.2500000000000,36.9050000000000;0.890000000000000,1.40620000000000;66.8300000000000,140.343000000000;0.0800000000000000,0.196000000000000;1.07000000000000,2.98530000000000;42,128.520000000000;2.37000000000000,7.51290000000000;2.44000000000000,9.56480000000000;32.7500000000000,168.662500000000];

%calculate the weights (upfront purchase costs) and values (impact of
%technology on total lifetime vehicle costs) of the lightweighting and
%efficiency technology
eff_measures={'NoEff'	'downsize/turbo'	'start-stop'	'full hybrid'	'plug-in hybrid'	'electric veh.'};
lw_measures={'NoLW','Transmission System','Body System (Group -A-) BIW & Closures','Body System (Group -D-) Glazing and Body Mechatronics','Lighting System','Driveline System','Frame and Mounting System','Fluid & Misc','Exhaust System','Fuel System','Steering System','Engine System','Electrical Dis. And Electronic Control System','Suspension System','Info, Gage, and Warning System','In-Vehicle Entertainment System','Body System (Group -B-) Interior','Body System (Group -C-) Exterior','Climate Control System','Brake System'};
eff_measures_short={'None'	'DS/T'	'S-S'	'Full'	'Plug-in'	'EV'};

dE=eff(1,1)+(0.000606*(BE-eff(1,1))+0.000708 )*0; % change in efficiency based on the change in mass selected, using the sensitivity calculated between hybrids, ICE vehicles, and EV's to 1 kg loss in weight
baseline_TCO=PP+(BE-dE)*FP*km*life/100+0+0; %the standard vehicle's total lifetime costs, without additional technology costs

%run through all of the maximum up-front costs
for s=1:length(max_upfront_vector)
    %run through all of the efficiency options available
    for i = 1:length(eff)
        % set the maximum that the consumer is willing to pay for fuel consumption reduction
        max_upfront_cost=max_upfront_vector(s);


        for j=1:length(lw)
            weight(j)=ceil(lw(j,2));%make the weights integer values
            wt_raw(j)=lw(j,2); %store the un-rounded up-front cost
            %handle cases where weight is negative
            dE_choice=(0.000606*(BE-eff(i,1))+0.000708 )*lw(j,1); % add the direct efficiency reduction from the powertrain to the sensitivity-adusted efficiency reduction from the weight reduction
            reduced_TCO(j)=PP+(BE-dE_choice)*FP*km*life/100+wt_raw(j); 
            value(j)=baseline_TCO-reduced_TCO(j);     
        end
        %run the efficiency choice
        weight(j+1)=ceil(eff(i,2));
        wt_raw(j+1)=eff(i,2);
        dE_choice=eff(i,1); % add the direct efficiency reduction from the powertrain to the sensitivity-adusted efficiency reduction from the weight reduction
        reduced_TCO(j+1)=PP+(BE-dE_choice)*FP*km*life/100+wt_raw(j+1); 
        value(j+1)=baseline_TCO-reduced_TCO(j+1);   

        red_TCO{i}=reduced_TCO;
        raw_values{i}=value;
        %LW problem 
        weights_raw{i}=weight; %store the rounded weights before removal of negative values
        
        measures{i}=[lw_measures eff_measures(i)];
        
        %sort by wt, lowest to highest
        [weight, ind_sort]=sort(weight);
        value=value(ind_sort);
        measures{i}=measures{i}(ind_sort);
        wt_unrounded{i}=wt_raw(ind_sort); %store the unrounded true weights in sorted order


        ind_selected=find(weight<=0); %remove weights which are less than or = to zero, and automatically include them in the selected items
        weight(ind_selected)=[];
        value(ind_selected)=[];

        
        edit_weights{i}=weight; 
        edit_values{i}=value;

    %     figure
    %     plot(values{i})
        [best(i), amount{i}] = LWcase_knapsack(edit_weights{i},edit_values{i}, max_upfront_cost);
        items{i} = find(amount{i});
        items{i}=[ind_selected items{i}+length(ind_selected) ]; %re-add the items which were removed for being below 0, which were sorted, so add the length of the selected values

        TCO_reduction(i)=sum(raw_values{i}(items{i}));
        upfront_cost(i)=sum(wt_unrounded{i}(items{i}));

        amount_combined{i}=[ones(1,length(ind_selected)) amount{i}]; %add the cases < 0 to the amount selected vector

    end


    
    %check to make sure the algorithm didn't take advantage of the
    %efficiency measures which were found without selecting the efficiency
    %measures for actual application

    tco_edit=zeros(1,length(TCO_reduction));
    for g=1:length(TCO_reduction)
        selected{g}=measures{g}(items{g}); %see which tech was selected
        ind_sel_eff=find(ismember(selected{g},eff_measures{g})); %check if the efficiency measure actually was selected
       if ~isempty(ind_sel_eff) %i.e. the tech was selected
            tco_edit(g)=TCO_reduction(g); %include the TCO reduction for that technology. Otherwise, it remains zero
       end
        
    end
    %find the maximum reduction run in the problem after considering
    %infeasible cases where the efficiency measure was kicked out
    [maxval ,maxind]=max(tco_edit);
    
    
    %compile the lw measures string
    lw_str=[];
    for t=1:length(items{maxind})
        lw_str=[lw_str ', ' measures{maxind}{items{maxind}}  ];
    end
    fprintf(['\tTotal lifetime cost of ' num2str(baseline_TCO-maxval) ' USD at this upfront cost: ' num2str(upfront_cost(maxind)) '  found for with this vehicle and this set of measures: '  lw_str   '\n'])
    toc
    
    lt_cost(s)=baseline_TCO-maxval;
    real_up_front(s)=upfront_cost(maxind);
    eff_sel{s}=eff_measures{maxind};
    eff_sel_short{s}=eff_measures_short{maxind};
    num_lw_measures(s)=length(items{maxind});
    lw_tech{s}=lw_str;
    selected_amounts{s}=amount_combined;
%     keyboard
    if plot_scen
    %values for different scenario plots
    figure
    hold
    cc=cool(length(edit_values));
    for m=1:length(edit_values)
        plot(1:length(edit_values{m}),baseline_TCO-edit_values{m},'color',cc(m,:));
    end
    legend(eff_measures);
    ylabel('Total cost of ownership USD');
    % title('Reduction values');
    set(gca,'XTick',1:length(edit_values{m}))
    set( gca(), 'XTickLabel', measures{maxind} )
    set(gca, 'XTickLabelRotation', 45);
%     rotateXLabels( gca(), 45 )
    % keyboard

    %values for different scenario plots
    figure
    hold
    for m=1:length(wt_unrounded)
        plot(1:length(wt_unrounded{m}),wt_unrounded{m},'color',cc(m,:));
    end
    legend(eff_measures);
    % title('Up-front cost values values');
    ylabel('Up-front scenario cost (USD)');
    set(gca,'XTick',1:length(wt_unrounded{m}))
    set( gca(), 'XTickLabel', measures{maxind} )
    set(gca, 'XTickLabelRotation', 45);
%     rotateXLabels( gca(), 45 )
    end
    
    
    figure
    cc=cool(length(edit_values));
    subplot(3,1,1)
    plot(baseline_TCO-best);
    ylabel('Best total cost of ownership USD');
    axis([1 length(items) 0 max(baseline_TCO-best)+1000])
    set(gca,'XTick',1:length(items))
    set( gca(), 'XTickLabel', eff_measures )
%     rotateXLabels( gca(), 45 )
set(gca, 'XTickLabelRotation', 45);
    
    %best case

    subplot(3,1,2)
    hold

    stem(amount_combined{maxind},':ok')

    axis([0 length(wt_unrounded{maxind})+1 -.1 1.1]) %set the axis to 
    legend(eff_measures{maxind});
    ylabel('Technology Selected?');
    title('Best case')
    set(gca,'XTick',1:length(wt_unrounded{maxind}))
    set( gca(), 'XTickLabel', [])
    set(gca,'YTick',[0 1])


    %all cases
    subplot(3,1,3)
    hold
    for n=1:length(items)
        stem(amount_combined{n}+n/100,':ok','color',cc(n,:))
    end
    axis([0 length(wt_unrounded{maxind})+1 -.1 1.1]) %set the axis to 
    legend(eff_measures);
    ylabel('Technology Selected?');
    title('All cases')
    set(gca,'XTick',1:length(wt_unrounded{maxind}))
    set(gca,'YTick',[0 1])
    set( gca(), 'XTickLabel', measures{maxind} )
%     rotateXLabels( gca(), 45 )
    set(gca, 'XTickLabelRotation', 45);
    
%     keyboard
end

%plot the values for the maximum scenario

figure
hold
plot(0:length(max_upfront_vector)-1,max_upfront_vector,'ok','linewidth',2)
plot(0:length(real_up_front)-1,real_up_front,'k','linewidth',2)
text([0:length(max_upfront_vector)-1]+0.2,max_upfront_vector+500,eff_sel_short, 'FontSize',8)
set(gca,'xtick',[])
legend({'max up-front' 'actual up-front'});
xlabel('Scenaro');
ylabel('Up-front cost USD');

figure
% axis([1 length(max_upfront_vector)  0 max(max_upfront_vector)+100])
[ax,h1,h2]=plotyy(0:length(lt_cost)-1,lt_cost,0:length(lt_cost)-1,num_lw_measures,'plot','bar');
set(ax(2),'Ylim',[0,100],'YTick',[0 10 20])
ylabel(ax(2),'# of LW measures');
% set(ax(1),'Ylim',[-2500,5e4])
% axis([1 length(real_up_front)  0 max(max_upfront_vector)+100])
set(h2,'FaceColor','none')
set(h1,'color','k','Marker','s','linewidth',2)

set(ax(1),'xtick',[])
ylabel(ax(1),'Total cost of ownership USD');
set(ax(2),'xtick',[])
legend({ 'optimum TCO' '# of LW measures'});
xlabel('Scenaro');

catch
    a=lasterror
    keyboard
end