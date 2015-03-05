%%% Knapsack 0-1 Lightweighting

%E.Wilhelm SUTD 2014

try
    
    
% add NPV calculations, and sensitivity analysis

clear
clc
close all

tic

m_baseline=1711; %kg baseline vehicle weight from EPA Lotus Engineering report page 141
PP=23935; %USD MSRP 2014 Ford Fusion
FP=0.95; %USD/L
km=21562; %km
life=11; %years

% life_vec=ceil(normrnd(11,3,1000,1)); %slifetime sensitivity

fuel_vec=normrnd(0.95,.2,1000,1); %fuel price sensitivity

%steps for the life vector
steps=length(fuel_vec);


BE=9; %L/100km baseline fusion fuel consumption

dr=0.05; %for sensitivity analysis
% dr=[0.5 .25 .05 0]; %discount rate vector
% scen_leg_A={'max up-front' 'actual up-front'}; %single case
scen_leg_A={'max up-front' 'actual @ r=0.5' 'actual @ r=0.25' 'actual @ r=0.05' 'actual @ r=0'};

% scen_leg_B={ 'optimal TCO' '# LW measures'}; %single case
scen_leg_B={ 'opt. TCO @ r=0.5' 'opt. TCO @ r=0.25' 'opt. TCO @ r=0.05' 'opt. TCO @ r=0' };

% max_upfront_vector=[10000 5000 2000 1000 750 500 400 300 200 100 50 1]; %maximum up-front cost
max_upfront_vector=[1000 10000]; %for sensitivitiy analysis

plot_scen=0; %flag to plot scenarios on a case-by-case basis

verbose=0; %flag to output results from one run of full upfront cost vectors

dr_scen=0; %flag to plot discount rate scenario output

life_sens=0; %flag to plot monte-carlo life sensitivity results

fuel_sens=1;

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



%for 0% discount rate
% baseline_TCO=PP+(BE-dE)*FP*km*life/100+0+0; %the standard vehicle's total lifetime costs, without additional technology costs
h = waitbar(0,'Please wait...');


%run through all of the other scenarios


% for z=1:length(dr) %discount rate sensitivity case
for z=1:length(fuel_vec) %life sensitivity case
    
%for real NPV calcs
Bl_yearly=(BE-dE)*fuel_vec(z)*km/100+0+0;
% baseline_TCO=pvvar([PP ones(1,life)*Bl_yearly],dr(z)); %for discount rate sensitivity case
baseline_TCO=pvvar([PP ones(1,life)*Bl_yearly],dr(end)); %for lifetime sensitivity case

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
                %0% discount rate
    %             reduced_TCO(j)=PP+(BE-dE_choice)*FP*km*life/100+wt_raw(j); 
                yearly=(BE-dE_choice)*fuel_vec(z)*km/100;
%                 reduced_TCO(j)=pvvar([PP+wt_raw(j) ones(1,life)*yearly],dr(z));% for the discount rate sensitivity case
                reduced_TCO(j)=pvvar([PP+wt_raw(j) ones(1,life)*yearly],dr(end)); %for the life sensitivity case
                value(j)=baseline_TCO-reduced_TCO(j);     
            end
            %run the efficiency choice
            weight(j+1)=ceil(eff(i,2));
            wt_raw(j+1)=eff(i,2);
            dE_choice=eff(i,1); % add the direct efficiency reduction from the powertrain to the sensitivity-adusted efficiency reduction from the weight reduction
            %0% discount rate
    %         reduced_TCO(j+1)=PP+(BE-dE_choice)*FP*km*life/100+wt_raw(j+1); 
            %real discount rate
            yearly=(BE-dE_choice)*fuel_vec(z)*km/100;
%             reduced_TCO(j+1)=pvvar([PP+wt_raw(j+1) ones(1,life)*yearly],dr(z));% discount rate sensitivity case
            reduced_TCO(j+1)=pvvar([PP+wt_raw(j+1) ones(1,life)*yearly],dr(end)); %for the life sensitivity case
            
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

        if verbose
            fprintf(['\tTotal lifetime cost of ' num2str(baseline_TCO-maxval) ' USD at this upfront cost: ' num2str(upfront_cost(maxind)) '  found for with this vehicle and this set of measures: '  lw_str   '\n'])
        end
        toc

        lt_cost(s)=baseline_TCO-maxval; %take the baseline TCO, and subtract the savings above baseline to yield the actual TCO
        
        
        real_up_front(s)=upfront_cost(maxind); 
        eff_sel{s}=eff_measures{maxind};
        eff_sel_short{s}=eff_measures_short{maxind};
        num_lw_measures(s)=length(items{maxind});
        lw_tech{s}=lw_str;
        selected_amounts{s}=amount_combined;
%         keyboard
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
        rotateXLabels( gca(), 45 )
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
        rotateXLabels( gca(), 45 )
        end


        if verbose %plot the scenario results
        figure
        cc=cool(length(edit_values));
        subplot(3,1,1)
        plot(baseline_TCO-best);
        ylabel('Best total cost of ownership USD');
        axis([1 length(items) 0 max(baseline_TCO-best)+1000])
        set(gca,'XTick',1:length(items))
        set( gca(), 'XTickLabel', eff_measures )
        rotateXLabels( gca(), 45 )


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
        rotateXLabels( gca(), 45 )
        end

    %     keyboard
    end
    numvec{z}=num_lw_measures;
    ltvec{z}=lt_cost;
    realvec{z}=real_up_front;
    effsels{z}=eff_sel_short;
waitbar(z / steps)
end
close(h) 

if fuel_sens %plot the lifetime sensitivity results

    
    figure
    effsels_test=[effsels{:}];
    %odd numbers representing the first max up-front
    num_none=length(find(ismember(effsels_test(1:2:end),'None')));
    num_ss=length(find(ismember(effsels_test(1:2:end),'S-S')));
    num_dst=length(find(ismember(effsels_test(1:2:end),'DS/T')));
    num_full=length(find(ismember(effsels_test(1:2:end),'Full')));
    num_plug=length(find(ismember(effsels_test(1:2:end),'Plug-in')));
    num_ev=length(find(ismember(effsels_test(1:2:end),'EV')));
    bar_data(1,:)=[num_none num_dst num_ss num_full num_plug num_ev];
     %even numbers representing the second max up-front
    num_none=length(find(ismember(effsels_test(2:2:end),'None')));
    num_ss=length(find(ismember(effsels_test(2:2:end),'S-S')));
    num_dst=length(find(ismember(effsels_test(2:2:end),'DS/T')));
    num_full=length(find(ismember(effsels_test(2:2:end),'Full')));
    num_plug=length(find(ismember(effsels_test(2:2:end),'Plug-in')));
    num_ev=length(find(ismember(effsels_test(2:2:end),'EV')));
    bar_data(2,:)=[num_none num_dst num_ss num_full num_plug num_ev];
    
    hbar=bar(bar_data');
    set(gca, 'XTickLabel',eff_measures_short)
%     set(hbar,'FaceColor','k');
    legend({'1000' '10,000'});
    ylabel('Number of vehicles')

    figure
    lt=cell2mat(ltvec);
    plot(fuel_vec,lt(1:2:end),'.k','MarkerSize',11);
    hold
    plot(fuel_vec,lt(2:2:end),'.r','MarkerSize',11);
    
    xlabel('Fuel Price ($/L)');
    ylabel('Optimal total cost of ownership')
    legend({'1000' '10,000'});
    
    figure
    rv=cell2mat(realvec);
    plot(fuel_vec,rv(1:2:end),'.k','MarkerSize',11);
    hold
    plot(fuel_vec,rv(2:2:end),'.r','MarkerSize',11);
    xlabel('Fuel Price ($/L)');
    
    ylabel('Real marginal cost')
    legend({'1000' '10,000'});
    
    
    figure
    hold
    plot(rv(1:2:end),lt(1:2:end),'.k','MarkerSize',8)
    plot(rv(2:2:end),lt(2:2:end),'.r','MarkerSize',8)
    xlabel('Real marginal cost with max. 1000 USD up-front')
    ylabel('Optimal total cost of ownership')
    
    legend({'1000' '10,000'});
    
    figure
    hold
    nv=cell2mat(numvec);
    plot(nv(1:2:end),lt(1:2:end),'.k','MarkerSize',8)
    plot(nv(2:2:end),lt(2:2:end),'.r','MarkerSize',8)

    xlabel('Number of lightweighting measures selected')
    ylabel('Optimal total cost of ownership')
    xlim([0 22])
    
end

if life_sens %plot the lifetime sensitivity results
%     figure
% 
%     plot(cell2mat(realvec),cell2mat(ltvec),'.k','MarkerSize',8)
%     xlabel('Real marginal cost with max. 1000 USD up-front')
%     ylabel('Optimal total cost of ownership')
%     
%     figure
%     plot(cell2mat(numvec),cell2mat(ltvec),'.k','MarkerSize',8)
%     xlabel('Number of lightweighting measures selected')
%     ylabel('Optimal total cost of ownership')
%     xlim([0 22])
    
    figure
    effsels_test=[effsels{:}];
    %odd numbers representing the first max up-front
    num_none=length(find(ismember(effsels_test(1:2:end),'None')));
    num_ss=length(find(ismember(effsels_test(1:2:end),'S-S')));
    num_dst=length(find(ismember(effsels_test(1:2:end),'DS/T')));
    num_full=length(find(ismember(effsels_test(1:2:end),'Full')));
    num_plug=length(find(ismember(effsels_test(1:2:end),'Plug-in')));
    num_ev=length(find(ismember(effsels_test(1:2:end),'EV')));
    bar_data(1,:)=[num_none num_dst num_ss num_full num_plug num_ev];
     %even numbers representing the second max up-front
    num_none=length(find(ismember(effsels_test(2:2:end),'None')));
    num_ss=length(find(ismember(effsels_test(2:2:end),'S-S')));
    num_dst=length(find(ismember(effsels_test(2:2:end),'DS/T')));
    num_full=length(find(ismember(effsels_test(2:2:end),'Full')));
    num_plug=length(find(ismember(effsels_test(2:2:end),'Plug-in')));
    num_ev=length(find(ismember(effsels_test(2:2:end),'EV')));
    bar_data(2,:)=[num_none num_dst num_ss num_full num_plug num_ev];
    
    hbar=bar(bar_data');
    set(gca, 'XTickLabel',eff_measures_short)
%     set(hbar,'FaceColor','k');
    legend({'1000' '10,000'});
    ylabel('Number of vehicles')

    figure
    lt=cell2mat(ltvec);
    plot(life_vec,lt(1:2:end),'.k','MarkerSize',11);
    hold
    plot(life_vec,lt(2:2:end),'.r','MarkerSize',11);
    xlabel('Operating Life (years)');
    
    ylabel('Optimal total cost of ownership')
    legend({'1000' '10,000'});
    
    figure
    rv=cell2mat(realvec);
    plot(life_vec,rv(1:2:end),'.k','MarkerSize',11);
    hold
    plot(life_vec,rv(2:2:end),'.r','MarkerSize',11);
    xlabel('Operating Life (years)');
    ylabel('Real marginal cost')
    legend({'1000' '10,000'});
    
end


%plot the values for the standard scenarios and discount rate

if dr_scen
    cmap = hsv(length(realvec));
    figure
    hold
    plot(0:length(max_upfront_vector)-1,max_upfront_vector,'ok','linewidth',2)
    for f=1:length(realvec)
        real_up_front=cell2mat(realvec(f));
        eff_sel_short=effsels{f};

        plot(0:length(real_up_front)-1,real_up_front,'-s','Color',cmap(f,:),'linewidth',2)
        text([0:length(real_up_front)-1]+0.2,real_up_front+500,eff_sel_short, 'FontSize',8)
    end
    set(gca,'xtick',[])
    legend(scen_leg_A);
    xlabel('Scenaro');
    ylabel('Up-front cost USD');

    figure
    hold
    % axis([1 length(max_upfront_vector)  0 max(max_upfront_vector)+100])
    for g=1:length(ltvec)
        lt_cost=cell2mat(ltvec(g));
        num_lw_measures=cell2mat(numvec(g));
        [ax,h1,h2]=plotyy(0:length(lt_cost)-1,lt_cost,0:length(lt_cost)-1,num_lw_measures,'plot','bar');
        set(ax(2),'Ylim',[0,100],'YTick',[0 10 20])

        ylabel(ax(2),'# of LW measures');
        % set(ax(1),'Ylim',[-2500,5e4])
        % axis([1 length(real_up_front)  0 max(max_upfront_vector)+100])
        set(h2,'FaceColor','none','EdgeColor',cmap(g,:))
        set(h1,'color',cmap(g,:),'Marker','s','linewidth',2)
        set(ax(2),'xtick',[])
        set(ax(1),'xtick',[])
    end
    set(ax(1),'Ylim',[2e4 5e4],'YTick',[2e4 3e4 4e4 5e4])

    ylabel(ax(1),'Total cost of ownership USD');
    set(ax(2),'xtick',[])
    legend(scen_leg_B);
    xlabel('Scenaro');
end

catch
    a=lasterror
    keyboard
end