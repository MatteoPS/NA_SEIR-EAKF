function Plotting_paper(all_file_name)

load(all_file_name);

cumu_dailyIr_post_perc = cumu_dailyIr_post_rec./population;
cumu_dailyIu_post_perc = cumu_dailyIu_post_rec./population;

startday ='20-Jan-2020';
endday = '31-Mar-2021';

estday1='06-Jun-2020';
estday2='07-Sep-2020';
estday3='15-Mar-2021';


estdays = {estday1, estday2, estday3};
%estdays = {estday1}
loc_of_interest_days = [ "British Columbia" "New York" "Distrito Federal"; ...
                         "Ontario" "Florida" "Estado de Mexico"; ...
                         "Quebec" "California" "Jalisco"]
                         
                         %;...
                         %"Alabama" "North Dakota" "Ohio"];


for i=1:length(estdays)

    loc_of_interest = loc_of_interest_days(i,:);

    Rt=para_post(betamap,:,:).*para_post(2,:,:).*(para_post(alphamaps,:,:)+((1-para_post(alphamaps,:,:)).*para_post(3,:,:)));
    Rt_mean=mean(Rt,2);
    

    alpha_value = 0.3;

    black = [0 0 0];
    blue = [114 147 203]./255;
    red = [211 94 96]./255;
    gray = [128 133 133]./255; gray=gray-0.1;
    green = [132 186 91]./255;
    brown = [171 104 87]./255;
    purple = [144 103 167]./255;
    yellow = [0.9290 0.6940 0.1250];


    prefix_files= "Obs-para";



    figpaper=figure('visible','off');

    t = tiledlayout(length(loc_of_interest),5);
    t.TileSpacing = 'compact';

    for iname = 1:length(loc_of_interest)
        loc_plot_name = loc_of_interest{iname};
        fprintf("plotting %s %d/%d\n", loc_plot_name,iname,length(loc_of_interest))
        loc_plot = table2array(statecodes(strcmp(statecodes.Var2,loc_plot_name),1));

        % 1 - obs
        nexttile;
        hold on
        if iname == 1
            title('\rmModel fitting','FontSize',20)
        end
     
        lb_obs = prctile(squeeze(dailyIr_post_rec(loc_plot,:,:)), 2.5,1);
        ub_obs = prctile(squeeze(dailyIr_post_rec(loc_plot,:,:)), 97.5,1);
        air=area(datetime(startday):datetime(endday), ub_obs, 'EdgeColor',"w", 'FaceColor',purple+0.1, 'FaceAlpha',alpha_value);
        area(datetime(startday):datetime(endday), lb_obs, 'EdgeColor','w', 'FaceColor', 'w');
        %plot(datetime(startday):datetime(endday), dailyincidence(loc_plot,:)', '*','MarkerSize', 1, color=gray);
        lir=plot(datetime(startday):datetime(endday), squeeze(mean(dailyIr_post_rec(loc_plot,:,:),2)),'LineWidth',1.5, color=[purple alpha_value+0.2]);
        lobs=plot(datetime(startday):datetime(endday), obs_case(loc_plot,:),'LineWidth',1, color=[black alpha_value+0.2]); % 7-days smoothed dailyincidence
        
        %xline(reference_date,"--",LineWidth=2)
        xline(datetime(estday1),"--",LineWidth=2)
        xline(datetime(estday2),"--",LineWidth=2)
        xline(datetime(estday3),"--",LineWidth=2)
        hold off
        set(gca,'xtick',[])
        set(gca,'fontsize', 15)
        xlim([datetime(startday) datetime(endday)]);
        box off
        if iname == length(loc_of_interest)
            set(gca, 'xtick', datetime(startday):30:datetime(endday), ...
                'xticklabel', string(datestr(datetime(startday):30:datetime(endday), "mmm-yy")), ...
                'fontsize', 11);
            xtickangle(45); % Rotate x-tick labels
            lgd=legend([air,lobs],"Ir","7d smooth daily incidence", "Orientation","horizontal","Location","southoutside");
            legend('boxoff')
            fontsize(lgd,15,'points')
        end
        ylabel({['\rm\bf' regexprep(loc_plot_name,' ','\n')];['\rm\itPop: ' num2str(population(loc_plot)/1000000,3) ' M']},"Rotation",0,'fontsize',20)



        % 2 - state variables
        nexttile;
        if iname == 1
            title('\rmState variables','FontSize',20,'FontWeight','bold')
        end
        hold on

        lb_obs = prctile(squeeze(S_post(loc_plot,:,:)), 2.5,2);
        ub_obs = prctile(squeeze(S_post(loc_plot,:,:)), 97.5,2);
        aS=area(datetime(startday):datetime(endday), ub_obs, 'EdgeColor',"w", 'FaceColor',green, 'FaceAlpha',alpha_value);
        area(datetime(startday):datetime(endday), lb_obs, 'EdgeColor','w', 'FaceColor', 'w');
        lS=plot(datetime(startday):datetime(endday), mean(S_post(loc_plot,:,:),3),"-", 'LineWidth',1.5, 'Color', green -0.1);

        hold on
        lb_obs = prctile(squeeze(cumu_dailyIu_post_perc(loc_plot,:,:)), 2.5,2);
        ub_obs = prctile(squeeze(cumu_dailyIu_post_perc(loc_plot,:,:)), 97.5,2);
        aIu=area(datetime(startday):datetime(endday), ub_obs, 'EdgeColor',"w", 'FaceColor',brown+0.1, 'FaceAlpha',alpha_value);
        area(datetime(startday):datetime(endday), lb_obs, 'EdgeColor','w', 'FaceColor', 'w');
        lIu=plot(datetime(startday):datetime(endday), mean(cumu_dailyIu_post_perc(loc_plot,:,:),3),'-', 'LineWidth',1.5, 'Color', brown - 0.1);


        hold on
        lb_obs = prctile(squeeze(cumu_dailyIr_post_perc(loc_plot,:,:)), 2.5,2);
        ub_obs = prctile(squeeze(cumu_dailyIr_post_perc(loc_plot,:,:)), 97.5,2);
        aIr=area(datetime(startday):datetime(endday), ub_obs, 'EdgeColor',"w", 'FaceColor',purple+0.1, 'FaceAlpha',alpha_value);
        area(datetime(startday):datetime(endday), lb_obs, 'EdgeColor','w', 'FaceColor', 'w');
        lIr=plot(datetime(startday):datetime(endday), mean(cumu_dailyIr_post_perc(loc_plot,:,:),3),'-', 'LineWidth',1, 'Color', purple -0.1);


        %xline(reference_date,"--",LineWidth=2)
        xline(datetime(estday1),"--",LineWidth=2)
        xline(datetime(estday2),"--",LineWidth=2)
        xline(datetime(estday3),"--",LineWidth=2)


        xlim([datetime(startday) datetime(endday)]);
        ylim([0 inf])


        %fontsize(lgd,10,'points')
        %title("% State Variables");
        set(gca,'xtick',[])
        set(gca,'fontsize', 15)
        box off
        hold off
        if iname == length(loc_of_interest)
            set(gca, 'xtick', datetime(startday):30:datetime(endday), ...
                'xticklabel', string(datestr(datetime(startday):30:datetime(endday), "mmm-yy")), ...
                'fontsize', 11);
            xtickangle(45); % Rotate x-tick labels
            lgd=legend([aS,aIu,aIr],"%S","cumu %Iu","cumu %Ir", "Orientation","horizontal","Location","southoutside");
            legend('boxoff')
            fontsize(lgd,15,'points')
        end

        %  3 - alpha
        nexttile;
        if iname == 1
            title('\alpha','FontSize',20)
        end
        lb_alpha = prctile(squeeze(para_post(alphamaps(loc_plot), :,:)), 2.5,1);
        ub_alpha = prctile(squeeze(para_post(alphamaps(loc_plot), :,:)), 97.5,1);
        hold on
        aest=area(datetime(startday):datetime(endday), ub_alpha, 'EdgeColor',"w", 'FaceColor',red, 'FaceAlpha',alpha_value);
        area(datetime(startday):datetime(endday), lb_alpha, 'EdgeColor','w', 'FaceColor', 'w');
        plot(datetime(startday):datetime(endday), squeeze(mean(para_post(alphamaps(loc_plot,:),:,:),2)),'LineWidth',1.5, color=red-0.2);
        %plot(paramin_rec(alphamaps(loc_plot),:),'k--')
        hold off

        %xline(reference_date,"--",LineWidth=2)
        xline(datetime(estday1),"--",LineWidth=2)
        xline(datetime(estday2),"--",LineWidth=2)
        xline(datetime(estday3),"--",LineWidth=2)

        %title('\alpha')
        xlim([datetime(startday) datetime(endday)]);
        ylim([alphalow alphaup])
        set(gca,'xtick',[])
        set(gca,'fontsize', 15)
        box off
        hold off
        if iname == length(loc_of_interest)
            set(gca, 'xtick', datetime(startday):30:datetime(endday), ...
                'xticklabel', string(datestr(datetime(startday):30:datetime(endday), "mmm-yy")), ...
                'fontsize', 11);
            xtickangle(45); % Rotate x-tick labels
            lgd=legend(aest,"\alpha","Location","southoutside");
            legend('boxoff')
            fontsize(lgd,15,'points')
        end


        %4 - beta
        nexttile;
        if iname == 1
            title('\beta','FontSize',20)
        end
        hold on
        lb_beta = prctile(squeeze(para_post(betamap(loc_plot), :,:)), 2.5,1);
        ub_beta = prctile(squeeze(para_post(betamap(loc_plot), :,:)), 97.5,1);
        aest=area(datetime(startday):datetime(endday), ub_beta, 'EdgeColor',"w", 'FaceColor',blue, 'FaceAlpha',alpha_value);
        area(datetime(startday):datetime(endday), lb_beta, 'EdgeColor','w', 'FaceColor', 'w');
        plot(datetime(startday):datetime(endday),squeeze(mean(para_post(betamap(loc_plot,:),:,:),2)),'LineWidth',1.5, color=blue-0.2);
        set(gca,'fontsize', 15)
        set(gca,'xtick',[])
        xlim([datetime(startday) datetime(endday)]);
        ylim([0 inf])

        %xline(reference_date,"--",LineWidth=2)
        xline(datetime(estday1),"--",LineWidth=2)
        xline(datetime(estday2),"--",LineWidth=2)
        xline(datetime(estday3),"--",LineWidth=2)
        
        box off

        hold off
        ylim([betalow,betaup])
        if iname == length(loc_of_interest)
            set(gca, 'xtick', datetime(startday):30:datetime(endday), ...
                'xticklabel', string(datestr(datetime(startday):30:datetime(endday), "mmm-yy")), ...
                'fontsize', 11);
            xtickangle(45); % Rotate x-tick labels
            lgd=legend(aest,"\beta","Location","southoutside");
            legend('boxoff')
            fontsize(lgd,15,'points')
        end

        % % 5 - Rt
        nexttile;
        if iname == 1
            title('\rm\itRt','FontSize',20)
        end
        hold on
        lb_rt = prctile(squeeze(Rt(loc_plot,:,:)), 2.5,1);
        ub_rt = prctile(squeeze(Rt(loc_plot,:,:)), 97.5,1);
        aest=area(datetime(startday):datetime(endday), ub_rt, 'EdgeColor',"w", 'FaceColor',yellow-0.1, 'FaceAlpha',alpha_value);
        area(datetime(startday):datetime(endday), lb_rt, 'EdgeColor','w', 'FaceColor', 'w');
        yline(1, color=gray)
        plot(datetime(startday):datetime(endday),Rt_mean(loc_plot,:)','LineWidth',1.5, color=yellow-0.1);

        %xline(reference_date,"--",LineWidth=2)
        xline(datetime(estday1),"--",LineWidth=2)
        xline(datetime(estday2),"--",LineWidth=2)
        xline(datetime(estday3),"--",LineWidth=2)
        %box off
        hold off

        set(gca,'fontsize', 15)
        set(gca,'xtick',[])
        xlim([datetime(startday) datetime(endday)]);

        ylim([0 inf])

        if iname == length(loc_of_interest)
            set(gca, 'xtick', datetime(startday):30:datetime(endday), ...
                'xticklabel', string(datestr(datetime(startday):30:datetime(endday), "mmm-yy")), ...
                'fontsize', 11);
            xtickangle(45); % Rotate x-tick labels
            lgd=legend(aest,"\rm\itRt","Location","southoutside");
            legend('boxoff')
            fontsize(lgd,15,'points')
        end

    end

    file_figpaper = strjoin(['Output/' i prefix_files '.svg'],'');
    set(figpaper,'PaperUnits','inches','PaperPosition',[0 0 24 12])
    print(figpaper,file_figpaper,'-dsvg','-r450');
end






%%%legends

