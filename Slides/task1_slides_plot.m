function task1_slides_plot(IRF, collector, response, shock, VAR_config, FLAGcumsum, pic_config)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % task1_slides_plot
    %
    % Inputs
    % IRF
    % collector - draw info from numerical IRF
    % reponse - variable number of response
    % shock - variable number of shock
    % VAR_config - object generated by SVAR_config.m
    % FLAGcumsum - logical flag for cumulative sum
    % pic_config - object generated by pic_config.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    debug = 0;
    if debug == 1
        VAR_config = t1_VAR_config;
        pic_config = picture_config;
        IRF = t1_IRF;
        collector = t1_collector;
        FLAGcumsum = false;
        response = 1;
        shock = 1;
    end
 
    % Unspool
    drawMatrix = VAR_config.drawMatrix;
    nhorizon = VAR_config.nhorizon;
    names = VAR_config.names;
 
    fig_fontsize = pic_config.fontsize;
    fig_width = pic_config.width;
    fig_height = pic_config.height;
    pic_dir = pic_config.directory;
 
    colormat = {'red', 'blue', 'green'};
    
    for ii = 1:size(drawMatrix, 1)
     
        % No cumsum
        if FLAGcumsum == 0
            shotgun_plotin = squeeze(IRF(response, shock, 1:nhorizon, 1:collector(ii, 1)));
            % Cumsum
        elseif FLAGcumsum == 1
            shotgun_plotin = cumsum(squeeze(IRF(response, shock, 1:nhorizon, 1:collector(ii, 1))), 1);
        end
     
        figure('Name', ['Task1_' num2str(ii)])
     
        plot(0:nhorizon - 1, shotgun_plotin, 'LineWidth', 1, 'Color', 'blue')
     
        % For bands
        maxsg = max(shotgun_plotin, [], 2);
        minsg = min(shotgun_plotin, [], 2);

        hold on;
        plot(0:nhorizon - 1, [maxsg, minsg], 'LineWidth', 2, 'Color', 'black')
   
        ylabel(['Response of $' names{response} '$'], 'fontsize', 12, 'Interpreter', 'latex');
        xlabel('Horizon', 'fontsize', 12, 'Interpreter', 'latex');
        title(['Total = ' num2str(collector(ii, 3)) '; Inside = ' num2str(collector(ii, 1))], 'Interpreter', 'latex');
        latex_fig(fig_fontsize, fig_width, fig_height);
        tightfig();
        print(gcf, '-depsc2', fullfile(pic_dir, ['slides_t1_' num2str(ii) '.eps']))
     
    end
end