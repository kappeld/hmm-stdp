function sorted_hist_plot( data, num_bins, idx, field_name )
% Plot spike histrogram for given data field.
%
% David Kappel 11.08.2011
%

    data_field = data.(field_name)(1,:);

    h = hist( data_field, 1:num_bins );
    
    barh( 1:length(idx), h(idx)./(data.time(end)-data.time(1)) );
    
    min_y_dist = 0.02;

    y_dist = 0.5;

    if y_dist/length(idx) < min_y_dist
        y_dist = min_y_dist*length(idx);
    end

    ylim([1-y_dist,length(idx)+y_dist]);
    set(gca,'YDir','reverse')
end
