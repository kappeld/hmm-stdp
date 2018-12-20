function plot_mean_var( data, it, color )
% plot_mean_var( data, it, color )
%
% plot mean and variance (as shaded area around the mean)
%
% 01.12.2012
%

    hold on;
    h = area( [ it', it' ], [ mean(data,2)-std(data,1,2), 2*std(data,1,2) ], ...
              'LineStyle', 'none', 'FaceColor', min(1,color+0.7) );
    delete(h(1));
    plot( it, mean(data,2), 'Color', color );
end
