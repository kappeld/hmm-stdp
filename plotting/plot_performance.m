function loglik_all = plot_performance( src_path )
% plot performance of results in src_path

    [base_path,vu] = snn_process_options( default_options(), 'base_path', '' );
    source_dirs = dir( [base_path,src_path] );    
    source_dirs = source_dirs([source_dirs.isdir]);
    
    loglik_all = [];
    
    for i = 3:length(source_dirs)
        cur_file = [base_path,src_path,source_dirs(i).name,'/performance_kl.mat'];
        if ~exist(cur_file,'file')
            fprintf('skipping file %s!\n', cur_file);
            continue;
        end
        fprintf('reading file %s...\n', cur_file);
        data = load(cur_file,'loglik');
        loglik_all = [ loglik_all, mean(data.loglik,2) ];
    end
end
