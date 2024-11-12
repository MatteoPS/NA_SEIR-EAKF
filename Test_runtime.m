%Comment the last two line of MODEL_RUN()
%to avoid counting the time to save variable and plotting

startTime = datetime('now');
MODEL_RUN("speedtest"); % Run your model
endTime = datetime('now');
runtime = seconds(endTime - startTime);
disp(['Runtime: ', num2str(runtime), ' seconds']);
save('run_info.mat', 'runtime'); %save the run time (in seconds) in run_info.mat

