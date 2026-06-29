function queue = createParallelProgressBar(totalIterations)
    % createParallelProgressBar Initializes a progress bar for parallel
    % computations with dynamic color changing from dark orange to blue.
    %
    % Args:
    %     totalIterations (int): Total number of iterations for the
    %                            progress bar.
    %
    % Returns:
    %     queue (parallel.pool.DataQueue): DataQueue to receive progress
    %                                      updates.
    %
    % Example usage:
    %     numSamples = 100;
    %     queue = createParallelProgressBar(numSamples);
    %     parfor i = 1:numSamples
    %         pause(0.1); % simulate work
    %         send(queue, i);
    %     end

    % Create waitbar
    progressBar = waitbar(0, 'Processing...', ...
        'Name', 'Computation Progress');

    % Define colors for interpolation
    colorStart = [171, 94, 0] / 255; % orange-brown
    colorEnd   = [12, 123, 220] / 255; % blue

    % Initialize DataQueue
    queue = parallel.pool.DataQueue;

    % Persistent counter
    persistent count
    count = 0;

    % Nested update function
    function updateProgress(~)
        count = count + 1;
        shareComplete = count / totalIterations;

        % Update waitbar percentage
        waitbar(shareComplete, progressBar, ...
            sprintf('Progress: %d%%', round(shareComplete*100)));

        % Update bar color (access waitbar's patch object)
        h = findobj(progressBar, 'Type', 'Patch');
        if ~isempty(h) && isvalid(h)
            currentColor = (1 - shareComplete) * colorStart + ...
                            shareComplete * colorEnd;
            set(h, 'FaceColor', currentColor);
        end

        % Close when complete
        if count == totalIterations
            pause(0.2); % let user see it finish
            close(progressBar);
            count = [];
        end
    end

    % Attach listener
    afterEach(queue, @updateProgress);
end