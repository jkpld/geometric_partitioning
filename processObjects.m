function [cuts, Info] = processObjects(pixels, B, n, K, I, S, r0, isConvex, nRows, options)

% See also COMPUTEOBJECTSEEDPOINTS CREATEOBJECTIMAGES MODELPARTICLEDYNAMICS

    % Number of input objects
    N = numel(pixels);

    % If there is only one object, then turn off parallel computing.
    if N == 1
        options.Use_Parallel = false;
    end


    % If we are computing in parallel, then first convert the options class
    % element to a structure to prevent reinitiallization on transfer to
    % each worker.
    if options.Use_Parallel
        warning('off','MATLAB:structOnObject')
        options = struct(options);
        warning('on','MATLAB:structOnObject')
    end

    % Initialize sliced variables
    Info = cell(N,1);
    cuts = cell(N,1);

    if options.Use_Parallel
        parfor obj = 1:N
            if isConvex(obj)
                cuts{obj} = zeros(0,4);
                continue
            end
            % ===========================================================
            % This could be a good place to put code for detecting
            % seed-points
            % ===========================================================
            
            % Create mask and interior potential images for the object
            [objI, objS] = createObjectImages(pixels{obj}, nRows, I{obj}, S{obj});
            [cuts{obj}, Info{obj}] = geometric_partition(B{obj}, n{obj}, K{obj}, objI, objS, r0{obj}, options);

        end
    else
        % Display the progress of the calculation (we can do this since we are not computing the objects in parallel)
        generateDisplayAt = unique(round(linspace(1,N,7)));
        processTimes = nan(1,N);
        fprintf('Starting partitioning\n')

        for obj = 1:N
            if isConvex(obj)
                cuts{obj} = zeros(0,4);
                continue
            end
            procTime = tic;

            % ===========================================================
            % This could be a good place to put code for detecting
            % seed-points
            % ===========================================================
            
            % Create mask and interior potential images for the object
            [objI, objS] = createObjectImages(pixels{obj}, nRows, I{obj}, S{obj});
            [cuts{obj}, Info{obj}] = geometric_partition(B{obj}, n{obj}, K{obj}, objI, objS, r0{obj}, options);

            processTimes(obj) = toc(procTime);
            if any(obj == generateDisplayAt)
                fprintf('%s >> %d/%d (%0.2f/%0.2f)...\n',datestr(now,31),obj,N,sum(processTimes(1:obj),'omitnan')/60,mean(processTimes(1:obj),'omitnan')*N/60)
            end % if
        end % for
        fprintf('Finished!\n')
    end % if
end % function
