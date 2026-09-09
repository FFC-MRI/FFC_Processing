classdef FFCJobRunner < handle
%FFCJOBRUNNER One cancellable task, with results delivered on the UI thread.
% New requests cancel previous work. Worker code never accesses graphics.
    properties (Access=private)
        Future = []
        PollTimer = []
        Complete = []
        Failure = []
        Work = {}
    end
    methods
        function start(obj, kind, input, payload, complete, failure)
            obj.cancel();
            obj.Complete = complete; obj.Failure = failure;
            obj.Work = {kind,input,payload};
            try
                obj.Future = parfeval(backgroundPool, @ffc_compute_task, 1, ...
                    kind, input, payload);
                obj.PollTimer = timer('ExecutionMode','fixedSpacing', ...
                    'Period',0.1,'BusyMode','drop','TimerFcn',@(~,~) obj.poll());
                start(obj.PollTimer);
            catch ME
                obj.cancel();
                failure(ME);
            end
        end
        function cancel(obj)
            if ~isempty(obj.PollTimer)
                stop(obj.PollTimer); delete(obj.PollTimer); obj.PollTimer = [];
            end
            if ~isempty(obj.Future)
                cancel(obj.Future); delete(obj.Future); obj.Future = [];
            end
            obj.Complete = []; obj.Failure = [];
            obj.Work = {};
        end
        function value = isRunning(obj)
            value = ~isempty(obj.Future);
        end
        function delete(obj)
            obj.cancel();
        end
    end
    methods (Access=private)
        function poll(obj)
            if isempty(obj.Future) || ~strcmp(obj.Future.State,'finished'), return; end
            complete = obj.Complete; failure = obj.Failure;
            try
                result = fetchOutputs(obj.Future);
            catch ME
                work = obj.Work;
                obj.cancel();
                report = lower(getReport(ME,'extended','hyperlinks','off'));
                unsupported = contains(report,'not supported') || contains(report,'unsupported');
                if unsupported && (contains(report,'thread-based') || ...
                        contains(report,'thread based') || contains(report,'threadpool'))
                    % Some toolbox routines (notably legacy Wavelet calls)
                    % cannot run in threads. Preserve functionality on client.
                    warning('FFCJobRunner:ForegroundFallback', ...
                        'This task uses a function unsupported in background threads; running on the client.');
                    try
                        result = ffc_compute_task(work{:});
                    catch fallbackError
                        failure(fallbackError);
                        return;
                    end
                    complete(result);
                else
                    failure(ME);
                end
                return;
            end
            obj.cancel();
            complete(result);
        end
    end
end
