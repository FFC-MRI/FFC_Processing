classdef FFCFutureScope < handle
%FFCFUTURESCOPE Own outstanding futures independently of a caller's workspace.
% guard() constructs its callback here, not in the exiting nested workspace.
    properties (Access=private)
        Futures = {}
    end
    methods
        function cleanup = guard(obj)
            cleanup = onCleanup(@() obj.cancelAll());
        end
        function add(obj,future)
            obj.Futures{end+1} = future;
        end
        function release(obj,future)
            match = cellfun(@(value) isequal(value,future),obj.Futures);
            obj.Futures(match) = [];
        end
        function n = count(obj)
            n = numel(obj.Futures);
        end
        function cancelAll(obj)
            futures = obj.Futures;
            obj.Futures = {}; % idempotent, including error unwinding
            for k = 1:numel(futures)
                try, cancel(futures{k}); catch, end
                try, delete(futures{k}); catch, end
            end
        end
        function delete(obj)
            obj.cancelAll();
        end
    end
end
