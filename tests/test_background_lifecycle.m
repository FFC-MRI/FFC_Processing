function tests = test_background_lifecycle
% Exercise queued result cancellation and replacement on MATLAB's event loop.
    tests = functiontests(localfunctions);
end
function setupOnce(testCase)
    testCase.TestData.oldPath = path;
    addpath(fileparts(fileparts(mfilename('fullpath')))); ffc_setup();
end
function teardownOnce(testCase)
    path(testCase.TestData.oldPath);
end
function testCancelledJobNeverPublishes(testCase)
    runner = FFCJobRunner(); cleanup = onCleanup(@() delete(runner));
    published = false; failure = [];
    runner.start('preview',ones(4,4,'single'),struct('filter','None','params',struct()), ...
        @completed,@failed);
    runner.cancel();
    pause(0.15); drawnow;
    verifyFalse(testCase,runner.isRunning());
    verifyFalse(testCase,published);
    verifyEmpty(testCase,failure);
    function completed(~), published = true; end
    function failed(err), failure = err; end
end
function testReplacementPublishesOnlyLatestRequest(testCase)
    runner = FFCJobRunner(); cleanup = onCleanup(@() delete(runner));
    outputs = {}; failure = []; payload = struct('filter','None','params',struct());
    runner.start('preview',ones(4,4,'single'),payload,@completed,@failed);
    runner.start('preview',2*ones(4,4,'single'),payload,@completed,@failed);
    started = tic;
    while isempty(outputs) && isempty(failure) && toc(started) < 30
        pause(0.02); drawnow;
    end
    verifyEmpty(testCase,failure);
    assertEqual(testCase,numel(outputs),1);
    verifyEqual(testCase,outputs{1},2*ones(4,4,'single'));
    verifyFalse(testCase,runner.isRunning());
    function completed(value), outputs{end+1} = value; end
    function failed(err), failure = err; end
end
