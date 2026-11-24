function tests = test_future_scope
    tests = functiontests(localfunctions);
end
function setupOnce(tc)
    tc.TestData.oldPath = path;
    addpath(fileparts(fileparts(mfilename('fullpath'))));
end
function teardownOnce(tc)
    path(tc.TestData.oldPath);
end
function testNormalNestedWorkspaceExit(tc)
    lastwarn('');
    scope = normalExit();
    verifyEqual(tc,scope.count(),0);
    [message,~] = lastwarn;
    verifyEmpty(tc,message);
end
function scope = normalExit()
    jobs = {}; scope = FFCFutureScope();
    cleanup = scope.guard(); %#ok<NASGU>
    touch();
    function touch()
        jobs{end+1} = 'shared nested variable';
    end
end
function testErrorUnwindOwnsFutureIndependently(tc)
    assumeTrue(tc,license('test','Distrib_Computing_Toolbox'));
    future = parfeval(backgroundPool,@pause,0,30);
    scope = FFCFutureScope();
    lastwarn('');
    verifyError(tc,@() failWithFuture(scope,future),'fixture:Expected');
    verifyEqual(tc,scope.count(),0);
    verifyFalse(tc,isvalid(future));
    [message,~] = lastwarn; verifyEmpty(tc,message);
end
function failWithFuture(scope,future)
    jobs = {future}; scope.add(future);
    cleanup = scope.guard(); %#ok<NASGU>
    fail();
    function fail()
        assert(numel(jobs)==1);
        error('fixture:Expected','Intentional failure');
    end
end
function testReleaseLeavesUnownedFutureUntouched(tc)
    assumeTrue(tc,license('test','Distrib_Computing_Toolbox'));
    future = parfeval(backgroundPool,@pause,0,30);
    cleanup = onCleanup(@() delete(future)); %#ok<NASGU>
    scope = FFCFutureScope(); scope.add(future); scope.release(future);
    scope.cancelAll(); scope.cancelAll();
    verifyTrue(tc,isvalid(future)); verifyEqual(tc,scope.count(),0);
    cancel(future);
end
function testWriteDisableShortCircuitsSave(tc)
    had = ispref('FFCProcessing','DiskCacheWriteEnabled');
    old = []; if had, old = getpref('FFCProcessing','DiskCacheWriteEnabled'); end
    cleanup = onCleanup(@() restoreWritePref(had,old)); %#ok<NASGU>
    setpref('FFCProcessing','DiskCacheWriteEnabled',false);
    verifyEmpty(tc,ffc_disk_cache('save','nonexistent-folder',[]));
end
function restoreWritePref(had,old)
    if had, setpref('FFCProcessing','DiskCacheWriteEnabled',old);
    else, rmpref('FFCProcessing','DiskCacheWriteEnabled'); end
end
