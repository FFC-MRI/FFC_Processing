# TEST29 — future cleanup fix and cache-write control

Close the older app, extract separately, and run ffc_start_cvi from this folder.

The reported warning occurred after processing, while MATLAB destroyed the
nested function workspace. onCleanup(@cancelJobs) referenced the shared jobs
variable after that variable had been destroyed. Both raw loading and batch
DICOM export now register futures with FFCFutureScope. Its cleanup callback is
created in a separate class-method workspace and accesses only the owner's
future list. Completed futures are released before deletion. Cancellation and
error unwinding cancel/delete only this operation's outstanding futures, leaving
the parallel pool available. Cleanup is idempotent.

The warning is a cleanup/lifecycle fault; it does not by itself indicate that
reconstructed images or already completed DICOM files are invalid.

Your reported 100.0 s of cache writes is a substantial part of the 228.4 s
opening time. Disk writes currently save the whole scan object using compressed
MAT v7.3, and can be expensive. Their duration alone cannot distinguish compression
cost, data volume and storage performance. The default cache path is under
prefdir; on managed PCs that location can be redirected. Check the actual path:

```matlab
getpref('FFCProcessing','DiskCacheDirectory', ...
    fullfile(prefdir,'FFC_Reconstruction_Cache_TEST26'))
```

TEST29 adds a read-existing-caches-only option. To avoid the cost of new cache
writes while retaining existing cache hits:

```matlab
setpref('FFCProcessing','DiskCacheWriteEnabled',false)
```

Re-enable writes when desired (the default remains true):

```matlab
setpref('FFCProcessing','DiskCacheWriteEnabled',true)
```

This is independent of DiskCacheEnabled. Reprocessing still bypasses existing
cache reads. With writes disabled, uncached scans will need reconstruction on
each reopening. No cached files are deleted by changing this setting.

Opening summaries now print pool startup separately. Worker timings overlap;
adding all the reported stage totals does not equal wall time. Remaining time
also includes transfers, waiting, graphics updates and other overhead.

All TEST28 features are retained. Numerical reconstruction and DICOM metadata,
geometry and writer code are unchanged. MATLAB runtime validation is unavailable
in the build environment; static checks and regression test sources are included.
Run in MATLAB:

```matlab
ffc_setup
results = runtests('tests');
disp(table(results)); assertSuccess(results)
```

New tests cover normal nested-workspace exit, error unwinding with an active
future, release/idempotence, and disabled cache writes. Recheck a normal load,
cancelled load and all-participant export to confirm the destructor warning has
gone in your MATLAB version.
