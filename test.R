mtx.dir <- "/Users/bioturing/Code/HuyMarkerBenchmark/data.2/pbmc4k_batch2/sim_batch2_logfc_2_checkpoints_1.7_2.2"
obj <- CreateSignacObject(mtx.dir, "mtx")
system.time(Signac::Harmony(mtx = obj@raw.data,
                ordering = 6,
                col_idx = seq(1, 500)))