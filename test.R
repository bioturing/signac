download.file("https://www.dropbox.com/s/jukonxw5g04n8ra/sim2.row.h5?dl=1", "sim2.row.h5")
cluster <- rep(2, 2499)
cluster[1:500] <- 1
system.time(Signac:::HarmonyMarkerH5("sim2.row.h5", cluster))
