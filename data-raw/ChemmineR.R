library("ChemmineR") # Loads the package
data(sdfsample)
sdfset <- sdfsample
sdfset # Returns summary of SDFset
job1 <- launchCMTool("pubchemID2SDF", 2244)
status(job1)
result1 <- result(job1)

job2 <- launchCMTool("Fingerprint Search", result1, "Similarity Cutoff" = 0.95, "Max Compounds Returned" = 200)
result2 <- result(job2)
job3 <- launchCMTool("pubchemID2SDF", result2)
result3 <- result(job3)
job4 <- launchCMTool("OpenBabel Descriptors", result1)
result4 <- result(job4)
