import casatasks

imagename = "raw"

d = casatasks.imstat(imagename+".image")
print(d["sigma"])
print(d["rms"], "Jy/beam")

# [0.02557265] Jy/beam for single channel and robust = 0.0