setwd("/code/rgeomstats")
library(reticulate)
use_python('/usr/local/bin/python3', required = T)
py_config()

gs <- import("geomstats")
dimension = 2L

sphere <- gs$hypersphere$Hypersphere(dimension=dimension)


data = sphere$random_uniform(10L)

sphere_visu = gs$visualization$Sphere()
sphere_visu$add_points(data)

gs$visualization$plot(data, 'S2')

mean = sphere$metric$mean(data)
tangent.pca = sphere$metric$tangent_pca(data)

regression.line = sphere$metric$geodesic(initial_point = mean, initial_tangent_vec = tangent.pca[2][1])
point = regression.line(1)

gs$visualization$plot(regression.line)
py_install("geomstats", envname = "python3env")
gs <- import("geomstats")

use_virtualenv("python3env", required = TRUE)
virtualenv_install("python3env", "geomstats")
