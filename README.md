This is the ReadMe file for the source code and datasets contained in the publication:
Holdo, R.M. and Onderdonk. D.A. Savanna tree abundance and spatial patterns are strongly associated
with river networks in Serengeti National Park, Tanzania. Landscape Ecology 2022.

Contact authors for all files:
Ricardo Holdo (rholdo@uga.edu)

File Final_Analysis.R
This file replicates all analyses and figures in the paper
The code draws from the following data files:

File "Final_30m_data.csv"
Raw tree cover data at 30-m resolution (after cubic interpolation). Columns are as follows:
ID: Plot ID
TC: Tree cover fraction
X: Eastings (Datum: WGS84, Projection: UTM 36S)
Y: Northings
driver: Distance to river (km)

File "Dist_to_riv_spatial_regressions_poly.csv"
Polynomial fits for regressions of TC vs. distance to river. Columns are as follows:
ID: Plot ID
AIC.null: AIC values for null (intercept) model
AIC.lin: AIC values for linear model
AIC.quad: AIC values for quadratic model
AIC.cub: AIC values for cubic model

File "Distance_to_riv_cutoff.csv"
Breakpoint for first segment of segmented regression of TC vs driver. Columns are as follows:
ID: Plot ID
psi: Breakpoint distance (km)
Slope: slope of first segment (1/km), based on gls regression of log TC vs. driver
p.val: P-value of gls regression

File "Aggregated_lacunarity_Data.csv"
Lacunarity data at three spatial scales. Columns are as follows:
ID: Plot ID
Dbin: Distance bin (m). Value represents midpoint in 50-m intervals from river.
Lac10: Normalized lacunarity at 10-m scale
Lac25: Normalized lacunarity at 25-m scale
Lac50: Normalized lacunarity at 50-m scale
dr: Distance to river (km)
dr2: dr ^ 2
dr3: dr ^ 3
