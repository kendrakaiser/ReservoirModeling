# Lucky Peak System Reservoir Modeling
Model of the Boise River Basin reservoir operations, coded to replicate the water control manual as closely as possible.

Monthly forecast data is used to determine what the storage is and bacn-calculate the minimum discharge on that day. This should be updated to be a daily forecast. Currently, it results in topping the dam or going beyond the diiscahrge limits. To deal with this the model enters an iterative loop where high flows are distributed over n days prior. The higher the n the better the discahrge looks, but also result in variable reservoir storage levels that aren't reasonable. 

