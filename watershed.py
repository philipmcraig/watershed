#Code to write file to define watershed release points

import numpy as np
import itertools as it

f = open('watershed.txt','w')
f.write('File defining watershed release points \n \n')

# Start with straight line down one longitude

# latitude range = 30N to 30S
# longitude = 70W

lat = np.linspace(30,-30,101)
lon = np.linspace(70,70,101)

latlon = zip(lat,lon)

#for i in xrange(p.shape[0]):
#	f.write(p[i])

f.write('Latitudes  Longitudes \n \n')

np.savetxt(f, latlon, fmt='%9.5f')

