import numpy as np

max_raidius = 10**6
              
num_points = 10**3

theta = np.random.default_rng().random(size=num_points, dtype=np.float32)
np.multiply(theta, 2*np.pi, out=theta)
radius = np.random.default_rng().random(size=num_points, dtype=np.float32)
np.sqrt(radius, out=radius)

print("Generated " + str(num_points) + " random point as radius and theta coordinates")

x = np.empty(num_points, dtype=np.float32)
np.cos(theta, out=x)
np.multiply(radius, x, out=x)
print("Coordinate conversion halfway finished")
np.sin(theta, out=theta)
np.multiply(radius, theta, out=radius)
y = radius
print("Coordinate conversion finished")

del theta, radius

x_min, y_min = x.min(), y.min()
np.subtract(x, x_min, out=x)
np.subtract(y, y_min, out=y)

print("Print first two elements to check if C code reads correctly(mostly because of endianity compatibility)")
print(x[:2])

f = open("round_" + str(num_points), "wb")
x = x.tobytes()
f.write(x)
del x
y = y.tobytes()
f.write(y)
del y
f.close()

print("All finished Correctly")
