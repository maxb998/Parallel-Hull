import numpy as np
import sys

def main():

    num_string = sys.argv[1]
    num_points = eval(num_string)
    print("num_points = " + str(num_points))
    fname = num_string.replace("**","e").replace("*", "x")
    fname = "round_" + fname
    print("fname = " + fname)
    
    max_radius = np.power(num_points, 1/4)

    theta = np.random.default_rng().random(size=num_points, dtype=np.float32)
    np.multiply(theta, 2*np.pi, out=theta)
    radius = np.random.default_rng().random(size=num_points, dtype=np.float32)
    np.sqrt(radius, out=radius)
    np.multiply(radius, max_radius, out=radius)

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

    np.add(x, max_radius, out=x)
    np.add(y, max_radius, out=y)

    print("Print first two elements to check if C code reads correctly(mostly because of endianity compatibility)")
    print(x[:2])

    f = open(fname, "wb")
    x = x.tobytes()
    f.write(x)
    del x
    y = y.tobytes()
    f.write(y)
    del y
    f.close()

    print("All finished Correctly")


if __name__ == "__main__":
    main()