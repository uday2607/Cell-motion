from ellipses import *
import time

from ellipses import ellipse_intersection, polygon_area

def timing(f):
    def wrap(*args, **kwargs):
        times = []
        for i in range(1000):
            time1 = time.time()
            ret = f(*args, **kwargs)
            time2 = time.time()
            times.append(time2 - time1)
        print('{:s} function took {:.3f} ms'.format(f.__name__, np.mean(times)*1000.0))

        return ret
    return wrap

@timing
def test():
    intersc = ellipse_intersection(20, 20, 6, 3, 0, 18, 14, 6, 3, pi/4)
    area = polygon_area(intersc)
    return area

print(test())