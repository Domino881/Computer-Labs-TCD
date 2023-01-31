# Python laboratory - week 7
# Monte Carlo method (task 1)
# author - Dominik Kuczynski

import matplotlib.pyplot as plt
import numpy as np

circle_x = np.linspace(0, 1, 100)
circle_y = np.sqrt(1-circle_x**2)
# lists of x and y coordinates:
# circle_x -> 100 evenly spaced numbers in (0,1)
# circle_y -> respective y-coordinates of points on the circle with center
#             at (0,0) and radius 1.


fig1 = plt.figure()
fig1.set_size_inches(5, 5)
# a square figure is created

times = 2000
# a variable for determining the amount of points scattered on the graph and
# used for approximation

x = np.random.rand(times)
y = np.random.rand(times)
# x,y - lists of random coordinates of the scattered points

plt.plot(x, y, 'o', markersize=1)
plt.plot(circle_x, circle_y)
plt.xlabel("x")
plt.ylabel("y")
# the points are plotted on a new graph

dist_squared = x**2 + y**2 
incircle = dist_squared <= 1
# a boolean determining whether each point is inside the circle.
# This is equivalent to np.sqrt(x**2+y**2)<=1 but faster computationally.

incircle_ratio = float(np.sum(incircle))/float(times)
# the ratio of points in the circle to all points

pi = incircle_ratio * 4
# as we're only looking at a quarter of the circle, the achieved value has
# to be multiplied by 4.
print('The estimated value of pi is: ', pi)

cumulative_incircle = np.cumsum(incircle)
# these are the numbers of points inside the circle, depending on the number
# of points scattered. 
# cumulative_incircle[x] = how many out of x points landed in the circle

cumulative_ratios = (cumulative_incircle / 
                     np.arange(1, times+1, dtype=np.float))
# the ratios of points inside the circle to all points, depending on number
# of points.

pis = cumulative_ratios * 4
# the respective predicted values of pi

plt.figure()
approx_pis, = plt.plot(pis)
pi, = plt.plot(np.repeat(np.pi, times))

plt.ylim(3.1, 3.3)
plt.xlabel("Sample size")
plt.title("Approximation vs Sample Size")
plt.legend([approx_pis, pi], ["Approximation", "Exact value"])
# the predicted values of pi are plotted against the total number of points
# scattered. The plot is scaled from 3.1 to 3.3 on the y-axis to better show
# the improvement of predictions.

def cummean(arr):
    return np.cumsum(arr) / np.arange(1, len(arr)+1, dtype=np.float)
# a function for calculating the cumulative means of elements in an array up to
# every index

def cumstd(arr):
    return np.sqrt(cummean(arr**2) - cummean(arr)**2)
# a function for calculating the cumulative standard deviation of elements in
# an array up to every index

plt.figure()
stdevs, = plt.plot(cumstd(pis))
plt.legend([stdevs], ["Standard deviation"])
plt.xlabel("Sample size")
plt.title("Standard Deviation vs Sample Size")
# the cummulative standard deviation is plotted against the number of samples