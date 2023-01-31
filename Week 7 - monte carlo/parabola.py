# Python laboratory - week 7
# Monte Carlo method (task 2)
# author - Dominik Kuczynski

import matplotlib.pyplot as plt
import numpy as np

parabola_x = np.linspace(0, 1, 100)
parabola_y = parabola_x**2
# lists of x and y coordinates:
# parabola_x -> 100 evenly spaced numbers in (0,1)
# parabola_y -> respective y-coordinates of points on graph of f(x)=x**2


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
plt.plot(parabola_x, parabola_y)
plt.xlabel("x")
plt.ylabel("y")
# the points are plotted on a new graph

test_function = y <= x**2
# a boolean determining whether each point is below the graph.

below_ratio = float(np.sum(test_function))/float(times)
# the ratio of points below the graph to all points

print('The estimated value of the integral is: ', below_ratio)

cumulative_below = np.cumsum(test_function)
# these are the numbers of points below the graph, depending on the number
# of points scattered. 
# cumulative_below[x] = how many out of x points landed below the graph

cumulative_ratios = (cumulative_below / 
                     np.arange(1, times+1, dtype=np.float))
# the ratios of points below the graph to all points, depending on number
# of points.
# As the area of the square in which all points are scattered is equal to 1,
# these are also the approximated integrals.

plt.figure()
approx_integrals, = plt.plot(cumulative_ratios)
integral, = plt.plot(np.repeat(1.0/3.0, times))

plt.ylim(0.2, 0.43)
plt.xlabel("Sample size")
plt.title("Approximation vs Sample Size")
plt.legend([approx_integrals, integral], ["Approximation", "Exact value"])
# the predicted values of the integral are plotted against the total number of 
# points scattered. The plot is scaled from 0.2 to 0.43 on the y-axis to better
# show the improvement of predictions.

def cummean(arr):
    return np.cumsum(arr) / np.arange(1, len(arr)+1, dtype=np.float)
# a function for calculating the cumulative means of elements in an array up to
# every index

def cumstd(arr):
    return np.sqrt(cummean(arr**2) - cummean(arr)**2)
# a function for calculating the cumulative standard deviation of elements in
# an array up to every index

plt.figure()
stdevs, = plt.plot(cumstd(cumulative_ratios))
plt.legend([stdevs], ["Standard deviation"])
plt.xlabel("Sample size")
plt.title("Standard Deviation vs Sample Size")
# the cummulative standard deviation is plotted against the number of samples