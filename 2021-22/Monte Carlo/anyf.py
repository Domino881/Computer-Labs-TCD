# Python laboratory - week 7
# Monte Carlo method (task 3)
# author - Dominik Kuczynski

# For this task, I have taken a flexible approach - this means that 
# (as far as I know), basically any function as well as any interval may be 
# used to approximate with the monte carlo method.

import matplotlib.pyplot as plt
import numpy as np

##############################################################################
def f(x):
    return (x<=10)*20

xlower_bound = 0
xupper_bound = 20

times = 2000
expected_value = 0.0

# Above are values to be input by the user: the function f of which an interval
# is to be approximated, the bounds of the interval as well as the number of
# samples and expected value.
##############################################################################

def sign(x):
    return x/np.abs(x)
# this is a function returning -1 if x<0 and 1 if x>0

def cummean(arr):
    return np.cumsum(arr) / np.arange(1, len(arr)+1, dtype=np.float)
def cumstd(arr):
    return np.sqrt(cummean(arr**2) - cummean(arr)**2)
# these are two functions described in the previous task, returning the 
# cummulative mean and standard deviation of an array.


x_plot = np.linspace(xlower_bound, xupper_bound, 
                    int(xupper_bound-xlower_bound+1)*100)
y_plot = f(x_plot)
# x_plot - 100 points spaced evenly between the bounds of the integral
# y_plot - the corresponding values of the function f


fig1 = plt.figure()
fig1.set_size_inches(5, 5)
# a square figure is created


ylower_bound = np.min(y_plot)
yupper_bound = np.max(y_plot)
if(ylower_bound > 0 and yupper_bound > 0):
    ylower_bound = 0
if(ylower_bound < 0 and yupper_bound < 0):
    yupper_bound = 0
# bounds for the y coordinates are created from the maximum and minimum
# values of the funtion on the interval, making sure to include
# the x-axis i.e. y=0.

x_random = np.random.rand(times)*(xupper_bound-xlower_bound) + xlower_bound
y_random = np.random.rand(times)*(yupper_bound-ylower_bound) + ylower_bound
# "times" random points are generated with random coordinates starting at 
# the lower bounds and up to the upper bounds.

plt.plot(x_random, y_random, 'o', markersize=1)
plt.plot(x_plot, y_plot)
plt.xlabel("x")
plt.ylabel("y")
# the random points as well as the function are plotted


box_size = (xupper_bound-xlower_bound)*(yupper_bound-ylower_bound)
# to correctly approximate the area below the graph,
# we need to know the overall area of the square, i.e. box_size.

below = sign(y_random) * (np.abs(y_random) <= np.abs(f(x_random)))
# a value indicating whether a point is closer to zero than a value of f
# if not closer to zero than f         --> 0
# if closer to zero than f & above y=0 --> 1
# if closer to zero than f & below y=0 --> -1
# this allows for calculating the area below y=0 as negative

below_ratio = np.sum(below) / times
# a ratio of points included in integral to all scattered points
approximation = box_size * below_ratio
print("The estimated value of the integral is: ", approximation)
# from the ratio, the area below the graph is calculated and printed


cummulative_ratio = np.cumsum(below) / np.arange(1, times+1, 
                                                 dtype=np.float)
cummulative_integral = box_size * cummulative_ratio
# the cummulative ratios of points below the graph to all points are 
# calculated, and from them, the respective approximated areas below the graph

fig2 = plt.figure()
fig2.set_size_inches(5,5)
# another square figure is created

achieved, = plt.plot(cummulative_integral)
expected, = plt.plot(np.repeat(expected_value, times))
# the achieved approximations and the expected value are plotted against
# sample size

plt.ylim(expected_value-cumstd(cummulative_integral)[int(0.05*times)],
         expected_value+cumstd(cummulative_integral)[int(0.05*times)])
# to make the plot legible for an arbitrary function, its y-axis is scaled
# from expected_value-delta to expected_value+delta, with delta being the
# standard deviation of the first 5% of approximations.

plt.legend([achieved, expected], ["Approximation", "Expected value"])
plt.xlabel("Sample size")
plt.title("Approximation vs Sample Size")

fig3 = plt.figure()
fig3.set_size_inches(5, 5)
# a third square figure is created

stdevs, = plt.plot(cumstd(cummulative_ratio))
plt.legend([stdevs], ["Standard deviation"])
plt.xlabel("Sample size")
plt.title("Standard Deviation vs Sample Size")
# the cummulative standard deviation is plotted against sample size
