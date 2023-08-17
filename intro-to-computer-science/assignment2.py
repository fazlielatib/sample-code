# Name : Latib, Mohamad Fazlie Bin
# Student Number : 30067991
# The application will draw graphs when mathematical expressions were entered.

from SimpleGraphics import *

# Set Background Colour
background ("white")

# Draw x and y axis line 
width = 800
half_width = 400
minimum_width = 0
height = 600
half_height = 300
minimum_height = 0
delta = 30

line(minimum_width,half_height, width,half_height)
line(half_width,minimum_height, half_width,height)

# Add ticks and scale to x axis
x_axis_tick_start = 295
x_axis_tick_end = 305
x_axis_scale_position = 315

for positive_x in range (1,14,1):
    line(half_width + delta * positive_x,x_axis_tick_start, half_width + delta * positive_x,x_axis_tick_end)
    setFont("Times", "10", "bold") 
    text(half_width + delta * positive_x, x_axis_scale_position, positive_x)

for negative_x in range (-13,0,1):
    line(half_width + delta * negative_x,x_axis_tick_start, half_width + delta * negative_x,x_axis_tick_end)
    setFont("Times", "10", "bold" ) 
    text(half_width + delta * negative_x, x_axis_scale_position, negative_x)
    
# Add ticks and scale to y axis
y_axis_tick_start = 395
y_axis_tick_end = 405
y_axis_scale_position = 385

for negative_y in range (-10,0,1):
        line(y_axis_tick_start, half_height - delta * negative_y, y_axis_tick_end, half_height - delta * negative_y)
        setFont("Times", "10", "bold" ) 
        text(y_axis_scale_position, half_height - delta * negative_y, negative_y)
    
for positive_y in range (1,11,1):
        line(y_axis_tick_start, half_height - delta * positive_y, y_axis_tick_end, half_height - delta * positive_y)
        setFont("Times", "10", "bold" ) 
        text(y_axis_scale_position, half_height - delta * positive_y, positive_y)

# Getting input from the user   
from math import *

expression = input("Enter an arithmetic expression (blank line to quit): \n y = ") 
x = -13
count = 0

# Making lines of different colors
while expression != " ": 
    if count%3 == 0:
        setColor("red")
    elif count%3 == 1:
        setColor("green")
    else:
        setColor("blue")
            
    while x >= -13 and x <= 13:
        first_y_coordinate = eval(expression)
        x = x + 0.1
        
        second_y_coordinate = eval(expression)
        x = x - 0.1
        x_shift = 0.1
        line(x * delta + half_width,first_y_coordinate * (-delta) + half_height, (x + x_shift) * delta + half_width,second_y_coordinate * (-delta) + half_height)

# Draw a point mark on the maximum and minimum points
        current_y_coordinate = eval(expression)
        x = x + 0.1
        next_y_coordinate = eval(expression)
        x = x - 0.2
        previous_y_coordinate = eval(expression)
        x = x + 0.1
        
        current_slope = next_y_coordinate - current_y_coordinate
        previous_slope = current_y_coordinate - previous_y_coordinate
        point_height = 5
        point_width = 5
        point_shift = 2
        
        if(current_slope * previous_slope <= 0):
            # Maximum point
            if (previous_slope > 0):
                setColor("purple")
                ellipse(x * delta + half_width - point_shift,first_y_coordinate * - delta + half_height - point_shift, point_width,point_height)
                
            # Minimum point
            else:
                setColor("orange")
                ellipse(x * delta + half_width - point_shift,first_y_coordinate * - delta + half_height - point_shift, point_width,point_height)
                
        if count%3 == 0:
            setColor("red")
        elif count%3 == 1:
            setColor("green")
        else:
            setColor("blue")
                            
        x = x + 0.1            
            
    expression = input("Enter an arithmetic expression (blank line to quit): \n  y = ")
    count = count + 1
    x = -13