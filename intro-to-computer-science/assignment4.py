# Name : Latib, Mohamad Fazlie Bin
# Student Number : 30067991
# CPSC 217 Assignment 4 : Sankey's Diagram

from SimpleGraphics import *
import sys
import math

sys.argv = ["Assignment4.py", "California_Electricity.txt"]

# Set background to white
background("white")

# Set a list of different colors
color = [[255, 110, 199], [0,255,255], [0,255,0], [255,255,0], [50,205,50], [30,144,255], [255, 199, 0], [170, 255, 195], [60, 180, 75], [70, 240, 240], [0, 130, 200], [145, 30, 180]]
    
# Part 2

# Read the lines in the file excluding the title and the source name and insert it into a dictionary.
#
# Parameter: 
#   rest_of_lines: the lines in the file excluding the title and the source name.
#
# Return:
#   dictionary: dictionary consists of keys which are the destination name and values which are the destination amount
#
def dictionary_data(rest_of_lines):
    # Initiate an empty dictionary
    dictionary = {}
    # 
    for lines in rest_of_lines:
        # Seperate each line 
        file_line = lines.strip()
        # Split the line into two part by removing the comma
        segment = lines.split(",")
        # Indicate that the first part is the dictionary key or the destinations
        key = segment[0]
        # Indicate that the second part is the dictionary value or the amount of flow and set it as a floating number
        value = float(segment[1])
        dictionary[key] = value
    
    # Return the dictionary
    return dictionary

# Part 3,4 & 5

# Draw the source bar, the destination bar and the paths connecting them
# Put the destination name at the right side if the destination bars
# Insert the color in the source bar and the destination bars
# Insert the color gradient in the paths connecting the source bar and the destination bars
# Draw the black outline surrounding the built Sankey's diagram 
#
# Parameter:
#   dictionary: dictionary consists of keys which are the destination name and values which are the destination amount
#
# Returns:
#   source bar: colored source bar
#   destination bar: colored destination bars and their destination names
#   gap between destination bar: constant gap between destination bars
#   paths connecting source and destination bar: gradient color in the paths
#   outline: black outline for the whole drawn diagram
#   
def plot_all(dictionary):
    # Initiate amount of flow to 0
    total_flow = 0
    
    # Set the destination total as the length of the dictionary
    dest_total = len(dictionary)
    
    # Set the destination total height
    dest_total_height = 450
    # Set the gap in between the destination bar
    dest_gap = 10
    # Calculate the number of pixel available for the source and destination bar without the gap
    total_pixel = dest_total_height - (dest_total - 1) * dest_gap
    
    # Calculate total amount of flow by keep adding the value in the dictionary
    for key, value in dictionary.items():
        # Add value from dictionary to amount of flow
        total_flow = total_flow + value
    
    # Calculate the number of pixel available for each amount of flow
    pixel_per_flow = total_pixel / total_flow
    
    # Calculate the source bar height
    source_height = total_flow * pixel_per_flow
    
    # Set the source bar width
    source_width = 25
    
    # Calculate and set the y coordinate to start drawing the source bar
    source_y_left = 330 - source_height / 2
    
    # Set the x coordinate to start drawing the source bar
    source_x_left = 130
    
    # Initiate the colour index from the beginning of the list
    color_index = 0
    
    # Set the color of the source bar
    source_color = color[color_index]
    setColor(source_color[0],source_color[1],source_color[2])
    
    # Draw the source bar
    rect(source_x_left,source_y_left, source_width,source_height)
    
    # Use different color in the color lists for different destination bar
    color_index += 1
    
    # Set the outline of the source bar to be black
    setColor("black")
    
    # Draw the left side outline of the source bar
    line(source_x_left,source_y_left, source_x_left,(source_y_left - 1) + source_height)
    
    # Draw the top outline of the source bar
    line(source_x_left,source_y_left, source_x_left + source_width,source_y_left)
    
    # Draw the bottom outline of the source bar
    line(source_x_left,(source_y_left - 1) + source_height, source_x_left + source_width,(source_y_left - 1) + source_height)
    
    # Set the x coordinate to start drawing the destination bar
    dest_x_left = 645
    
    # Set the y coordinate to start drawing the destination bar
    dest_y_left = 330 - dest_total_height / 2

    # Iniatiate the destination total height to 0
    dest_sum_height = 0
    
    for key, value in dictionary.items():
        # Calculate the height for each destination bar
        dest_height = value * pixel_per_flow
        
        # Set the width of the destination bar
        dest_width = 25
        
        # Set the x coordinate to put the destination name
        dest_x_text = 680
        
        # Calculate and set the coordinate to put the destination name
        dest_y_text = dest_y_left + dest_height / 2
        
        # Set the colour of the destination bar
        dest_color = color[color_index]
        setColor(dest_color[0],dest_color[1],dest_color[2])
        
        # Draw the destination bar
        rect(dest_x_left,dest_y_left, dest_width,dest_height)
        
        # Set the color of the destination name
        setColor("black")
        
        # Draw the name of the destination
        text(dest_x_text,dest_y_text, key, "w")
        
        # Set the outline of the destination bar to be black
        setColor("black")
        
        # Draw the right side outline of the destination bar
        line((dest_x_left - 1) + dest_width,dest_y_left, (dest_x_left - 1) + dest_width,dest_y_left + dest_height) 
        
        # Draw the top outline of the destinantion bar
        line(dest_x_left,dest_y_left, (dest_x_left - 1) + dest_width,dest_y_left)  
        
        # Draw the bottom outline of the destinantion bar
        line(dest_x_left,(dest_y_left - 1) + dest_height, (dest_x_left - 1) + dest_width,(dest_y_left - 1) + dest_height)
        
        # Calculate and set the y coordinate of the source bar connecting to the destination bar
        source_y_right = source_y_left + dest_sum_height
        #polygon(source_x_left + source_width,source_y_right, 650,dest_y_left, 650,(dest_y_left - 1) + dest_height, source_x_left + source_width,(source_y_right - 1) + dest_height)


        for x in range(source_x_left + source_width, dest_x_left):
            # Use ratio with 0 representing being at the source, increasing to 1 when the destination is reached
            p = (x - (source_x_left + source_width)) / (dest_x_left - (source_x_left + source_width))
            
            # Set the path from the source bar to the destination bar as a curve
            p = p * math.pi - math.pi / 2
            p = (math.sin(p) + 1) / 2 
            
            # Calculate the r, g and b value according to the ratio at certain x coordinate
            r = source_color[0] * (1 - p) + dest_color[0] * p
            g = source_color[1] * (1 - p) + dest_color[1] * p
            b = source_color[2] * (1 - p) + dest_color[2] * p
            
            # Calculate the y coordinate for the top path following the ratio at x coordinate
            y_top = source_y_right * (1 - p) + dest_y_left * p
            
            # # Calculate the y coordinate for the bottom path following the ratio at x coordinate
            y_bottom = ((source_y_right - 1) + dest_height) * (1 - p) + ((dest_y_left - 1) + dest_height) * p
            
            # Set color of the line path
            setColor(r,g,b)
            
            # Draw a line vertically connecting the top and bottom y coordinate
            line(x,y_top, x,y_bottom)
            
            # Set the outline of the path to be black
            setColor("black")
            
            # Draw the outline of the path at top and bottom of the path
            line(x,y_top, x,y_top + 1)
            line(x,y_bottom, x,y_bottom + 1)      
            
        # Calculate new y coordinate of the destination bar
        dest_y_left += (dest_gap + dest_height)
        
        # Calculate the new height total of the destination bar
        dest_sum_height += dest_height
        
        # Use the next color position in the color list
        color_index += 1
    
    # Print to check the total amount of flow
    print(total_flow)
    
# Part 1

# Read, display the title, and display the label for the source on the left side of the graph
# If no command line argument is provided then program should prompt the user and read the name of the file
# Program should quit with an appropriate error message if the user provides more than one file name on the command line, or if the file specified by the user (either on the command line or read using the input function) doesn’t exist
# Call the other functions to draw the diagram
#
# Parameter: 
#   (None)
#
# Return:
#   title: the title is drawn at the top centre of the background
#   source title: the source title is drawn at the left side of the diagram
#   the Sankey's Diagram: the drawn diagram if there's no error in the beginning
#
def main():
    if len(sys.argv) == 1:
        print("No file name was provided.")
        file_name = input("Enter a file name: ")
        
    elif len(sys.argv) == 2:
        file_name = sys.argv[1]
        
    else:       
        print("Too many command lines.")
        quit()
    
    try:    
        inf = open(file_name, "r")
    
    except FileNotFoundError:
        print("Invalid file is entered.")
        quit()
    
    # Read the first line of the file   
    file_line = inf.readline()
    
    # Set the x and y coordinate for the title
    text_x = 400    
    text_y = 65  
    
    # Draw the first line of the file as the title
    setFont("Cambria", "20", "bold")
    text(text_x,text_y, file_line)
    
    # Read the second line of the file
    file_line = inf.readline()
    
    # Set x and y coordinate for the source name
    source_x_text = 115
    source_y_text = 330
    
    # Draw the second line of the file as the source name
    setFont("Cambria", "12")
    text(source_x_text,source_y_text, file_line, "e")
        
    # Read the rest of the lines in the file 
    rest_of_lines = inf.readlines()        
        
    # Call the dictionary function
    dictionary = dictionary_data(rest_of_lines)
    
    # Print the dictionary to verify that all of the data was loaded and stored into the dictionary correctly
    print(dictionary)    
    
    # Call the plot function
    plot_all(dictionary) 

  
main()

