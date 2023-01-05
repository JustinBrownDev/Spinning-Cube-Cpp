# Spinning-Cube-Cpp
This is a program that displays a spinning cube in the terminal

It does this by projecting the points of a cube onto a 2d plane, tracing the lines between the connected points, and printing that image to the terminal.
The shadow effect is done using the same process, by first projecting the points onto a plane which represents the floor surface, and then projecting those points onto the view window.

You can use the values in the dispconf struct to fine tune the program to work with your specific computer, if you set it to go too fast or too slow for your hardware you'll see a gross flickering effect.
