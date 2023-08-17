# Name : Latib, Mohamad Fazlie Bin
# Student Number : 30067991

# The purpose of my program is to create a smiley that looks like a girl when the program is executed.

from SimpleGraphics import*

background("light sky blue")

x = int(input("Enter X: "))
y = int(input("Enter Y: "))

# Face
setOutline("yellow")
setFill("yellow")
ellipse(x-200,y-200, 400,400)

# Cheeks
setOutline("gold")
setFill("gold")
pieSlice(x-170,y+10, 50,50, 0,359)

setOutline("gold")
setFill("gold")
pieSlice(x+120,y+10, 50,50, 0,359)

# White Eyes
setOutline("white")
setFill("white")
ellipse(x-100,y-100, 50,100)

setOutline("white")
setFill("white")
ellipse(x+50,y-100, 50,100)

# Black Eyes
setOutline("black")
setFill("black")
ellipse(x-90,y-60, 30,40)

setOutline("black")
setFill("black")
ellipse(x+60,y-60, 30,40)

# Eyebrows
setOutline("black")
setFill("black")
polygon(x-105,y-120, x-95,y-130, x-30,y-120, x-40,y-110)

setOutline("black")
setFill("black")
polygon(x+105,y-120, x+95,y-130, x+30,y-120, x+40,y-110)

# Hair
setOutline("pink")
setFill("pink")
blob(x-200,y+150, x-205,y, x-180,y-200, x,y-240, x+180,y-200, x+205,y, x+200,y+150, x+210,y+305, x+220,y+300, x+220,y, x,y-260, x-220,y, x-220,y+300, x-210,y+305)

# Hairbands
setOutline("red")
setFill("red")
rect(x-210,y+28, 15,15)

setOutline("red")
setFill("red")
rect(x+195,y+28, 15,15)

# Flowers
setOutline("blue")
setFill("blue")
ellipse(x-143,y-161, 30,30)

setOutline("blue")
setFill("blue")
ellipse(x-130,y-180, 30,30)

setOutline("blue")
setFill("blue")
ellipse(x-145,y-200, 30,30)

setOutline("blue")
setFill("blue")
ellipse(x-165,y-190, 30,30)

setOutline("blue")
setFill("blue")
ellipse(x-165,y-168, 30,30)

setOutline("light blue")
setFill("light blue")
ellipse(x-152,y-182, 34,34)

# Mouth
setOutline("wheat")
setFill("wheat")
blob(x-140,y+90, x,y+170, x+140,y+90, x,y+90)

# Teeth
setOutline("white")
setFill("white")
blob(x-130,y+90, x,y+120, x+130,y+90, x,y+90)

# Tongue
setOutline("salmon")
setFill("salmon")
blob(x-80,y+134, x,y+154, x+80,y+134, x,y+114)