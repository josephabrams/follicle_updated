
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

bounding_box_1=np.array([[ -85, -85],[ -75, -85],[ -65, -85],[ -55, -85],[ -45, -85],[ -35, -85],[ -25, -85],[ -15, -85],[ -5, -85],[ 5, -85],[ 15, -85],[ 25, -85],[ 35, -85],[ 45, -85],[ 55, -85],[ 65, -85],[ 75, -85],[ -85, -75],[ -75, -75],[ -65, -75],[ -55, -75],[ -45, -75],[ -35, -75],[ -25, -75],[ -15, -75],[ -5, -75],[ 5, -75],[ 15, -75],[ 25, -75],[ 35, -75],[ 45, -75],[ 55, -75],[ 65, -75],[ 75, -75],[ -85, -65],[ -75, -65],[ -65, -65],[ -55, -65],[ -45, -65],[ -35, -65],[ -25, -65],[ -15, -65],[ -5, -65],[ 5, -65],[ 15, -65],[ 25, -65],[ 35, -65],[ 45, -65],[ 55, -65],[ 65, -65],[ 75, -65],[ -85, -55],[ -75, -55],[ -65, -55],[ -55, -55],[ -45, -55],[ -35, -55],[ -25, -55],[ -15, -55],[ -5, -55],[ 5, -55],[ 15, -55],[ 25, -55],[ 35, -55],[ 45, -55],[ 55, -55],[ 65, -55],[ 75, -55],[ -85, -45],[ -75, -45],[ -65, -45],[ -55, -45],[ -45, -45],[ -35, -45],[ -25, -45],[ -15, -45],[ -5, -45],[ 5, -45],[ 15, -45],[ 25, -45],[ 35, -45],[ 45, -45],[ 55, -45],[ 65, -45],[ 75, -45],[ -85, -35],[ -75, -35],[ -65, -35],[ -55, -35],[ -45, -35],[ -35, -35],[ -25, -35],[ -15, -35],[ -5, -35],[ 5, -35],[ 15, -35],[ 25, -35],[ 35, -35],[ 45, -35],[ 55, -35],[ 65, -35],[ 75, -35],[ -85, -25],[ -75, -25],[ -65, -25],[ -55, -25],[ -45, -25],[ -35, -25],[ -25, -25],[ -15, -25],[ -5, -25],[ 5, -25],[ 15, -25],[ 25, -25],[ 35, -25],[ 45, -25],[ 55, -25],[ 65, -25],[ 75, -25],[ -85, -15],[ -75, -15],[ -65, -15],[ -55, -15],[ -45, -15],[ -35, -15],[ -25, -15],[ -15, -15],[ -5, -15],[ 5, -15],[ 15, -15],[ 25, -15],[ 35, -15],[ 45, -15],[ 55, -15],[ 65, -15],[ 75, -15],[ -85, -5],[ -75, -5],[ -65, -5],[ -55, -5],[ -45, -5],[ -35, -5],[ -25, -5],[ -15, -5],[ -5, -5],[ 5, -5],[ 15, -5],[ 25, -5],[ 35, -5],[ 45, -5],[ 55, -5],[ 65, -5],[ 75, -5],[ -85, 5],[ -75, 5],[ -65, 5],[ -55, 5],[ -45, 5],[ -35, 5],[ -25, 5],[ -15, 5],[ -5, 5],[ 5, 5],[ 15, 5],[ 25, 5],[ 35, 5],[ 45, 5],[ 55, 5],[ 65, 5],[ 75, 5],[ -85, 15],[ -75, 15],[ -65, 15],[ -55, 15],[ -45, 15],[ -35, 15],[ -25, 15],[ -15, 15],[ -5, 15],[ 5, 15],[ 15, 15],[ 25, 15],[ 35, 15],[ 45, 15],[ 55, 15],[ 65, 15],[ 75, 15],[ -85, 25],[ -75, 25],[ -65, 25],[ -55, 25],[ -45, 25],[ -35, 25],[ -25, 25],[ -15, 25],[ -5, 25],[ 5, 25],[ 15, 25],[ 25, 25],[ 35, 25],[ 45, 25],[ 55, 25],[ 65, 25],[ 75, 25],[ -85, 35],[ -75, 35],[ -65, 35],[ -55, 35],[ -45, 35],[ -35, 35],[ -25, 35],[ -15, 35],[ -5, 35],[ 5, 35],[ 15, 35],[ 25, 35],[ 35, 35],[ 45, 35],[ 55, 35],[ 65, 35],[ 75, 35],[ -85, 45],[ -75, 45],[ -65, 45],[ -55, 45],[ -45, 45],[ -35, 45],[ -25, 45],[ -15, 45],[ -5, 45],[ 5, 45],[ 15, 45],[ 25, 45],[ 35, 45],[ 45, 45],[ 55, 45],[ 65, 45],[ 75, 45],[ -85, 55],[ -75, 55],[ -65, 55],[ -55, 55],[ -45, 55],[ -35, 55],[ -25, 55],[ -15, 55],[ -5, 55],[ 5, 55],[ 15, 55],[ 25, 55],[ 35, 55],[ 45, 55],[ 55, 55],[ 65, 55],[ 75, 55],[ -85, 65],[ -75, 65],[ -65, 65],[ -55, 65],[ -45, 65],[ -35, 65],[ -25, 65],[ -15, 65],[ -5, 65],[ 5, 65],[ 15, 65],[ 25, 65],[ 35, 65],[ 45, 65],[ 55, 65],[ 65, 65],[ 75, 65],[ -85, 75],[ -75, 75],[ -65, 75],[ -55, 75],[ -45, 75],[ -35, 75],[ -25, 75],[ -15, 75],[ -5, 75],[ 5, 75],[ 15, 75],[ 25, 75],[ 35, 75],[ 45, 75],[ 55, 75],[ 65, 75],[ 75, 75]])
bounding_box_2=([ -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75 ],[ -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75 ] )
a=np.diag(range(15))
#print(a)
#print(len(bounding_box))
#print(bounding_box.size)
#print(bounding_box.shape)
#print(bounding_box.reshape(2,289))
#plt.matshow(bounding_box.reshape(2,289))
#(nrows, ncols) = bounding_box.shape
#pixel_plot=plt.figure() #make plot
#plt.title("bounding box")
fig,ax=plt.subplots()
plt.plot([ -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75 ],[ -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75 ], 'ro')
for item in range(289):
    recy=mpatches.Rectangle((bounding_box_1[item][0]-5,bounding_box_1[item][1]-5), 10, 10, ec="red")  
    ax.add_patch(recy) 
    #recy=mpatches.Rectangle((bounding_box_1[item][0]-5,bounding_box_1[item][1]-5), 10, 10, ec="none")
    #ax.add_patch(recy)

circ=mpatches.Circle((0, 0), 84, ec="black", fc="none")
plt.plot([0],[0],'ro',color='black','ro')
ax.add_patch(circ)
#ax.grid()
plt.axis([-100, 100, -100, 100])
plt.show()

