
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


bounding_box_2=np.array([[ -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85 ],[ -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85 ]])

bounding_box_3=np.array([[ -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, -55, -45, 45, 55, -75, -65, -55, 55, 65, 75, -75, 75, -85, -75, 75, 85, -95, -85, 85, 95, -95, 95, -95, 95, -95, 95, -95, 95, -95, 95, -95, 95, -95, 95, -95, 95, -95, -85, 85, 95, -85, -75, 75, 85, -75, 75, -75, -65, -55, 55, 65, 75, -55, -45, 45, 55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45],[ -95, -95, -95, -95, -95, -95, -95, -95, -95, -95, -85, -85, -85, -85, -75, -75, -75, -75, -75, -75, -65, -65, -55, -55, -55, -55, -45, -45, -45, -45, -35, -35, -25, -25, -15, -15, -5, -5, 5, 5, 15, 15, 25, 25, 35, 35, 45, 45, 45, 45, 55, 55, 55, 55, 65, 65, 75, 75, 75, 75, 75, 75, 85, 85, 85, 85, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95] ])

cell_rep=np.array([[ -35, -25, -15, -5, 5, 15, 25, 35, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, -35, -25, -15, -5, 5, 15, 25, 35 ],[ -75, -75, -75, -75, -75, -75, -75, -75, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 75, 75, 75, 75, 75, 75, 75, 75  ]])
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
plt.plot([ -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, -85, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85 ],[ -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -85, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85 ], 'wo')
for item in range(bounding_box_3.shape[1]):
    recy=mpatches.Rectangle((bounding_box_3[0][item]-5,bounding_box_3[1][item]-5), 10, 10, ec="red")  
    ax.add_patch(recy) 
    #recy=mpatches.Rectangle((bounding_box_1[item][0]-5,bounding_box_1[item][1]-5), 10, 10, ec="none")
    #ax.add_patch(recy)

circ0=mpatches.Circle((0, 0), 83.6505, ec="black", fc="none")
ax.plot([0],[0],'ko')
ax.add_patch(circ0)
circ1=mpatches.Circle((5, 20), 10, ec="green", fc="none")
circ2=mpatches.Circle((-7, 28), 10, ec="green", fc="none")
circ3=mpatches.Circle((40, 3.5), 10, ec="green", fc="none")
circ4=mpatches.Circle((100, -7.8), 10, ec="green", fc="none")
circ5=mpatches.Circle((17, -37), 10, ec="green", fc="none")
circ6=mpatches.Circle((-5, -62), 10, ec="green", fc="none")
circ7=mpatches.Circle((-20, 2), 10, ec="green", fc="none")
circ8=mpatches.Circle((-100, -6,0), 10, ec="green", fc="none")
ax.add_patch(circ1)
ax.add_patch(circ2)
ax.add_patch(circ3)
ax.add_patch(circ4)
ax.add_patch(circ5)
ax.add_patch(circ6)
ax.add_patch(circ7)
ax.add_patch(circ8)

#ax.plot([ -25, -15, -5, 5, 15, 25, -45, -35, -25, 25, 35, 45, -55, -45, 45, 55, -65, -55, 55, 65, -75, -65, 65, 75, -75, 75, -85, -75, 75, 85, -85, 85, -85, 85, -85, 85, -85, 85, -85, -75, 75, 85, -75, 75, -75, -65, 65, 75, -65, -55, 55, 65, -55, -45, 45, 55, -45, -35, -25, 25, 35, 45, -25, -15, -5, 5, 15, 25 ],[ -85, -85, -85, -85, -85, -85, -75, -75, -75, -75, -75, -75, -65, -65, -65, -65, -55, -55, -55, -55, -45, -45, -45, -45, -35, -35, -25, -25, -25, -25, -15, -15, -5, -5, 5, 5, 15, 15, 25, 25, 25, 25, 35, 35, 45, 45, 45, 45, 55, 55, 55, 55, 65, 65, 65, 65, 75, 75, 75, 75, 75, 75, 85, 85, 85, 85, 85, 85 ],'bo')
#ax.plot([ -25, -15, -5, 5, 15, 25, -45, 45, -55, 55, -65, 65, -75, 75, -85, 85, -85, 85, -85, 85, -85, 85, -85, 85, -85, 85, -75, 75, -65, 65, -55, 55, -45, 45, -25, -15, -5, 5, 15, 25, ],[ -85, -85, -85, -85, -85, -85, -75, -75, -65, -65, -55, -55, -45, -45, -25, -25, -15, -15, -5, -5, 5, 5, 15, 15, 25, 25, 45, 45, 55, 55, 65, 65, 75, 75, 85, 85, 85, 85, 85, 85,  ], 'go' )
#ax.plot([ -35, -25, -15, -5, 5, 15, 25, 35, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, -35, -25, -15, -5, 5, 15, 25, 35 ],[ -75, -75, -75, -75, -75, -75, -75, -75, -65, -65, -65, -65, -65, -65, -65, -65, -65, -65, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 75, 75, 75, 75, 75, 75, 75, 75  ],'yo')
#ax.grid()
#for item in range(cell_rep.shape[1]):
#    recy2=mpatches.Rectangle((cell_rep[0][item]-5,cell_rep[1][item]-5), 10, 10, ec="yellow",fc="black")  
#    ax.add_patch(recy2) 
    #recy=mpatches.Rectangle((bounding_box_1[item][0]-5,bounding_box_1[item][1]-5), 10, 10, ec="none")
    #ax.add_patch(recy)

plt.axis([-150, 150, -150, 150])
plt.show()

