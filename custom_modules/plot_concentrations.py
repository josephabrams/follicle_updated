
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import csv
def fit_voxel_plot_from_file(filename,time,x_column,y_column,z_column,concentration, lines_to_remove_from_bottom, save_figure=False):
    data_arr=np.genfromtxt(filename,delimiter=",")
    number_of_rows=16
    x_array=np.zeros(number_of_rows,dtype=float)
    y_array=np.zeros(number_of_rows,dtype=float)
    z_array=np.zeros(number_of_rows,dtype=float)
    concentration_array=np.zeros(number_of_rows,dtype=float)
    for i in range(number_of_rows):
        x_array[i]= data_arr[i][x_column]
        y_array[i]= data_arr[i][y_column]
        z_array[i]= data_arr[i][z_column]
        concentration_array[i]= data_arr[i][concentration]
    print(x_array) 
    fig, ax = plt.subplots()
    title="concentration fit from file: "+ filename.strip(".csv")+"/"
    print("Plotting ", title)
    X, Y= np.meshgrid(x_array,y_array)
    Z = X^2+Y^2#concentration_array
    ax.contour(X,Y,Z)
    plt.show()
#    if save_figure==True:
#        plt.savefig(filename.strip(".csv")+"_plot")
#
def main():
    fit_voxel_plot_from_file("./virial_eq/concentrations.csv",0,1,2,3,4,0)
if __name__ == '__main__':
    main()
