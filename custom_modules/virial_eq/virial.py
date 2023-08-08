# Virial and Polynomial fitting CPA solutions
# Author: Joseph E.S. Abrams
#The following is code for converting between molarity, molality, and osmolality of solutions for cryopreservation, the virial equations and coefficients come from 
# Prickett, R. C., Elliott, J. A. W. & McGann, L. E. Application of the osmotic virial equation in cryobiology. Cryobiology 60, 30–42 (2010)
# and Elliott, J. A. W., Prickett, R. C., Elmoazzen, H. Y., Porter, K. R. & McGann, L. E. A multisolute osmotic virial equation for solutions of interest in biology. J. Phys. Chem. B 111, 1775–1785 (2007)
# if you do use or adopt this code please cite the github citation 

import numpy as np
import matplotlib.pyplot as plt
import csv
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
#for one electrolyte m2 and one cpa m3 from Prickett, Elliott and McGann
#ternary solutions osmolality is calculated with the virial_eq_electrolyte using binary molalities
def virial_eq_electrolyte(kdiss2,kdiss3,m2,m3,B2,B3,C2,C3):
    # by convention m1 kdiss1 etc is for the solvent    
    m2=kdiss2*m2
    m3=kdiss3*m3
    osmolality=m2+m3+B2*((m2)**2)+B3*(m3**2)+(B2+B3)*m2*m3+C2*((m2)**3)+C3*(m3**3)+3*(((C2**2)*C3)**(1/3))*((m2)**2)*m3+3*((C2*(C3**2))**(1/3))*m2*(m3**2)
    return osmolality
def single_electrolyte_virial(kdiss,m2,B2,C2):
    osmolality=kdiss*m2+B2*((kdiss*m2)**2)+C2*((kdiss*m2)**3)
    return osmolality
#volume of the holding media
def holding_media_volume(m1,m2):
    # mole
    pass
def percent_volume_EG(m1,m2,m3,psi):
    pass
# brute force guess with constraints
#a=virial_eq_electrolyte(1.678,.125,2,0.044,0.037,0,-0.001)

virial_dictionary ={
        "NaCl":{
            "kdiss": 1.678,
            "B": 0.044,
            "C": 0.00
        },
        "KCI":{
            "kdiss": 1.772,
            "B": 0.00,
            "C": 0.00

        },
        "Me2SO":{
            "kdiss": 1,
            "B": 0.108,
            "C": 0.00
        },
        "Glycerol":{
            "kdiss": 1.00,
            "B": 0.023,
            "C": 0.00
        },
        "GLY":{
            "kdiss": 1.00,
            "B": 0.023,
            "C": 0.00
        },
        "PG":{
            "kdiss": 1.00,
            "B": 0.039,
            "C": 0.00
        },
        "EG":{
            "kdiss": 1.00,
            "B": 0.037,
            "C": -0.001
        },
        "Methanol":{
            "kdiss": 1.00,
            "B": 0.004,
            "C": 0.00
        },
        "Mannitol":{
            "kdiss": 1.00,
            "B": 0.00,
            "C": 0.00
        },
        "Sucrose":{
            "kdiss": 1.00,
            "B": 0.125,
            "C": 0.00
        },
        "Dextrose":{
            "kdiss": 1.00,
            "B": 0.044,
            "C": 0.00,
        },
        "Trehalose":{
            "kdiss": 1.00,
            "B": -0.394,
            "C": 0.388
        },
    #Note the following are for very small molalities less than 0.0123,0.00972 and 0.0195 respectively   
        "Hemoglobin":{
            "kdiss": 1.00,
            "B": 49.252,
            "C": 3.07e4
        },
        "BSA":{
            "kdiss": 1.00,
            "B": 3.70e2,
            "C": 1.60e5
        },
        "OVL":{
            "kdiss": 1.00,
            "B": 3.78e2,
            "C": 0.00
        }
    }
### put in the target osmolality and grid search molalities to reach it within 0.0001, solute1 and solute2 should refer the names in the above virial_dictionary    
def grid_search_ternary_2D(known_molality_for_solute1,osmolality,m1_start_value,m1_end_value,m2_start_value,m2_end_value,solute1,solute2):
    for i in np.linspace(m1_start_value,m1_end_value,20000):
        for j in np.linspace(m2_start_value,m2_end_value,20000):
            test_m1=i
            test_m2=j
            kdiss1=virial_dictionary[solute1]["kdiss"]
            kdiss2=virial_dictionary[solute2]["kdiss"]
            B1=virial_dictionary[solute1]["B"]
            B2=virial_dictionary[solute2]["B"]
            C1=virial_dictionary[solute1]["C"]
            C2=virial_dictionary[solute2]["C"]
            solution_osmolality_guess=virial_eq_electrolyte(kdiss1,kdiss2,test_m1,test_m2,B1,B2,C1,C2)
            if abs(known_molality_for_solute1-test_m1 )<=0.00001 and abs(osmolality-solution_osmolality_guess)<=0.00001:
                print(solution_osmolality_guess, " Osmole/kg produced by ", test_m1, " mol/kg ",solute1,", and ", test_m2, " mol/kg ", solute2)
                return {test_m1,test_m2} 
def grid_search_binary(kdiss,osmolality,virial_coeff_B,virial_coeff_C, start_value, end_value):
    value=-1.00
    for i in np.linspace(start_value,end_value,200000):
        test_m=i
        guess=single_electrolyte_virial(kdiss,test_m,virial_coeff_B,virial_coeff_C)
       # print(" guess: ", guess,", test value: ",test_m)
        if abs(osmolality-guess)<=0.00001:
            value=test_m
            osmo_out=guess
    if value<0:
        print("failed to converge to a solution")
        return 0.00
    else:
        print("media at ", osmo_out, "osmole/kg has ", value, " mol/kg of single solute")
        return value 
        #note may have more than one solution depending on precision asked for
            
#combining the two if you have a known osmolality for one of the solutes may need to manually adjust intervals to get convergence
def grid_search_ternary_1D(osmolality1,solution_osmolality,m2_start_value,m2_end_value,solute1,solute2):
    value=-1
    kdiss1=virial_dictionary[solute1]["kdiss"]
    kdiss2=virial_dictionary[solute2]["kdiss"]
    B1=virial_dictionary[solute1]["B"]
    B2=virial_dictionary[solute2]["B"]
    C1=virial_dictionary[solute1]["C"]
    C2=virial_dictionary[solute2]["C"]
    wide_start=osmolality1/2
    wide_end=osmolality1+wide_start
    test_m1=grid_search_binary(kdiss1,osmolality1,B1,C1,wide_start,wide_end)
    for j in np.linspace(m2_start_value,m2_end_value,200000):
        test_m2=j
        solution_osmolality_guess=virial_eq_electrolyte(kdiss1,kdiss2,test_m1,test_m2,B1,B2,C1,C2)
        if abs(solution_osmolality-solution_osmolality_guess)<=0.0001:
            print(solution_osmolality_guess, " Osmole/kg produced by ", test_m1, " mol/kg ",solute1,", and ", test_m2, " mol/kg ", solute2)
            value=1
            return {test_m1,test_m2}
    if value<0:
        print("Failed to converge")
        return {-1,-1}
#could use some multithreading to speed up fitting #TODO: fix this code to newer function calls
# def func(t2, t3):
#     m2=0.1+float(10.0/float(t2+1))
#     m3=1.5+float(10.0/float(t3+1))
#     b=virial_eq_electrolyte(1.678,m2,m3,0.044,0.037,0,-0.001)
#     if m2<0.2 and abs(b-2.503)<=0.0001:
#         print("value ", b, " produced by ", m2, " and ", m3)
#     return b
##########################################################################################
# Polynomial fits to CRC data for conversion between molality and molarity

def fit_polynomial_from_file(filename,solute_name,x_name,y_name, x_column,y_column,polynomial_degree, lines_to_remove_from_bottom, plotfit=False,save_figure=False):
    data_arr=np.genfromtxt(filename,delimiter=",",names=True,skip_footer=lines_to_remove_from_bottom)
    number_of_rows=data_arr.shape[0]
    print(number_of_rows)
    x_array=np.zeros(number_of_rows,dtype=float)
    y_array=np.zeros(number_of_rows,dtype=float)
    for i in range(number_of_rows):
        x_array[i]= data_arr[i][x_column]
        y_array[i]= data_arr[i][y_column]
    print("Imported Data from ",filename, " is: ",x_array, " for x and ", y_array, " for y." )
    fit_x=x_array.flatten()
    fit_y=y_array.flatten()
    poly_coefficients=np.polynomial.polynomial.polyfit(fit_x,fit_y,polynomial_degree)
    print("Terms A,B,C,D.. for the polynomial Ax^0+Bx^1+Cx^2+Dx^3... are: ", poly_coefficients)

    if plotfit==True:
        fig, ax = plt.subplots()
        title="Osmolality plotted with polynomial fit from " + x_name+ " VS "+y_name+ " from file: ."+ filename.strip(".csv")+"/"
        equation="$y="
        for term in range(len(poly_coefficients)):
            equation+= str(f"{round(poly_coefficients[term],7):.6f}")+"x^"+str(term)+" + "
        equation=equation.strip(" + ")
        equation=equation+"$"
        print("Plotting ", title)
        y_coord_fitted= 0
        x_coord_fitted= np.linspace(0, int(np.max(fit_x)+1), 200)
        y_coord_osmol=np.ones(len(x_coord_fitted))
        for terms in range(len(poly_coefficients)):
            y_coord_fitted+=poly_coefficients[terms]*x_coord_fitted**terms
        for i in range(len(y_coord_osmol)): 
            y_coord_osmol[i]=single_electrolyte_virial(virial_dictionary[solute_name]["kdiss"],y_coord_fitted[i],virial_dictionary[solute_name]["B"],virial_dictionary[solute_name]["C"])
        ax.plot(x_coord_fitted, y_coord_fitted,color='#FFC107',linestyle='dotted',label='fitted line')
        ax.plot(x_array,y_array,color='black',marker='x',linestyle='none',label='CRC data points')
        x_pos=int(x_coord_fitted.size/2)
        y_pos=int(y_coord_fitted.size/2)
        x_note=x_coord_fitted[x_pos]
        y_note=y_coord_fitted[y_pos]
        ax.plot(x_coord_fitted,y_coord_osmol,color='#D81B60',linestyle='dashed', label='osmolality from OVE')
        ax.set_xlabel(x_name)
        ax.set_ylabel(y_name)
        ax.set_title(title)
        ax.annotate(equation,xy=(x_note,y_note), xytext=((2*np.max(x_coord_fitted)/5), (np.min(y_coord_fitted))),arrowprops=dict(facecolor='black', shrink=0.05),annotation_clip=False)
        if save_figure==True:
            plt.savefig(filename.strip(".csv")+"_plot")
        #plot test point fitted from below
       # if solute_name=="EG":
       #     plt.plot(2.5,2.9318185249999997,'bo')
       # if solute_name=="GLY":
       #     plt.plot(3.5, 4.6401493355067505,'bo')
       # if solute_name=="NaCl":
       #     plt.plot(4.5,4.1135133170092635,'bo')
        secax=ax.secondary_yaxis('right',facecolor="#D81B60")
        secax.set_ylabel("Osmolality")
        secax.set_color('#D81B60')
        ax.spines['right'].set_color('#D81B60')
        ax.legend()
        plt.show()
    return poly_coefficients
""" converting between molarity and molality using polynomial fits found with the above function"""
def EG_molality_to_molarity(molality):
    molarity=(4.94707483e-02)+(9.25567439e-01)*molality+(-3.30583531e-02)*molality**2+(5.19836044e-04)*molality**3
    return molarity

def EG_molarity_to_molality(molarity):
    molality=(-0.1020293 )+(1.22474978)*molarity+(-0.03903846)*molarity**2+(0.01382168)*molarity**3
    return molality

def GLY_molarity_to_molality(molarity):
    molality=(2.88553718e-04 )+(1.00173256e+00)*molarity+(7.03911431e-02)*molarity**2+(6.33248557e-03)*molarity**3
    return molality
def GLY_molality_to_molarity(molality):
    molarity=(-2.47406091e-04)+(9.98060423e-01)*molality+(-6.98861715e-02)*molality**2+(3.93428154e-03)*molality**3
    return molarity

def NaCl_molarity_to_molality(molarity):
    molality=(-6.57836263e-04 )+(1.00094560e+00)*molarity+(-1.96631754e-02 )*molarity**2+(8.88368189e-05)*molarity**3
    return molality
   
def NaCl_molality_to_molarity(molality):
    molarity=(1.41031837e-04 )+(1.00118966e+00 )*molality+(1.78231574e-02 )*molality**2+(1.15175694e-03)*molality**3
    return molarity

 
#print(b)
def main():
    """ Finding molality in ternary solutions when only 1 osmolality is known using grid search"""

    """Values measured in experiment: """
    """
        0.5x PBS: .187 mmol/kg
        1x PBS: .323 osmol/kg
        2x PBS: .643 osmol/kg
        5x PBS: 1.600 osmol/kg
        15% glycerol: 2.185 osmol/kg
        15% EG: 2.503 osmol/kg
        15% glycerol+ 15% EG: 4.788 osmol/kg
        Holding media: .255 osmol/kg

    """
#    osmolality_05x_PBS=0.187
#    osmolality_1x_PBS=0.323
#    osmolality_2x_PBS=0.643
#    osmolality_5x_PBS=1.600
#    osmolality_GLY=2.185
#    osmolality_EG=2.503
#    osmolality_EG_and_GLY=4.788 
#    osmolality_HM=0.255
#    grid_search_ternary_1D(osmolality_HM,osmolality_EG,1,3,"NaCl","EG")
#    grid_search_ternary_1D(osmolality_HM,osmolality_GLY,1,3,"NaCl","GLY")
#   
#    grid_search_binary(virial_dictionary["NaCl"]["kdiss"],osmolality_05x_PBS,virial_dictionary["NaCl"]["B"],virial_dictionary["NaCl"]["C"], osmolality_05x_PBS/2, osmolality_05x_PBS*1.5)
#    grid_search_binary(virial_dictionary["NaCl"]["kdiss"],osmolality_1x_PBS,virial_dictionary["NaCl"]["B"],virial_dictionary["NaCl"]["C"], osmolality_1x_PBS/2, osmolality_1x_PBS*1.5)
#    grid_search_binary(virial_dictionary["NaCl"]["kdiss"],osmolality_2x_PBS,virial_dictionary["NaCl"]["B"],virial_dictionary["NaCl"]["C"], osmolality_2x_PBS/2, osmolality_2x_PBS*1.5)
#    grid_search_binary(virial_dictionary["NaCl"]["kdiss"],osmolality_5x_PBS,virial_dictionary["NaCl"]["B"],virial_dictionary["NaCl"]["C"], osmolality_5x_PBS/2, osmolality_5x_PBS*1.5)
#    grid_search_binary(virial_dictionary["NaCl"]["kdiss"],osmolality_HM,virial_dictionary["NaCl"]["B"],virial_dictionary["NaCl"]["C"], osmolality_HM/2, osmolality_HM*1.5)
#    binary
    """ if EG+GLY-HM has 2.058 mol/kg EG"""
    osm_eg=single_electrolyte_virial(1.0,2.057952897644882,virial_dictionary["EG"]["B"],virial_dictionary["EG"]["C"])
    print("Possible EG osmolality is: " , osm_eg)
    grid_search_ternary_1D(osm_eg,4.233,1,7,"EG","GLY")

    """ if EG+GLY-HM has 1.822 mol/kg EG"""
    osm_gly=single_electrolyte_virial(1.0,1.8227411370568527,virial_dictionary["GLY"]["B"],virial_dictionary["GLY"]["C"])
    print("Possible GLY osmolality is: " , osm_gly)
    grid_search_ternary_1D(osm_gly,4.233,1,7,"GLY","EG")
 
    """
    code gives values:

        Holding media (NaCl)  0.2550080729272943 osmole/kg has  0.15030348901744509  mol/kg of single solute
        EG  2.5029883007704425  Osmole/kg produced by  0.15030348901744509 HM  and  2.057952897644882 EG
        GLY  2.184964787088076  Osmole/kg produced by  0.15030348901744509 HM and  1.8227411370568527 GLY
        0.5x PBS at  0.1870087019952007 osmolality has  0.11054513522567613  mol/kg of single solute
        1x PBS at 0.32300834420843616 osmolality has  0.1898353166765834  mol/kg of single solute
        2x PBS at  0.6430078599913145 osmolality has  0.3729306121530608  mol/kg of single solute
        5x PBS at 1.6000013296562823 osmolality has  0.8944484722423612  mol/kg of single solute
        4.788-.255=4.233 EG+GLY, if EG is 2.058 mol/kg ->
        4.232901376334969  Osmole/kg produced by  2.057957808365267  mol/kg  EG , and  1.7420537102685514  mol/kg  GLY

        alternative solution:
        4.232904658767443  Osmole/kg produced by  1.8227478221442377  mol/kg  GLY , and  1.979924899624498  mol/kg  EG
            
           err: 1.8227411370568527 -1.7420537102685514= 0.08068742678
           err: 2.057957808365267 -1.979924899624498 = 0.07803290874

            both solutions are very close to predicted molalities but alternative solution is closer, so this one is used for simulation

        
        
    """
    """ NOTE: for the following code, plots are designed for molarity vs molality, 
    the fitted line will be correct but the osmolality line will be wrong for molality vs molarity"""

    """ Fitting Ethylene Glycol (EG) conversion with data from CRC"""
    fit_polynomial_from_file("EG_CRC_data_20_deg_C.csv","EG","molarity","molality",2,1,3,1,plotfit=False,save_figure=False)
#    fit_polynomial_from_file("EG_CRC_data_20_deg_C.csv","EG","molality","molarity",1,2,3,1,plotfit=False,save_figure=False)
#    
#    #testing the conversion function
#    m1_test=EG_molarity_to_molality(2.5)
#    print("m1 is ", m1_test)
    """ Fitting Glycerol (GLY) conversion with data from CRC"""
    fit_polynomial_from_file("./Glycerol_CRC_data_20_deg_C.csv","GLY","molarity","molality",2,1,3,1,plotfit=False)
#    fit_polynomial_from_file("Glycerol_CRC_data_20_deg_C.csv","GLY","molality","molarity",1,2,3,1,plotfit=False,save_figure=False)
#    
#    #testing the conversion function
#    m2_test=GLY_molarity_to_molality(3.5)
#    print("m2 is ", m2_test)
    """ Fitting NaCl conversion with data from CRC"""
    fit_polynomial_from_file("./NaCl_CRC_data_20_deg_C.csv","NaCl","molarity","molality",2,1,3,1,plotfit=False)
#    fit_polynomial_from_file("./NaCl_CRC_data_20_deg_C.csv","NaCl","molality","molarity",1,2,3,1,plotfit=False)
#     
#    #testing the conversion function
#    m3_test=NaCl_molarity_to_molality(4.5)
#    print("m3 is ", m3_test)
#   using values computed above:
    least_sig_fig_decimal=3
    nacl_molarity=NaCl_molality_to_molarity(0.15030348901744509)
    eg_molarity=EG_molality_to_molarity(2.057952897644882)
    gly_molarity=GLY_molality_to_molarity(1.8227411370568527)
    gly_mixed_molarity=GLY_molality_to_molarity(1.8227478221442377)
    eg_mixed_molarity=EG_molality_to_molarity(1.979924899624498 )
    PBS_05x_molarity=NaCl_molality_to_molarity(0.11054513522567613 )
    PBS_1x_molarity=NaCl_molality_to_molarity(0.1898353166765834 )
    PBS_2x_molarity=NaCl_molality_to_molarity(0.3729306121530608)
    PBS_5x_molarity=NaCl_molality_to_molarity(0.8944484722423612 )

    print("Initial molarities are: ", round(nacl_molarity,3), " mol/L NaCl, ", round(eg_molarity,3)," mol/L EG, and ", round(gly_molarity,3), " mol/L GLY for the single permeating CPAs. ")
    print("For EG & GLY the intial molarities are ", round(nacl_molarity,3), " mol/L NaCl, ",round(eg_mixed_molarity,3), " mol/L EG, and ", round(gly_mixed_molarity,3), " mol/L GLY." )
    print("For PBS solutions the intial molarities are ", round(PBS_05x_molarity,3), " mol/L for 0.5x, ",round(PBS_1x_molarity,3), " mol/L for 1x, ", round(PBS_2x_molarity,3), " mol/L for 2x, and ", round(PBS_5x_molarity,3), " mol/L for 5x PBS." )
if __name__ == '__main__':
    main()
