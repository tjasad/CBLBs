import matplotlib.pyplot as plt 
from solver import * 
import numpy as np  
import os.path      
import pickle     
import seaborn as sns
import pandas as pd
import multiprocessing
from random import randint

from model_clb import model_clb, param_references
 


def getParamDistrib(number_points = 0, file_name = ""):    
    rand_samples = []    
    
    region = model_regions[0]   
    if number_points:
        samples = region.points[np.random.choice(region.points.shape[0], number_points, replace=False), :]
    else:
        samples = region.points

    rand_samples.append(samples)

    rand_samples = np.array(rand_samples)

    param_names = [r"$\delta_L$", 
                    r"$\gamma_L$", 
                    r"$n_y$", 
                    r"$\theta_L$", 
                    r"$\eta_x$", 
                    r"$\omega_x$", 
                    r"$m_x$", 
                    r"$\delta_x$",
                    #r"$\delta_y$", #to remove
                    r"$\gamma_x$",
                    r"$\theta_x$",
                    r"$\rho_x$"]
    
    df = pd.DataFrame(rand_samples[0])
    df.columns = param_names

    if file_name:
        df.to_csv(file_name, index=False)
    return df

def plotParamsdf(df=None, number_points = 0):
    if not type(df):
        df = getParamDistrib(number_points)
    
    param_names = [ r"$\gamma_L$", 
                    r"$\eta_x$", 
                    r"$\gamma_x$",
                    r"$\theta_L$", 
                    r"$\omega_x$", 
                    r"$\theta_x$",
                    r"$\delta_L$",  
                    r"$\delta_x$",
                    #r"$\delta_y$", #to remove                  
                    r"$\rho_x$",
                    r"$n_y$",                   
                    r"$m_x$"]

    fig, axes = plt.subplots(4,3)

    for i, param_name in enumerate(param_names):
        if param_name:
            ax = axes.flat[i]
            sns.violinplot(data=df[param_name], ax = ax) #,palette="Pastel1")
            ax.set_xticks([0])
            ax.set_xticklabels([param_name])        

    """
    for param_id in range(len(param_names)):
        ax = axes.flat[param_id]

        sns.violinplot(y = param_names[param_id], x="Model id", data=df[[param_names[param_id], "Model id"]], ax = ax) #,palette="Pastel1")
    """ 
    fig=plt.gcf()

    fig.set_size_inches([20,12])
    plt.savefig(os.path.join(base_path_robustness, 'params_distrib_sns.pdf'), bbox_inches = 'tight')
    plt.show()

def test_random_point()):
    points = model_regions[0].points
    candidate = tuple(points[randint(0,len(points)-1)])
    
    #out = model.simulate(candidate[:8]+candidate[9:], plot_on=True)
    out = model.simulate(candidate, plot_on=True)


if __name__ == "__main__":

    sns.set_style("white")
    #flatui = ['#d9d9d9','#bdbdbd','#969696','#636363']
    #sns.set_palette(flatui)

    #
    # SETTINGS
    #
    read_data = True
    ga_solutions = True
    local_solutions = False

    
    
    base_paths_opt = ["results"]
    base_path_robustness = "results"


    #
    # END OF SETTINGS
    #

    if read_data:


        #folders = [os.path.join(base_path, "one_bit_model"), os.path.join(base_path, "two_bit_model"), os.path.join(base_path, "three_bit_model")]   
        model = model_clb()
        solver = Solver(model)

        model_regions = []

        model_str = ""
        region_files =  []
        for base_path_opt in base_paths_opt:
            if ga_solutions:
                region_files.append(os.path.join(base_path_opt, model_str+"ViableSet_IterGA.p"))
            if local_solutions:
                for i in range(10):
                    region_files.append(os.path.join(base_path_opt, model_str+"Region0ViableSet_Iter" + str(i+1) + ".p"))

        viablePoints = []   
        for region_file in region_files: 
            viablePointsRegion = pickle.load(open(region_file, "rb"))   
            print(len(viablePointsRegion))   
            viablePoints.extend(viablePointsRegion)
        print("Number of points: ",len(viablePoints))
        region = Region(viablePoints, model, "region")              
        model_regions.append(region)                                                                        
        
    


 
    ##df = getParamDistrib(file_name="results_robustness\\params.csv")
    #df = pd.read_csv("results_robustness\\params.csv")
    #plotParamsdf(df)

    test_random_point()
