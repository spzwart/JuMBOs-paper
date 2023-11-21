from binary_evolution import FinalInitialProperties
from read_data import ReadData


read_data = ReadData()
bin_plot = FinalInitialProperties()

time_crop = [False, True]

for crop_ in range(1):
    crop_ = True
    if (crop_):
        print("...Reducing for 9% Snapshot...")

    bin_plot.process_final_data(crop_)
    bin_plot.initialise(crop_)
    bin_plot.obs_scatter(0, crop_)
    STOP
    if not (crop_):
        
        bin_plot.obs_scatter(0, crop_)


        model_choices = [[0,1,2], [3,4,5], 
                         [1,7], [1,10], [6]]
        for models_ in model_choices:
            print("Time Evol for ", models_)
            bin_plot.time_evol_nJumbo(models_)

        print("Processing orbital cdf and mass cdf plots")
        model_choices = [[0,1,2], [0,3], [3,4,5], [0,6], 
                         [1,7,10], [1,7], [1,10],
                         [1,4,10,11]]
        for models_ in model_choices:
            bin_plot.orb_cdf_plot(models_, crop_)
            bin_plot.mass_CDF(models_, crop_)

        models = [0,1,2,3,4,5,6,7,8,9,10,11]
        for model_ in models:
            print("Processing pop stats and mass params for ", model_, bin_plot.models[model_])
            bin_plot.population_statistics(model_, crop_)
            bin_plot.mass_params(model_, crop_)

        print("Processing event statistics")
        bin_plot.event_statistics(crop_) 

        models = [8,9,10,11]#7,
        for model_ in models:
            bin_plot.two_point_correlation(model_, crop_)
        
    else:
        bin_plot.obs_scatter(0, crop_)
        print("Processing orbital cdf and mass cdf plots")
        model_choices = [[0,1,2], [0,3], [3,4,5], [0,6], 
                         [1,7,10], [1,7], [1,10]]
        for models_ in model_choices:
            bin_plot.orb_cdf_plot(models_, crop_)
            bin_plot.mass_CDF(models_, crop_)

        models = [0,1,2,3,4,5,6,7,8,9,10]
        for model_ in models:
            print("Processing pop stats and mass params for ", bin_plot.models[model_])
            bin_plot.population_statistics(model_, crop_)
            bin_plot.mass_params(model_, crop_)
            #bin_plot.two_point_correlation(model_, crop_)
            
