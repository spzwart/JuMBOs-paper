
#include <fstream>
#include <iostream>
#include <filesystem>

float CO2Emission(float &cpu_time){
    /*Calculate CO2 emission*/

    const float cpuWattsPerHour = 12.0;   // Power consumption of CPU in Watt/hr [Portegies Zwart 2020]
    const float cpuEmitPerkWatt = 0.283;   // Emission intensity in kWh/kg [Wittmann et al. 2013]
    float cons_wH = (cpu_time*cpuWattsPerHour)/60.;
    float total_ems = cons_wH/(cpuEmitPerkWatt*1000);
    return total_ems;
}

int main(){
    const std::string configs[12] = {
                                     "Fractal_rvir0.5", 
                                     "Fractal_rvir0.5_FF", 
                                     "Fractal_rvir0.5_FF_10Myr", 
                                     "Fractal_rvir0.5_FF_Obs", 
                                     "Fractal_rvir0.5_FFOnly", 
                                     "Fractal_rvir0.5_Obs",
                                     "Fractal_rvir0.5_Obs_Circ", 
                                     "Fractal_rvir1.0", 
                                     "Plummer_rvir0.5", 
                                     "Plummer_rvir0.5_FF",
                                     "Plummer_rvir0.5_FF_Obs", 
                                     "Plummer_rvir1.0"
                                    };

    const int cpu_line = 1;
    float total_cpu = 0.0;

    //string fdirec = "../data/"+str(model)+"simulation_stats"
    for (const std::string f:configs){
        std::string fdirec = "../data/Simulation_Data/"+f+"/simulation_stats/";
        for (const auto & entry : std::filesystem::directory_iterator(fdirec)){
            std::cout << entry.path() << std::endl;
            std::ifstream data_file(entry.path());
            if (!data_file.is_open()){
                std::cerr << "Error! File " << entry.path() << " not open." << std::endl;
            }

            std::string line;
            std::string str_cpu_time;
            int nrow = 0;
            while (getline(data_file, line)){
                nrow++;
                if (nrow == cpu_line){
                    for (char c:line){
                        if (isdigit(c) || c == '.'){
                            str_cpu_time += c;
                        }
                    }
                    float cpu_time = stof(str_cpu_time);
                    total_cpu += cpu_time;
                }
            }
        }
    }
    std::cout << "Total cpu time in minutes: " << total_cpu << std::endl;
    float total_ems = CO2Emission(total_cpu);
    std::cout << "Total emission in kg: " << total_ems;
}