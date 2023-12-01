#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "data_files.h"

struct DataArr{
    vector<float> mass_prim;
    vector<float> mass_secn;
    vector<float> proj_dist;
};

void obs_file::file_read(string argFName){

    DataArr data;
    fstream myfile;
    myfile.open(argFName, ios::in);

    string line;
    if(!myfile){
        cout << "!!!Error: No file found!!!" << endl;
        return;
    }

    int row_num = 0;
    while(getline(myfile, line, ',')){
        row_num++;
        if(row_num%8 == 4){
            data.mass_prim.push_back(std::stof(line));
        }
        else if(row_num%8 == 6){
            data.mass_secn.push_back(std::stof(line));
        }
        else if(row_num%8 == 0){
            data.proj_dist.push_back(std::stof(line));
        }
    } 

    vector<float> mass_ratio_arr;
    float tot_m = 0.;
    float mass_ratio = 0.;
    float min_ratio = 1.;

    int no_entry = data.mass_prim.size();
    cout << "Total JuMBOs: " << no_entry << endl;
    for(int i = 0; i < no_entry; i++){
        float m1 = data.mass_prim[i];
        float m2 = data.mass_secn[i];
        tot_m += m1 + m2;

        float mratio = min(m1/m2, m2/m1);
        mass_ratio_arr.push_back(mratio);
        mass_ratio += mratio;
        if(mratio < mass_ratio_arr[i-1]){
            min_ratio = mratio;
        }
    }

    float data_entry = data.mass_prim.size();
    float mean_mass = tot_m/data_entry;
    float mean_mass_ratio = mass_ratio/data_entry;

    float sum_of_squared_diff_m = 0;
    float sum_of_squared_diff_q = 0;
    for(int i = 0; i < no_entry; i++){
        float diff1 = mean_mass - data.mass_prim[i];
        float diff2 = mean_mass - data.mass_secn[i];
        sum_of_squared_diff_m += diff1 * diff1 + diff2 * diff2;

        float diffq = mean_mass_ratio - mass_ratio_arr[i];
        sum_of_squared_diff_q += diffq * diffq;
    }
    float standard_dev_m = sqrt(sum_of_squared_diff_m / (data_entry - 1));
    float standard_dev_q = sqrt(sum_of_squared_diff_q / (data_entry - 1));

    cout << "Mean Mass:                " << mean_mass << endl;
    cout << "Mass Standard Dev.:       " << standard_dev_m << endl;
    cout << "Mean Mass Ratio:          " << mean_mass_ratio << endl;
    cout << "Min. Mass Ratio:          " << min_ratio << endl;
    cout << "Mass Ratio Standard Dev.: " << standard_dev_q << endl;


}