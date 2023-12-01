#include <iostream>
#include <string>
#include <cmath>

#include "src/data_files.h"


using namespace std;

int main(){
    obs_file obs_obj;
    string file_name = "data/observations/src/obs_data.txt";

    obs_obj.file_read(file_name);
}