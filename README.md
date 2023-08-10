# Find POTATOs

## Unresolved issues
- need requirements file
- need to update how to use 

[![Build Status](https://img.shields.io/static/v1.svg?label=CSL&message=software%20against%20climate%20change&color=green?style=flat&logo=github)](https://img.shields.io/static/v1.svg?label=CSL&message=software%20against%20climate%20change&color=green?style=flat&logo=github)


## Introduction
Congratulations! You are now the happy(?) user of findPOTATOs, an asteroid detection linking software. This is software created by N. Tan and C.R. Nugent. If you'd like more information on this code, feel free to read N. Tan's thesis at https://repository.wellesley.edu/object/ir1199 

### Installation
1. Clone this repository
2. Ensure Python3 is installed
3. Navigate to this folder using the terminal
4. Run python3 -m venv venv and source venv/bin/activate to set up and activate a virtual environment.
5. Use the `requirements.txt` file to install the requirements to run the script using `pip3 install -r requirements.txt`


## How to use FindPOTATOs.py
1. Within the directory containing findPOTATOs.py, there should be a sub-directory called "inputs". Inside this directory, there should be a plain text file called "linking_parameters.ini". Open it and edit it to toggle settings for findPOTATOs.py

2. Navigate to the directory containing findPOTATOs.py in a Terminal window, and run >> python findPOTATOs.py

3. If all goes well, the directory containing findPOTATOs.py should now have a new sub-directory called "outputs". In it you should find "tracklets.txt", containing all matched tracklets found. 


