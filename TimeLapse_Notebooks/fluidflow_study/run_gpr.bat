@echo off 
cd /d "C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\fluidflow_study" 
conda run --no-capture-output -n gprMax python -m gprMax "shift_0p5lambda\resolution_0p5lambda.in" -n 380 -gpu 0 
