# analytical-gradient-validation-for-shape-optimization

This repository showcases how to use MATLAB's API `checkGradients` to validate analytically derived expressions of gradients with application to sensitivity analysis for shape optimization problems in structural mechanics. 

Author: Deha Şen Köse, dehasenkose@gmail.com

## Prerequisites
The coding language used is MATLAB, and MATLAB2024a release is used throughout the development. Thus, it is advised that either MATLAB2024a or a newer version should be used.
The following MATLAB toolboxes or newer releases are required to run the available code.
```
Optimization Toolbox                                  Version 24.1        (R2024a)
```
## How to run?
1. First, please extract the repository from the zipped folder. Then open Matlab and navigate to the same folder using the directory tab that lies from left to right, close to the top of the page.
2. Please run the project file called **"analytical-gradient-validation-for-shape-optimization"** to add the necessary files to the directory.
3. Go to the **applications** folder and then to the **"cantileverBeam"**. Here, you will find the main file called **"main_CAD_Shape_Optimization.mlx"**. Click on the run to run the code.
4. The gradient validation using the checkGradients function can be seen with the following output:
   ```
    ____________________________________________________________
    
    Objective function derivatives:
    Maximum relative difference between supplied 
    and finite-difference derivatives = 1.08938e-08.
    
    checkGradients successfully passed.
    ____________________________________________________________
   ```  
5. The bottom part of the code shows two examples of different optimization methods applied to the plane stress IGA plate in membrane action. Just so you know, there are different constraints for these problems.
6. Please refer to the explanations of the well-documented functions in case of an unclear point.

## Troubleshooting
If you encounter any bug that prevents you from running the live scripts, please use the correct MATLAB version with the toolboxes mentioned above. You can go ahead and follow Matlab's output to solve the issue.
If you have any questions, please feel free to contact me using the email above.
