# gradient-validation-msc-thesis-optimization
This code is an integral part of the M.Sc. The thesis is titled "Sensitivity Analysis and Shape Optimization in Structural Mechanics."

Below, one can find the validation of the gradients used in the optimization procedure and two examples of optimization using two different methods.
Please note that the thesis and the code are the intellectual property of the original writer. Commercial use in any way has to be permitted by the Faculty of Design and Engineering of the Technical University of Munich and the author himself.
Please refer to the corresponding function scripts and/or the thesis for the resources used in the code.

Author: Deha Şen Köse, deha.koese@tum.de

## Prerequisites
The coding language used is Matlab, and Matlab2024a release is used throughout the development. Thus, it is hereby advised that either Matlab2024a or a newer version should be used.
The following Matlab toolboxes or newer releases are required to run the available codes.
```
Global Optimization Toolbox                           Version 24.1        (R2024a)
Optimization Toolbox                                  Version 24.1        (R2024a)
Symbolic Math Toolbox                                 Version 24.1        (R2024a)
```
## How to run?
1. First, please extract the repository from the zipped folder. Then open Matlab and navigate to the same folder using the directory tab that lies from left to right, close to the top of the page.
2. Please run the project file called **"Gradientvalidationmscthesisoptimization"** to add the necessary files to the directory.
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
5. Two examples of different optimization methods applied to the plane stress IGA plate can be seen in the bottom part of the code. Please note that there are different constraints for these problems.
6. Please refer to the explanations of the well-documented functions and the thesis in case of an unclear point.

## Troubleshooting
If you encounter any bug that prevents you from running the live scripts, please use the correct Matlab version with the toolboxes mentioned above. You can go ahead and follow Matlab's output to solve the issue.
If you have any questions, please do not hesitate to contact me using the email above.
