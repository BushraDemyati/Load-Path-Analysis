## Load Path Analysis Using the U* Index Method

_Created with MATLAB R2024b. Compatible with R2024b and later releases._ 

## Description

The codes provided utilize the [Partial Differential Equations (PDE) Toolbox](https://de.mathworks.com/products/pde.html) in MATLAB to solve a static structural problem, namely, plates in membrane action.

Using the information derived from the solver, the U* index field is calculated and the load paths are derived.

## Owner

Bushra Demyati, M.Sc. ([bushra.demyati@tum.de](mailto:bushra.demyati@tum.de))

## Contents
The repository contains the following Live Scripts which solve and find the load paths for different setups:
- **``Geometry.stl``**

This is the Standard Triangle Language (STL) file that contains the necessary information about the geometry.

- **``PinnedVertices_CompressiveEdgeLoad.mlx``**

This is the file that contains the load path analysis for the plate in membrane action shown below:

![image](https://github.com/user-attachments/assets/056f2138-c156-45d3-83e1-07681e222be2)

- **``PinnedVertices_PointLoadMiddleofEdge.mlx``**

This is the file that contains the load path analysis for the plate in membrane action shown below:

![image](https://github.com/user-attachments/assets/d328916c-f8be-48bb-8d74-4f1c890cd0c8)

- **``PinnedVertices_PointLoadYDirection.mlx``**

This is the file that contains the load path analysis for the plate in membrane action shown below:

![image](https://github.com/user-attachments/assets/6f7545c2-7538-47c8-a086-d37eb7ab8f90)

- **``ClampedEdge_CompressiveEdgeLoad.mlx``**

This is the file that contains the load path analysis for the plate in membrane action shown below:

![image](https://github.com/user-attachments/assets/8b45c89c-d675-4498-b7da-214c55995f28)

- **``ClampedEdge_PointLoadYDirection.mlx``**

This is the file that contains the load path analysis for the plate in membrane action shown below:

![image](https://github.com/user-attachments/assets/0b72dd0a-7833-40dd-97ef-3008bef86816)

- **``ClampedEdge_EdgeLoadYDirection.mlx``**

This is the file that contains the load path analysis for the plate in membrane action shown below:

![image](https://github.com/user-attachments/assets/7f90fb79-0d39-450d-902c-4bd488aff21c)

- **``SpringSupports.mlx``**

This is the file that contains the load path analysis for a plate in membrane action that is supported by springs that are only activated when under compression. The process to deactivate the springs under tension is also performed using the PDE Toolbox in MATLAB. The schematic sketch of the analyzed plate is shown below:

![image](https://github.com/user-attachments/assets/6edc7e38-27ff-4478-9a12-62f7ad8b0526)

## Results
The codes solve the static problems using the PDE Toolbox, allowing for the stresses, strains, and displacements to be visualized using the ['Visualize PDE Results'](https://de.mathworks.com/help/pde/ug/visualizepderesults.html?searchHighlight=visualize+pde+results&s_tid=srchtitle_support_results_1_visualize+pde+results) Live Task in MATLAB.
The U* index field is calculated based on the solutions of the system. The ['streamline'](https://de.mathworks.com/help/matlab/ref/streamline.html) function in MATLAB is used to find the steepest gradients whcih form the load paths. An example of the U* index field and corresponding load paths is shown below:

![image](https://github.com/user-attachments/assets/6181bc07-198c-4222-854d-94838501d03e)

Moreover, the variation of the U* index values along the load path (the uniformity condition) is also calculated in the codes. An example of that is shown below:

![image](https://github.com/user-attachments/assets/9ecd5d06-323c-4661-9daa-0fd1752193ac)

To get more information about how the gradients of the U* index field vary along the x- and y-directions, the streamlines in the entire span of the plates are also plotted using the ['streamslice'](https://de.mathworks.com/help/matlab/ref/streamslice.html?searchHighlight=streamslice&s_tid=srchtitle_support_results_1_streamslice) function in MATLAB. An example is shown below:

![image](https://github.com/user-attachments/assets/4838cd63-2b4b-4aa8-a858-cda97c76300f)


## Concepts
Finite Element Methods, structural analysis, plates in membrane action, U* index method, load path analysis, [Live Tasks](https://de.mathworks.com/help/matlab/matlab_prog/add-live-editor-tasks-to-a-live-script.html), [PDE Toolbox](https://de.mathworks.com/products/pde.html), plotting, [Streamlines](https://de.mathworks.com/help/matlab/ref/streamline.html).

## Suggested Audience
All engineering disciplines, such as, civil engineers, mechanical engineers, etc.

## Workflow
Open the Live Script of the desired structural probelm and hit Run! 

Please note that the values of the load, or even the boundary conditions, could be changed. The inly thing that should be taken into consideration is that the load and support nodes are updated accordingly. The rest of the code should function without any other changes.

Moreover, depending on the loading and, in case these codes are used as reference for more complex problems, it is important to increase the number of interpolation points as necessary after the U* index is calculated for the triangular mesh created using the PDE Toolbox. 

There is also the possibility to assess the uniformity condition of different load paths. Simply, change the number of the retrieved line from the 'lineobj' variable. 

## Release last tested
R2024b

## References

K. Marhadi and S. Venkataraman, “Comparison of Quantitative and Qualitative Information Provided by Different Structural Load Path Definitions,” International Journal for Simulation and Multidisciplinary Design Optimization, vol. 3, pp. 384–400, Jul. 2009, doi: 10.1051/ijsmdo/2009014.

H. Hoshino, T. Sakurai, and K. Takahashi, “Vibration reduction in the cabins of heavy-duty trucks using the theory of load transfer paths,” JSAE Review, vol. 24, no. 2, pp. 165–171, Apr. 2003, doi: 10.1016/S0389-4304(03)00005-5.
