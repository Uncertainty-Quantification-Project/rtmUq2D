# rtmUq2D

Reverse Time Migration for Uncertainty Quantification.

### MAIN FEATURES OF THE FORMULATION

- Physical and mathematical characteristics:

      -- Second-order wave equation

      -- Isotropic acoustic medium
      
      -- Invariant density

- Numerical characteristics:

      -- Finite Difference Method (Second-order in time and Eighth-order accuracy in space)
    
      -- Convolutional Perfectly Matched Layer for wavefield attenuation on the boundary

- Image Condition:

      -- Cross-correlation: correlation between source wavefield and receivers wavefield. 
      
                            This process produces amplitude squared.
____________________________________________________________


### COMPILATION AND RUNNING INSTRUCTIONS

- Compilation Requeriments:

      -- ICC compiler
      
      -- OpenMP Library
      
      -- MPI Library


- Compilation:

      -- Use make file typing:
      
         make clean
      
         make

- Local Running:

      -- Use the file "run_rtm_simulations.local":
         ./run_rtm_simulations.local

- Cluster Running:
   
      -- Configure the file "jobRtmSimulations.job" as desired
____________________________________________________________


### EXTRA INFORMATION

- For more details about input and output requirements see:

      ++ ./INPUTS/0_inputs_readme.INFO
      
      ++ ./OUTPUTS/0_outputs_readme.INFO
   
- To run this code, make download of the last version commited to
   a local repository.

- The program is setting to migrate one seismogram on two rtm simulations.

      ++ See ./INPUTS/0_inputs_readme.INFO for more details.
      
 ### CONTACT
 
 - E-mail: c.barbosa@nacad.ufrj.br
