# Thesis_ERiQ_Model
Computational Matlab model simulating stress response pathways in an energy restricted quiescent cell experimental platform. All work was done for my PhD dissertation and was published December, 2017.

Alfego, D., & Kriete, A. (2017). Simulation of Cellular Energy Restriction in Quiescence (ERiQ)—A Theoretical Model for Aging. Biology, 6(4), 44.

ERiQ Computational Model Instructions
**All .m files need to be in the same folder directory in Matlab

•	ERiQ.m
-This file is not run, it is only called by the ode solver
-This calls separate functions for each equation. These are each their own .m files labeled ‘f_.m’
    -Separated for organization and ease of optimization

•	ERiQ_event.m
  -Does not run, but needs to be called in the ODE solver as options to terminate the system if MFUNCT reaches 0.5.

•	ERiQ_ODE.m
  -Full aging simulation
  -Use this file to run the odesolver on ERiQ.m
    -It uses initial conditions established via homeostasis (our original conditions should be in a comment)
    -The code will establish variables for you – i.e. MFUNCT, ATP, AKT, etc.
    -The code will auto-plot these major nodes
  -Global constants will be colored in light blue – these affect ALL functions in the folder
    -P53_Act and MDR can be changed to affect changes
    -To Run the Local Sensitivity Analysis
      -Each ‘f_.m’ function is multiplied by a global “SA” value that is set at 1 for normal conditions within the ERiQ_ODE.m code
        -Labeled ‘_SA’ (i.e. NFKB_SA = 1)
      -Changing this SA value will manipulate the function the name corresponds to.
        -i.e. NFKB_SA = 1.1  a 10% increase in NFKB activity
        -i.e. NFKB_SA = 0.9  a 10% decrease in NFKB activity
      -To change ROS levels, please use ROS_SA2 and NOT ROS_SA
      -Rerun ERiQ_ODE.m for new values
        -‘Parameters’ table will list new values for you
        -Record values in Excel for ease of SOF calculation
        -Don’t forget to reset the ‘_SA’ global value back to 1 before looking at the next sensitivity analysis.

•	ERiQ_Pulse
  -Simulation of ERiQ pulses in MDAMAGE and MFUNCT (not full aging simulation)
  -Will plot two figures representing the responses to these pulses, simply run the whole script
  -Should you wish to change strength of the pulse
    -Rerun individual assessment of the simulation between t = [0 3000] and t = [3000 6000] where the new pulse occurs in the second time range; plot these on the same figure
      -Use conditions found in Y(end) of the first time range as base for the initial conditions used in the odesolver in the second time range
      -The pulse is reflected by reducing or increasing the value of a specific initial condition in the second range
        -I.e. at end of t = [0 3000], MDAMAGE = 0.0724. To reflect a pulse increase of 0.5, the initial condition for MDAMAGE during t = [3000 6000] would be 0.0724+0.5 = 0.5724.
        -You can change these as you see fit, but you need to update the initial conditions to reflect new pulse

•	ERiQ_homeostasis.m
  -Generates plot of varying MDR values and effect on lifespan/MFUNCT
  -Can run as is; change MDR values (for Ys 1 through 7) for different simulations

•	ERiQ_simulation.m
  -Compilation of all .m files in one script that can be published as PDF or html
  -Does not run – will yield error. For supplemental file production only

•	ERiQ_SA.m
  -Used to run multi-parameter local sensitivity analyses. This file is NOT run, similar to how ERiQ.m is not run – it just establishes the equations
  -P53 vs MDR
    -Run script ERiQ_SA_p53_MDR to use loop and generate 3d plot
      -Can edit the range of values used here for parameter analysis – but if changed, make sure you also update the plot axis ticks accordingly
  -NFKB vs MDR
    -Run script ERiQ_SA_NFKB_MDR to use loop and generate 3d plot
      -Can edit the range of values used here for parameter analysis – but if changed, make sure you also update the plot axis ticks accordingly

•	ERiQ_Global_Analysis_RUN.m
  -Script to run a global sensitivity analysis
  -Script calls the function ERiQ_Global_Analysis.m. Do not run that function, only use the file labeled RUN.
  -Run this script to perform a global Monte Carlo sensitivity analysis that perturbs each parameter in the same simulation
    -N = 1000, but can be changed in code
    -Gaussian distribution, mean = 1, std = 0.033
      -Represents max 10% increase or decrease
      -For higher increase, change standard deviation
  -Saves 2 tables
    -GSA = actual random SA parameters chosen
    -GSA_percent = % change in parameter change and final values
    -Can copy these tables into Excel for further analysis

•	ERiQ_treatments.m
  -Generates plots for lifetime inhibition of NFKB, AKT, mTOR and Autophagy
  -Follow comments to change these inhibition concentrations
    -Follows the idea of sensitivity parameter changes

•	hline.m and vline.m
  -Funtions from MathWorks to place vertical line at specific point on x-axis. Needed for pulse plots
  -Simply need to be in the same folder directory, do not change
