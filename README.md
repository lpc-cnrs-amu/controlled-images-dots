# Project COMPDOT
**Generating nonsymbolic number stimuli.**

(*from T. Gebuis &amp; B. Reynvoet, published in Behavioral Research Methods, 43, 2011*)

----

*The texts and files of the actual repository were hosted on a server that not functional anymore (as of July 2019; http://titiagebuis.eu). As the associated Matlab code was useful for generating images used as stimuli in cognitive psychology research studies, I made these scripts dowloadable again. If you are one of the original author, please consider managing this repository under your own account.*

---

**Abstract**

Studies investigating nonsymbolic numbers (e.g., dot arrays) are confronted with the problem that changes in numerosity are always accompanied by changes in the visual properties of the stimulus. It is therefore debated whether the visual properties of the stimulus rather than number can explain the results obtained in studies investigating nonsymbolic number processing. In this report, we present a program (available at http://titiagebuis.eu/Materials.html ; note that the program is designed to work with the Psychophysics Toolbox in MATLAB) that exports information about the visual properties of stimuli that co-vary with number (area extended, item size, total surface, density, and circumference). Consequently, insight into the relation between the visual properties of the stimulus and numerical distance can be achieved, and post hoc analyses can be conducted to directly reveal whether numerical distance or (some combinations of) the visual properties of a stimulus could be the most likely candidate underlying the results. Here, we report data that demonstrate the program's usefulness for research on nonsymbolic number stimuli.

----

**Original webpage information**

Script to create non-symbolic number stimuli Script is used to create the stimuli presented in: 
Gebuis, T. & Reynvoet, B. (2011). Generating non-symbolic number stimuli. *Behavior and Research Methods*.

 - Version January 18th 2012: comp_dots_version180112.m
To control for the visual cues the previous files created 4 images per number (larger number has a larger average dot diameter and a larger convex hull, (2) larger number has a larger average dot diameter but smaller convex hull, (3) larger number has a smaller average dot diameter and a smaller convex hull and (4) larger number has a smaller average dot diameter but larger convex hull. The current file randomly chooses one of the four options for each trial. Consequently only one image for each number is created. 
 - Version May 20th 2011: comp_dots_version200511.m
This version generates an additional output file called congr.txt. Here information is given about the congruency between the visual parameters and target number. For details see script.
 - Original version used in the manuscript: comp_dots_version090511.m

