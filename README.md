# Description
This is a python script to analyse nanopore data generated from the use of polymer electrolyte. Due to the high signal to noise ratio generated from the polymer electrolyte, some of the steps become trivial, for example a detail calculation on the thresholds can be omitted as translocation events are very visible. At the moment, the translocation events shape does not encode further information, thus detail tracing on the event shape is redundant. It utilises not complicated method to calculate the basic characteristics of the translocation events.


For details on polymer electrolyte please read below, the data from these paper can be analysed through this python script:

Macromolecular crowding enhances the detection of DNA and proteins by a solid-state nanopore

[https://pubs.acs.org/doi/10.1021/acs.nanolett.0c02246]

Probing RNA Conformations Using a Polymer–Electrolyte Solid-State Nanopore

[https://pubs.acs.org/doi/10.1021/acsnano.2c08312]

Mechanistic Study of the Conductance and Enhanced Single-Molecule Detection in a Polymer–Electrolyte Nanopore

[https://doi.org/10.1021/acsnanoscienceau.2c00050]

Nanopore fingerprinting of supramolecular DNA nanostructures

[https://doi.org/10.1016/j.bpj.2022.08.020]

Next-Generation Nanopore Sensors Based on Conductive Pulse Sensing for Enhanced Detection of Nanoparticles

[https://doi.org/10.1002/smll.202305186]


# Packages
  - tkinter
  - pyabf [https://pypi.org/project/pyabf/]
  - numpy
  - pandas
  - matplotlib
  - pybaselines [https://pypi.org/project/pybaselines/]
  - seaborn
  - math
  - scipy
  - palettable [https://pypi.org/project/palettable/]

# File Input
By default it is to use with .abf file (Molecular Devices). However the downstream processing does not rely on .abf file as the pyabf packages will read and then convert into arrays of time and current datapoints, thus, with other means of data conversion, one can load two arrays of data points then use this script for analysis.

# Files Output
  - Excel file for all the detected events, with the following column outputs:
    - Event count (from 1)
    - Event positions (raw)
    - Event positions (in seconds after adjustment with sampling frequency)
    - Event Start position (for translocation peak plotting)
    - Event End position (for translcoation peak plotting)
    - Event height (in pA and nA)
    - Event width (in ms)
    - Event area (pA⋅ms and nA⋅ms, in terms of Area under the curve, through summation of all the events height across a defined time within the boundary defined by the event Start and event End (i.e. integration))


# Examples
The script will output the following figure for the users and save it the .abf file location in the computer.

Figure 1: Baseline Statistics-sensitive Non-linear Iterative Peak-clipping (SNIP) fitting 
![Example trace_01_baseline_fitting](https://github.com/chalmers4c/Nanopore_event_detection/assets/97460687/bd2bcbcb-9c73-4b9f-9134-b8e2e50f7307)

Figure 2: Baseline adjustment to move the trace through the SNIP fit
![Example trace_02_adjusted_baseline](https://github.com/chalmers4c/Nanopore_event_detection/assets/97460687/51abad1d-e31d-4b37-bc18-8014fa083705)

Figure 3: 1x2 plot to show the raw trace and the baseline adjusted plot
![Example trace_03_baseline_and_adjusted](https://github.com/chalmers4c/Nanopore_event_detection/assets/97460687/927e3c31-3421-401a-8f47-a41844b1866e)

Figure 4: Translocation events calling map on top of the adjusted baseline
![Example trace_04_events_calling](https://github.com/chalmers4c/Nanopore_event_detection/assets/97460687/4c3f27b8-759f-451f-915d-9038a24c8194)

Figure 5: 20 random translocation events
![Example trace_05_20_random_events](https://github.com/chalmers4c/Nanopore_event_detection/assets/97460687/964b86df-57b7-427a-858f-8eeac66ff5a8)

Figure 6: All detectable translocation events overlaid
![Example trace_06_all_events_overlay](https://github.com/chalmers4c/Nanopore_event_detection/assets/97460687/c509f766-6803-4d33-87ec-987b9a0c2a2e)

Figure 7: 50 random 20 translocation events overlay, this is a screening map to find a seed number that best represent your data
![Example trace_07_seed_screening_for_overlay_plot](https://github.com/chalmers4c/Nanopore_event_detection/assets/97460687/380c0e7a-346c-491b-b3f7-77f5acc80e46)

Figure 8: 20 translocation events overlay, this is to showcase translocation events that best represents your dataset
![Example trace_08_20_events_overlay](https://github.com/chalmers4c/Nanopore_event_detection/assets/97460687/135c4651-1877-4278-98fb-78a79be7ebe2)

# Acknowledgement
The script is largely based on spikesandburst blog post after heavily editing. For the original tutorial: [https://spikesandbursts.wordpress.com]
