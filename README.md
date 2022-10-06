# BREWPITOPES
Suite of tools to manage epitope prediction results from linear and structural origin and to integrate a pipeline of prioritization filters to curate B-cell epitopes specific for neutralizing antibody recognition.

## INSTALLATION (DOCKER IMAGE)
0. Download the Dockerfile from this [repository](https://github.com/AlbertCS/Brewpitopes_container). To compile the Dockerfile you will need to have docker installed. 

2. Create docker image from Dockerfile (may take a while):
```
sudo docker build -t brewpitopes PATH/TO/Dockerfile
``` 

2. Create a shared folder between Brewpitopes docker image and your local machine.

```
sudo docker run -it --volume /host/your/path/to/brewpitopes_projects:/home/Projects brewpitopes 
```

Change /host/your/path/to/brewpitopes_projects for your desired local directory.

Once you add files/folders in this directory, they will appear automatically within the brewpitopes docker image at /home/Projects.

3. Explore the folder structure of the Docker image named "brewpitopes".
```
ls
```
You should see the folders: brewpitopes, example, Projects

3.1 To execute the scripts of the pipeline go to:
```      
cd brewpitopes
```

3.2 To find an example of a brewpitopes project and the required files go to
```
cd ../example
```

3.3 To create new project, do so within the Projects folder:
```
cd ../Projects
mkdir your_project
```

3.4 Once you have created "your_project" folder, to start running the pipeline move back to the folder "brewpitopes" and you are ready to go !!

## PIPELINE
** All the steps shoud be run within the Docker image except when it indicates "(Locally)".

4. Use directories.R to create the folder environment.

```
Rscript directories.R --path ../Projects/your_project
```

5. (Locally) Download the FASTA file of the target protein at [Uniprot](https://www.uniprot.org/).    
      
      Save at /host/your/path/to/brewpitopes_projects/your_project/brewpitopes/Z_fasta  
      
6. (Locally) Use the FASTA to predict linear epitopes using [Bepipred 2.0] (https://services.healthtech.dtu.dk/service.php?BepiPred-2.0) server and export results as csv (default parameters).  
       
      Save at /host/your/path/to/brewpitopes_projects/your_project/brewpitopes/A_linear_predictions/bepipred/bepipred_results.csv  
      
7. Extract epitopes from Bepipred results using epixtractor_linear_bebipred.py.  
```
python3 epixtractor_linear_bepipred.py
```
```
Add path to bepipred results: ../Projects/your_project/brewpitopes/A_linear_predictions/bepipred/bepipred_results.csv
Add path to output folder: ../Projects/your_project/brewpitopes/C_epixtractor    
```
      
8. (Locally) Use the FASTA to predict linear epitopes using [ABCpred](https://webs.iiitd.edu.in/raghava/abcpred/ABC_submission.html) server.

      Predict using all the epitope windows (10,12,14,16,18,20) and overlapping filter ON.  
      Copy results from the webpage table to a .csv  
      Save as: /host/your/path/to/brewpitopes_projects/your_project/brewpitopes/A_linear_predictions/abcpred/abcpred_10mers.csv 
      
9. Extract epitopes from ABCpred results using epixtractor_linear_abcpred.R  

````
Rscript epixtractor_linear_abcpred.R --outpath your/path/to/brewpitopes/C_epixtractor --input_10mers ../Projects/your_project/brewpitopes/A_linear_predictions/abcpred/abcpred_10mers.csv --input_12mers ../Projects/your_project/brewpitopes/A_linear_predictions/abcpred/abcpred_12mers.csv --input_14mers ../Projects/your_project/brewpitopes/A_linear_predictions/abcpred/abcpred_14mers.csv --input_16mers ../Projects/your_project/brewpitopes/A_linear_predictions/abcpred/abcpred_16mers.csv --input_18mers ../Projects/your_project/brewpitopes/A_linear_predictions/abcpred/abcpred_18mers.csv --input_20mers ../Projects/your_project/brewpitopes/A_linear_predictions/abcpred/abcpred_20mers.csv
````
      
10. (Locally) Download the PDB file of the target protein at PDB DB. 
       Save at /host/your/path/to/brewpitopes_projects/your_project/brewpitopes/B_structural_predictions/pdb  
       
11. (Locally) Use [PDBrenum](http://dunbrack3.fccc.edu/PDBrenum/) server to renumerate the PDB residues according to its corresponding FASTA file in Uniprot.  
      Download results as .pdb  
      Save at /host/your/path/to/brewpitopes_projects/your_project/brewpitopes/B_structural_predictions/pdbrenum 
      
12. (Locally) Use the renumbered PDB to predict structural epitopes using [Discotope 2.0](https://services.healthtech.dtu.dk/service.php?DiscoTope-2.0) server and export the results as csv.    
      Default threshold.  
      Select chain A by default.
      Save at ../Projects/your_project/brewpitopes/B_structural_predictions/discotope  
      
13. Extract epitopes from Discotope results using epixtract_structural.py
```
python3 epixtract_structural.py  
```
```
Add path to discotope results: ../Projects/your_project/brewpitopes/B_structural_predictions/discotope/discotope_results.csv  
Add path to output folder: ../Projects/your_project/brewpitopes/C_epixtractor
```

14. Merge the epitopes extracted from Bepipred, ABCpred and Discotope results using epimerger.R
```
Rscript epimerger.R --abcpred ../Projects/your_project/brewpitopes/C_epixtractor/abcpred_results_extracted.csv --bepipred ../Projects/your_project/brewpitopes/C_epixtractor/abcpred_results_extracted.csv --discotope ../Projects/your_project/brewpitopes/C_epixtractor/discotope_results_extracted.csv --outdir ../Projects/your_project/brewpitoeps/D_epimerger
```

Take steps 15, 16 and 17.1 if you want to predict protein topology using CCTOP. Otherwise, if the topology of your protein is already described and you want to add it manually go directly to step 17.2.

15. (Locally) Predict the protein topology using [CCTOP](http://cctop.enzim.ttk.mta.hu/?_=/jobs/submit) server.  
      Donwload results as .xml.
      Save at /host/your/path/to/brewpitopes_projects/your_project/brewpitopes/E_topology/CCTOP/cctop.xml
      
16. Extract the topological domains using xml_cctop_parser.R  (ONLY IF 
```
Rscript xml_cctop_parser.R --xml ../Projects/your_project/brewpitopes/E_epitopology/CCTOP/cctop.xml --outdir ../Projects/your_project/brewpitopes/E_epitopology/CCTOP
```

17. Label the epitopes based on their topology (intracellular, membrane or extracellular) using epitopology.R
17.1 Using CCTOP predictions --> use epitopology_cctop.R
```
Rscript epitopology_cctop.R --input_CCTOP ../Projects/your_project/brewpitopes/E_epitopology/CCTOP/cctop_domains.csv --input_epitopes ../Projects/your_project/brewpitopes/D_epimerger/merged.csv --outdir ../Projects/your_project/brewpitopes/E_epitopology
```

17.2 Manual topology annotation --> use epitopology_manual.R. Add manually the starting and ending positions of your extracellular domains at --start_pos and end_pos.
```
Rscript epitopology_manual.R --start_pos 1,12,22 --end_pos 8,18,28 --input_epitopes ../Projects/your_project/brewpitopes/D_epimerger/merged.csv --outdir ../Projects/your_project/brewpitopes/E_epitopology
```
      
18. (Locally) Predict the glycosilation profile of the protein using the FASTA file.  
      N-GLYCOSILATIONS at [NetNGlyc 1.0](https://services.healthtech.dtu.dk/service.php?NetNGlyc-1.0) server.    
      COPY MANUALLY THE DATAFRAME HEADED: SeqName	Position	Potential	Jury_agreement	NGlyc_result	Prediction  
      SAVE AS CSV at /host/your/path/to/brewpitopes_projects/your_project/brewpitopes/F_epiglycan/netnglyc  
      
      O-GLYCOSILATIONS AT [NetOGlyc 4.0](https://services.healthtech.dtu.dk/service.php?NetOGlyc-4.0) server.
      COPY MANUALLY THE DATAFRAME HEADED: seqName  	source	feature	start 	end	score strand      frame       comment  
      SAVE AS CSV at /host/your/path/to/brewpitopes_projects/your_project/brewpitopes/F_epiglycan/netoglyc
      
19. Extract the glycosilated positions from both N-glyc and O-glyc outputs using epiglycan_extractor.R
```
Rscript epiglycan_extractor.R --oglyc ../Projects/your_project/brewpitopes/F_epiglycan/netoglyc/oglyc.csv --nglyc ../Projects/your_project/brewpitopes/F_epiglycan/netnglyc/nglyc.csv --outdir ../Projects/your_project/brewpitopes/F_epiglycan/
```

20. Use epiglycan.py to label the glycosilated epitopes.  
```
python3 epiglycan.py
```
```
Add path to input epitopes: ../Projects/your_project/brewpitopes/E_epitopology/topology_extracted.csv
Add path to output folder: ../Projects/your_project/brewpitopes/F_epiglycan
Add path to extracted glycosilated positions: ../Projects/your_project/brewpitopes/F_epiglycan/glycosilated_positions.csv  
```

21. (Locally) Use ICM_browser (MOLSOFT) to extract the RSA values for accessibility calculation.  
      Download ICM_browser from [http://www.molsoft.com/icm_browser.html](http://www.molsoft.com/icm_browser.html)
      Move the PDB renumbered file of the corresponding protein (step 11) to a local folder.
      Move the script Compute_ASA.icm to the same local folder.
      Execute in the command line of the programme ICM Browser the code in Compute_ASA.icm  
      Save results at /host/your/path/to/brewpitopes_projects/your_project/brewpitopes/G_episurf
      
22. Extract the buried positions using icm_extractor.R  
```
Rscript icm_extractor.R --icm ../Projects/your_project/brewpitopes/G_episurf/icm/rsa.csv --outdir ../Projects/your_project/brewpitopes/G_episurf/
```

21. Label the epitopes based on their buried positions using episurf.py
```
python3 episurf.py  
```
```
Add path to input epitopes: ../Projects/your_project/brewpitopes/F_epiglycan/glycan_extracted.csv
Add path to output folder: ../Projects/your_project/brewpitopes/G_episurf
Add path to extracted buried positions: ../Projects/your_project/brewpitopes/G_episurft/buried_positions_list.csv
```

22. Select the epitopes that are extraviral, non-glycosilated, exposed and length >= 5 using epifilter.R  
```
Rscript epifilter.R --data ../Projects/your_project/brewpitopes/G_episurf/access_extracted.csv --outdir ../Projects/your_project/brewpitopes/I_final_candidates
```

23. Extract the epitope regions using epiregions.py
```
python3 epiregions.py
```
```
Add path to input epitope dataframe: ../Projects/your_project/brewpitopes/I_final_candidates/brewpitopes_results_df.csv
Add path to output folder: ../Projects/your_project/brewpitopes/K_epitope_regions
```
24. Plot the yield results of the pipeline using yield_plot.R
```
Rscript yield_plot.R --data ../Projects/your_project/brewpitopes/H_epifilter/brewpitopes_unfiltered_df.csv --merged ../Projects/your_project/brewpitopes/D_epimerger/merged.csv --eregs ../Projects/your_project/brewpitopes/K_epitope_regions/epitope_regions_extracted.csv --outdir ../Projects/your_project/brewpitopes/J_plots
```

## APPENDIX FOR VARIANTS (OF CONCERN)
1. Generate a MUTANT FASTA using fasta_mutator.R  
      (Locally) Download reference FASTA from the target protein (in this case SARS-CoV-2 Spike) from UniprotKB and save at 
      /host/your/path/to/brewpitopes_projects/your_project/Z_fasta
      
      Save the mutations of the corresponding mutant protein or VOC as .csv at 
      /host/your/path/to/brewpitopes_projects/your_project/brewpitopes/Z_fasta
      
      Ensure the format is equal to the example file 20211203_spike_gamma_vocs.csv.
```
Rscript fasta_mutator.R --fasta ../Projects/your_project/brewpitopes/Z_fasta/target_protein.fasta --mut ../Projects/your_project/brewpitopes/Z_fasta/target_mutations.csv --mut_header your_header_without_> --sample yourfile.fasta --outdir ../Projects/your_project/brewpitopes/Z_fasta
```      
