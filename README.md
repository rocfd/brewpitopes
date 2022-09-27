# BREWPITOPES
Set of tools to manage epitope prediction results from linear and structural origin and to integrate a pipeline of prioritization filters to curate B-cell epitopes specific for neutralizing antibody recognition.

## INSTALLATION (DOCKER IMAGE)
To compile the Dockerfile you will need to have docker installed. And use the following commands:
1. Create docker image from Dockerfile (may take a while):
```
      sudo docker build -t brewpitopes PATH/TO/Dockerfile
```
2. Run the docker image:
```
      sudo docker run -it brewpitopes
```      
3. Create a shared folder between Brewpitopes docker image and your local machine.
```
      sudo docker run -it brewpitopes --volume /your/machine/directory:/home/Projects
```

## PIPELINE
4. Use directories.R to create the folder environment.

      ```
      Rscript directories.R --path /your/desired/folder
      ```

5. Download the FASTA file of the target protein at [Uniprot](https://www.uniprot.org/).    
      
      Save at /Z_fasta  
      
6. Use the FASTA to predict linear epitopes using [Bepipred 2.0] (https://services.healthtech.dtu.dk/service.php?BepiPred-2.0) server and export results as csv (default parameters).  
       
      Save at /A_linear_predictions/bepipred/bepipred_results.csv  
      
7. Extract epitopes from Bepipred results using epixtractor_linear_bebipred.py.  
```
python3 epixtractor_linear_bebipred.py
```
```
Add path to bepipred results: your/path/to/A_linear_predictions/bepipred/bepipred_results.csv
Add path to output folder: your/path/to/C_epixtractor    
```
      
6. Use the FASTA to predict linear epitopes using [ABCpred](https://webs.iiitd.edu.in/raghava/abcpred/ABC_submission.html) server.

      Predict using all the epitope windows (10,12,14,16,18,20) and overlapping filter ON.  
      Copy results from the webpage to a .csv  
      Save at: path/to/brewpitopes/A_linear_predictions/abcpred/abcpred_10mers.csv 
      
7. Extract epitopes from ABCpred results using epixtractor_linear_abcpred.R  

````
Rscript epixtractor_linear_abcpred.R --outpath your/path/to/brewpitopes/C_epixtractor --input_10mers your/path/to/brewpitopes/A_linear_predictions/abcpred/abcpred_10mers.csv --input_12mers your/path/to/brewpitopes/A_linear_predictions/abcpred/abcpred_12mers.csv --input_14mers your/path/to/brewpitopes/A_linear_predictions/abcpred/abcpred_14mers.csv --input_16mers your/path/to/brewpitopes/A_linear_predictions/abcpred/abcpred_16mers.csv --input_18mers your/path/to/brewpitopes/A_linear_predictions/abcpred/abcpred_18mers.csv --input_20mers your/path/to/brewpitopes/A_linear_predictions/abcpred/abcpred_20mers.csv
````
      
8. Download the PDB file of the target protein at PDB DB. 
       Save at /brewpitopes/B_structural_predictions/pdb  
       
9. Use [PDBrenum](http://dunbrack3.fccc.edu/PDBrenum/) server to renumerate the PDB residues according to its corresponding FASTA file in Uniprot.  
      Download results as .pdb  
      Save at /brewpitopes/B_structural_predictions/pdbrenum 
      
10. Use the renumbered PDB to predict structural epitopes using [Discotope 2.0](https://services.healthtech.dtu.dk/service.php?DiscoTope-2.0) server and export the results as csv.    
      Default threshold.  
      Select chain A by default.
      Save at /brewpitopes/B_structural_predictions/discotope  
      
11. Extract epitopes from Discotope results using epixtract_structural.py
```
python3 epixtract_structural.py  
```
```
Add path to discotope results: brewpitopes/B_structural_predictions/discotope/discotope_results.csv  
Add path to output folder: brewpitopes/C_epixtractor
```

12. Merge the epitopes extracted from Bepipred, ABCpred and Discotope results using epimerger.R
```
Rscript epimerger.R --abcpred your/path/to/brewpitopes/C_epixtractor/abcpred_results_extracted.csv --bepipred your/path/to/brewpitopes/C_epixtractor/abcpred_results_extracted.csv --discotope your/path/to/brewpitoeps/C_epixtractor/discotope_results_extracted.csv --outdir your/path/to/brewpitoeps/D_epimerger
```

13. Predict the protein topology using [CCTOP](http://cctop.enzim.ttk.mta.hu/?_=/jobs/submit) server.  
      Donwload results as .xml.
      Save at your/path/to/brewpitopes/E_topology/CCTOP/cctop.xml
      
14. Extract the topological domains using xml_cctop_parser.R  
```
Rscript xml_cctop_parser.R --xml path/to/brewpitopes/E_epitopology/CCTOP/cctop.xml --outdir path/to/brewpitopes/E_epitopology/CCTOP
```

15. Label the epitopes based on their topology (intracellular, membrane or extracellular) using epitopology.
Using CCTOP predictions --> use epitopology_cctop.R
```
Rscript epitopology_cctop.R --input_CCTOP path/to/brewpitopes/E_epitopology/CCTOP/cctop_domains.csv --input_epitopes path/to/brewpitopes/D_epimerger/merged.csv --outdir path/to/brewpitopes/E_epitopology
```

Using manual annotation --> use epitopology_manual.R
```
Rscript epitopology_manual.R --start_pos 1,12,22 --end_pos 8,18,28 --input_epitopes path/to/brewpitopes/D_epimerger/merged.csv --outdir path/to/brewpitopes/E_epitopology
```
      
16. Predict the glycosilation profile of the protein using the FASTA file.  
      N-GLYCOSILATIONS at [NetNGlyc 1.0](https://services.healthtech.dtu.dk/service.php?NetNGlyc-1.0) server.    
      COPY MANUALLY THE DATAFRAME HEADED: SeqName	Position	Potential	Jury_agreement	NGlyc_result	Prediction  
      SAVE AS CSV at brewpitopes/F_epiglycan/netnglyc  
      
      O-GLYCOSILATIONS AT [NetOGlyc 4.0](https://services.healthtech.dtu.dk/service.php?NetOGlyc-4.0) server.
      COPY MANUALLY THE DATAFRAME HEADED: seqName  	source	feature	start 	end	score strand      frame       comment  
      SAVE AS CSV at brewpitopes/F_epiglycan/netoglyc
      
17. Extract the glycosilated positions from both N-glyc and O-glyc outputs using epiglycan_extractor.R
```
Rscript epiglycan_extractor.R --oglyc /your/path/to/brewpitopes/F_epiglycan/netoglyc/oglyc.csv --nglyc /your/path/to/brewpitopes/F_epiglycan/netnglyc/nglyc.csv --outdir brewpitopes2/F_epiglycan/
```

18. Use epiglycan.py to label the glycosilated epitopes.  
```
python3 epiglycan.py
```
```
Add path to input epitopes: brewpitopes/E_epitopology/topology_extracted.csv
Add path to output folder: brewpitopes/F_epiglycan
Add path to extracted glycosilated positions: brewpitopes/F_epiglycan/glycosilated_positions.csv  
```

19. Use ICM_browser (MOLSOFT) to extract the RSA values for accessibility calculation.  
      Download ICM_browser from [http://www.molsoft.com/icm_browser.html](http://www.molsoft.com/icm_browser.html)
      Open the PDB renumbered file of the corresponding protein (step 9).  
      Execute in the command line of the programme the code in Compute_ASA.icm  
      Save results at /G_episurf
      
20. Extract the buried positions using icm_extractor.R  
```
Rscript icm_extractor.R --icm /your/path/to/brewpitopes/G_episurf/icm/rsa.csv --outdir your/path/to/brewpitopes/G_episurf/
```

21. Label the epitopes based on their buried positions using episurf.py
```
python3 episurf.py  
```
```
Add path to input epitopes: brewpitopes/F_epiglycan/glycan_extracted.csv
Add path to output folder: brewpitopes/G_episurf
Add path to extracted buried positions: brewpitopes/G_episurft/buried_positions_list.csv
```

22. Use epifilter.R to retain the epitopes that are extraviral, non-glycosilated, exposed and length >= 5.  
      Follow the R file instructions.  
      Save at /I_final_candidates.  
23. Use epicontig.ipynb (Jupiter Notebook) to extract the epitopic regions / contigs.  
      Upload the candidates_df.csv generated at step 22.
      Follow the instructions in the Notebook.  
24. Use yield_plot.R to plot the results of the pipeline.  
      Follow the instructions in the R file.  
      
## APPENDIX FOR VARIANTS OF CONCERN
1. Generate FASTA using fasta_mutator.R  
      Download reference FASTA from Spike protein from UniprotKB.  
      Upload where indicated at script instructions.  
      Upload the mutations of the corresponding VOC found as attached files in this Github. (ie Gamma = 20211203_spike_gamma_vocs.csv)
      Execute the script and save the VOC Fasta file.  
      Once saved, remove "" from the file to obtain a properly formatted FASTA.  
      Start the pipeline above with the mutated FASTA file.  
