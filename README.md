# BREWPITOPES: a pipeline to refine B-cell epitope predictions during public health emergencies

## ABSTRACT
The application of B-cell epitope identification for the development of therapeutic antibodies is well established but consuming in terms of time and resources. For this reason, in the last few years, the immunoinformatic community has developed several computational predictive tools. While relatively successful, most of these tools only use a few properties of the candidate region to determine their likelihood of being a true B-cell epitope. However, this likelihood is influenced by a wide variety of protein features, including the presence of glycosylated residues in the neighbourhood of the candidate epitope, the subcellular location of the protein region or the three-dimensional information about their surface accessibility in the parental protein.

In this study we created Brewpitopes, an integrative pipeline to curate computational predictions of B-cell epitopes by accounting for all the aforementioned features. To this end, we implemented a set of rational filters to mimic the conditions for the in vivo antibody recognition to enrich the B-cell epitope predictions in actionable candidates. To validate Brewpitopes, we analyzed the SARS-CoV-2 proteome. In the S protein, Brewpitopes enriched the initial predictions in 5-fold on epitopes with neutralizing potential (p-value < 2e-4). Other than S protein, 4 out of 16 proteins in the proteome contain curated B-cell epitopes and hence, have also potential interest for viral neutralization, since mutational escape mainly affects the S protein. Our results demonstrate that Brewpitopes is a powerful pipeline for the rapid prediction of refined B-cell epitopes during public health emergencies.

Available as preprint at: https://doi.org/10.1101/2022.11.28.518301


## INSTALLATION (DOCKER IMAGE)
0. Download the Dockerfile from this [repository](https://github.com/AlbertCS/Brewpitopes_container). To compile the Dockerfile you will need to have docker installed. 

2. Create docker image from Dockerfile (may take a while):
```
sudo docker build -t brewpitopes /host/your/path/to/Dockerfile
``` 

2. To run the Brewpitopes docker image while linking a local folder to a folder within the docker image do:

```
sudo docker run -it --volume /host/your/path/to/brewpitopes_projects:/home/Projects brewpitopes 
```

Change /host/your/path/to/brewpitopes_projects for your desired local directory.    
Once you add files/folders in this directory, they will appear automatically within the brewpitopes docker image at /home/Projects. 

3. Explore the folder structure of the Docker image named "brewpitopes".
```
cd ..
ls
```
You should see the folders: Brewpitopes ## revise

3.1 To execute the scripts of the pipeline go to:
```      
cd Brewpitopes
```

3.2 To find an example of a brewpitopes project and the required files go to
```
cd ../example
```

3.3 To create new project, do so within the Projects folder:
```
cd ../Projects
mkdir your_project
cd your_project
```

3.4 You are now ready to predict and refine your B-cell epitope candidates.

## PIPELINE
** All the steps should be run within the Docker image at your Terminal except when it indicates "(Locally)". **

4. Use directories.R to create the directory environment in your project folder.

```
cd Projects/your_project
Rscript ../../directories.R
cd brewpitopes
```

5. (Locally) Download the FASTA file of the target protein at [Uniprot](https://www.uniprot.org/).        
      Save at Z_fasta  ##(Check if it requires SUDO permissions when using Docker).       
      
6. (Locally) Use the FASTA to predict linear epitopes using [Bepipred 2.0](https://services.healthtech.dtu.dk/service.php?BepiPred-2.0) server and export results as csv (default parameters).  
       
      Save at A_linear_predictions/bepipred/bepipred_results.csv  
      If Bepipred2.0 does not predict any epitope in your target sequence see Note1.      
      
7. Extract epitopes from Bepipred results using epixtractor_linear_bepipred.py.  
```
python3 ../../../epixtractor_linear_bepipred.py
```
```
Add path to bepipred results: A_linear_predictions/bepipred/bepipred_results.csv
Add path to output folder: C_epixtractor    
```
      
8. (Locally) Use the FASTA to predict linear epitopes using [ABCpred](https://webs.iiitd.edu.in/raghava/abcpred/ABC_submission.html) server.

      Predict using all the epitope windows (10,12,14,16,18,20), overlapping filter ON and the default threshold at 0.51.     
      Copy results from the webpage table to a text editor (See Note 2 and 3).      
      Save as: A_linear_predictions/abcpred/abcpred_allmers.tsv         
      If ABCpred does not predict any epitope in your target sequence see Note1.    

      
9. Extract epitopes from ABCpred results using epixtractor_linear_abcpred.R  

````
Rscript ../../../epixtractor_linear_abcpred.R --outdir C_epixtractor --input_allmers A_linear_predictions/abcpred/abcpred_allmers.tsv 
````
      
10. (Locally) Download the PDB file of the target protein at PDB DB. 
       Save at ZZ_pdb/pdb  
       
11. (Locally) Use [PDBrenum](http://dunbrack3.fccc.edu/PDBrenum/) server to renumerate the PDB residues according to its corresponding FASTA file in Uniprot.  
      Download results as .pdb  (Selecting the .pdb options and deselect .mmCIF options)  
      Uncompress the .zip file donwloaded.      
      Uncompress the file your_pdb_id.pdb.gz    
      Save your_pdb_id.pdb at ZZ_pdb/pdbrenum         
      
12. (Locally) Use the renumbered PDB to predict structural epitopes using [Discotope 2.0](https://services.healthtech.dtu.dk/service.php?DiscoTope-2.0) server (input type 3).

      Select default score threshold (-3.7).    
      Select chain A by default or your target chain of interest.       
      Export the results as .txt. Remove the last line "Identified...". Then, save as .csv by replacing "\t" for commas.      
      Save as B_structural_predictions/discotope/discotope_results.csv        
      If Discotope2.0 does not predict any epitope in your target sequence see Note1.     
      
13. Extract epitopes from Discotope results using epixtract_structural.py
```
python3 ../../../epixtract_structural.py  
```
```
Add path to discotope results: B_structural_predictions/discotope/discotope_results.csv
Add path to output folder: C_epixtractor
```

14. Merge the epitopes extracted from Bepipred, ABCpred and Discotope results using epimerger.R
```
Rscript ../../../epimerger.R --abcpred C_epixtractor/abcpred_results_extracted.csv --bepipred C_epixtractor/bepipred_results_extracted.csv --discotope C_epixtractor/discotope_results_extracted.csv --outdir D_epimerger
```

Take steps 15, 16 and 17.1 if you want to predict protein topology using CCTOP. Otherwise, if the topology of your protein is already described and you want to add it manually go directly to step 17.2.

15. (Locally) Predict the protein topology using [CCTOP](http://cctop.enzim.ttk.mta.hu/?_=/jobs/submit) server.  
      Donwload results as .xml.
      Save at E_topology/CCTOP/cctop.xml
      
16. Extract the topological domains using xml_cctop_parser.R. (See Note 4)
```
Rscript ../../../xml_cctop_parser.R --xml E_epitopology/CCTOP/cctop.xml --outdir E_epitopology/CCTOP
```

17. Label the epitopes based on their topology (intracellular, membrane or extracellular) using epitopology.R

17.1 Using CCTOP predictions --> use epitopology_cctop.R
```
Rscript ../../../epitopology_cctop.R --input_CCTOP E_epitopology/CCTOP/cctop_domains.csv --input_epitopes D_epimerger/merged.csv --outdir E_epitopology
```

17.2 Manual topology annotation --> use epitopology_manual.R. Add manually the starting and ending positions of your extracellular domains at --start_pos and end_pos.
```
Rscript ../../../epitopology_manual.R --start_pos 1,12,22 --end_pos 8,18,28 --input_epitopes D_epimerger/merged.csv --outdir E_epitopology
```
      
18. (Locally) Predict the glycosilation profile of the protein using the FASTA file.  
      N-GLYCOSYLATIONS at [NetNGlyc 1.0](https://services.healthtech.dtu.dk/service.php?NetNGlyc-1.0) server.    
      Copy manually into a spreadsheet the output table headed: SeqName	Position	Potential	Jury_agreement	NGlyc_result          	 
      SAVE AS CSV as F_epiglycan/netnglyc/nglyc_results.csv 
      (See Note 5)      
      
      O-GLYCOSYLATIONS AT [NetOGlyc 4.0](https://services.healthtech.dtu.dk/service.php?NetOGlyc-4.0) server.
      Copy manually into a spreadsheet the output rows of table headed: seqName  	source	feature	start 	end	score strand      frame       comment. Do NOT include the header, only the data.           
      SAVE AS CSV at F_epiglycan/netoglyc/oglyc_results.csv       
      (See Note 5)
      
19. Extract the glycosylated positions from both N-glyc and O-glyc outputs using epiglycan_extractor.R
```
Rscript ../../../epiglycan_extractor.R --oglyc F_epiglycan/netoglyc/oglyc_results.tsv --nglyc F_epiglycan/netnglyc/nglyc_results.tsv --outdir F_epiglycan
```

20. Use epiglycan.py to label the glycosylated epitopes.  
```
python3 ../../../epiglycan.py
```
```
Add path to input epitopes: E_epitopology/topology_extracted.csv
Add path to output folder: F_epiglycan
Add path to extracted glycosylated positions: F_epiglycan/glycan_positions.csv  
```

21. (Locally) Use ICM_browser (MOLSOFT) to extract the RSA values for accessibility calculation.  
      Download ICM_browser from [http://www.molsoft.com/icm_browser.html](http://www.molsoft.com/icm_browser.html)      
      Move the PDB renumbered file of the corresponding protein (step 11) to a local folder.    
      Move the script Compute_ASA.icm to the same local folder.   
      Execute in the command line of the programme ICM Browser the code in Compute_ASA.icm        
      Save results at G_episurf     
      
22. Extract the buried positions using icm_extractor.R  
```
Rscript ../../../icm_extractor.R --icm G_episurf/icm/rsa.csv --outdir G_episurf
```

21. Label the epitopes based on their buried positions using episurf.py
```
python3 ../../../episurf.py
```
```
Add path to input epitopes: F_epiglycan/glycan_extracted.csv
Add path to output folder: G_episurf
Add path to extracted buried positions: G_episurf/buried_positions_list.csv
```

22. Select the epitopes that are extraviral, non-glycosilated, exposed and length >= 5 using epifilter.R  
```
Rscript ../../../epifilter.R --data G_episurf/access_extracted.csv --outdir I_final_candidates
```

23. Extract the epitope regions using epiregions.py
```
python3 ../../../epiregions.py
```
```
Add path to input epitope dataframe: I_final_candidates/brewpitopes_results_df.csv
Add path to output folder: K_epitope_regions
```
24. Plot the yield results of the pipeline using yield_plot.R
```
Rscript ../../../yield_plot.R --data H_epifilter/brewpitopes_unfiltered_df.csv --merged D_epimerger/merged.csv --eregs K_epitope_regions/epitope_regions_extracted.csv --outdir J_plots
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

## NOTES
1. In the case one of the predictor softwares (Bepipred, ABCpred or Discotope) does not identify any epitope in your target protein, do create the empty dataframe containing only the required headers. It will be used later in the pipeline.
2. When copying the tabular results from ABCpred into a text file make sure to remove all lines with no sequence information (see example). Otherwise, the script "epixtractor_linear_abcpred.R" will not be able to read the file and will prompt the following error: "Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec, : line XXX did not have 5 elements ". 
![ABCpred results processing](/abcpred_processing.png?raw=true "ABCpred processing")
3. When copying the tabular results from ABCpred into a text file make sure to put a tab at the end of the last line and follow it by a blank line below.
4. When predicting the protein topology using CCTOP you might encounter no extracellular regions and you should not continue the Brewpitopes pipeline with this target protein. In such, case you will get the error: "STOPPER!! Your target protein has no predicted extracellular domains. Hence, neutralizing antibodies will not recognize it. You should consider another protein from your target organism."
5. You can use spreadsheet editors like Microsoft Excel or Libre Office.

