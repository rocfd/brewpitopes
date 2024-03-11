# BREWPITOPES: a pipeline to refine B-cell epitope predictions during public health emergencies
Update based on WNV project
## ABSTRACT
The application of B-cell epitope identification for the development of therapeutic antibodies is well established but consuming in terms of time and resources. For this reason, in the last few years, the immunoinformatic community has developed several computational predictive tools. While relatively successful, most of these tools only use a few properties of the candidate region to determine their likelihood of being a true B-cell epitope. However, this likelihood is influenced by a wide variety of protein features, including the presence of glycosylated residues in the neighbourhood of the candidate epitope, the subcellular location of the protein region or the three-dimensional information about their surface accessibility in the parental protein.

In this study we created Brewpitopes, an integrative pipeline to curate computational predictions of B-cell epitopes by accounting for all the aforementioned features. To this end, we implemented a set of rational filters to mimic the conditions for the in vivo antibody recognition to enrich the B-cell epitope predictions in actionable candidates. To validate Brewpitopes, we analyzed the SARS-CoV-2 proteome. In the S protein, Brewpitopes enriched the initial predictions in 5-fold on epitopes with neutralizing potential (p-value < 2e-4). Other than S protein, 4 out of 16 proteins in the proteome contain curated B-cell epitopes and hence, have also potential interest for viral neutralization, since mutational escape mainly affects the S protein. Our results demonstrate that Brewpitopes is a powerful pipeline for the rapid prediction of refined B-cell epitopes during public health emergencies.

Available as preprint at: https://doi.org/10.1101/2022.11.28.518301


## INSTALLATION (DOCKER IMAGE)

1. Download the docker image from [Dockerhub](https://hub.docker.com/r/bsceapm/brewpitopes) (it might take a while). You need to have Docker installed.
```
docker pull bsceapm/brewpitopes:1.0
```

2. To run the Brewpitopes docker image while linking a local folder to a folder within the docker image do:

```
sudo docker run -it --volume /host/your/path/to/brewpitopes_projects:/home/brewpitopes/Projects brewpitopes
```

Change /host/your/path/to/brewpitopes_projects for your desired local directory.    
Once you add files/folders in this directory, they will appear automatically within the brewpitopes docker image at /home/brewpitopes/Projects.

3. Explore the folder structure of the Docker image named "brewpitopes".
```
cd ..
ls
```
You should see the folder: brewpitopes

3.1 To find an example of a brewpitopes project and the required files go to
```
cd brewpitopes/example
```

3.2 To create new project, do so within the Projects folder:
```
cd ../Projects
mkdir your_project
cd your_project
```

3.3 You are now ready to predict and refine your B-cell epitope candidates.

## PIPELINE
** All the steps should be run within the Docker image at your Terminal except when it indicates "(Locally)". **

4. Use directories.R to create the directory environment in your project folder.

```
cd Projects/your_project
Rscript ../../directories.R
cd brewpitopes
```

5. (Locally) Download the FASTA file of the target protein at [Uniprot](https://www.uniprot.org/).        
      Save at Z_fasta               

6. (Locally) Use the FASTA to predict linear epitopes using [Bepipred 3.0](https://services.healthtech.dtu.dk/services/BepiPred-3.0/) server and export results (as csv)zip file). Use default options (High-confidence, Thr: 0.152).  

      Save the file Bcell_linepitope_top_20pct_preds.fasta (within zip) at A_linear_predictions/bepipred/bepipred_linear_pred.fasta 
      If Bepipred3.0 does not predict any epitope in your target sequence see Note1.      

7. Extract epitopes from Bepipred results using epixtractor_linear_bepipred3.py. This will automatically save the results at C_epixtractor
```
python3 ../../../epixtractor_linear_bepipred3.py --path /pathto/project/brewpitope
```

8. (Locally) Use the FASTA to predict linear epitopes using [ABCpred](https://webs.iiitd.edu.in/raghava/abcpred/ABC_submission.html) server.

      Predict using all the epitope windows (10,12,14,16,18,20), overlapping filter ON and the default threshold at 0.51.     
      Copy results from the webpage table to a text editor (See Note 2 and 3).      
      Save as: A_linear_predictions/abcpred/abcpred_allmers.tsv         
      If ABCpred does not predict any epitope in your target sequence see Note1.    


9. Extract epitopes from ABCpred results using epixtractor_linear_abcpred.R  

```
Rscript ../../../epixtractor_linear_abcpred.R --outdir C_epixtractor --input_allmers A_linear_predictions/abcpred/abcpred_allmers.tsv
```

10. (Locally atm). Install EpitopeVec from its repo (https://github.com/hzi-bifo/epitope-prediction)


11. (Locally atm) Compute the epitope results from EpitopeVec for num_pep 10 to 20 (i.e 20 files)

     Save the outputs at A_linear_predictions/epitopeveclinear_pred_{len_peptide}.csv
     E.g A_linear_predictions/epitopeveclinear_pred_11.csv
         A_linear_predictions/epitopeveclinear_pred_17.csv

12. (Locally atm) Extract epitopes from EpitopeVec results using epixtractor_linear_epitopevec.R

```
Rscript ../../../epixtractor_linear_epitopevec.R --outdir C_epixtractor
```

13. (Locally) Download the PDB file of the target protein at PDB DB.
       Save at ZZ_pdb/pdb  

14. (Locally) Use [PDBrenum](http://dunbrack3.fccc.edu/PDBrenum/) server to renumerate the PDB residues according to its corresponding FASTA file in Uniprot.  
      Download results as .pdb  (Selecting the .pdb options and deselect .mmCIF options)  
      Uncompress the .zip file donwloaded.      
      Uncompress the file your_pdb_id.pdb.gz    
      Save your_pdb_id.pdb at ZZ_pdb/pdbrenum   
      In the case you are using an Alphafold model, you will not need to renumber the pdb.

15. (Locally) Use the renumbered PDB to predict structural epitopes using [Discotope 3.0](https://services.healthtech.dtu.dk/services/DiscoTope-3.0/) server.

      Select pertinent input PDB option (i.e AF or resolved).
      Select Higher Confidence option
      Download .zip from output
      Copy the .pdb and .csv of desired Chain to B_structural_predictions/discotope
      Rename files to B_structural_predictions/discotope/discotope_pred.csv
                      B_structural_predictions/discotope/discotope_pred.pdb 
      If Discotope3.0 does not predict any epitope in your target sequence see Note1.     

16. (Locally) Extract the residues predicted by Discotope3.0, using a surface-based proximity approach
```
     It is important to have the script prot_surface_cluster.py on same folder
     as epixtractor_structural_discotope3.py
     You will need to have MSMS installed (see script instructions and Note XX)
     
     python3 ../../../epixtractor_structural_discotope3.py --path /pathto/project/brewpitope
     
```

17. (Locally) Use the FASTA file to predict structural epitopes using [Bepipred 3.0](https://services.healthtech.dtu.dk/services/BepiPred-3.0/) 
server and export results (as csv)zip file). Use default options (High-confidence, Thr: 0.152).

      Save the file raw_output.csv  (within zip) at 
      B_structural_predictions/bepipred3_conf/bepipred3_conf_output.csv
      Copy the renumbered PDB at B_structural_predictions/bepipred3_conf/brepitope_renumpdb.pdb
      If Bepipred3.0 does not predict any epitope in your target sequence see Note1. 

18. (Locally) Extract the residues predicted by Bepipred3.0, using a surface-based proximity approach
```
     It is important to have the script prot_surface_cluster.py on same folder
     as epixtractor_structural_bepipred3.py
     You will need to have MSMS installed (see script instructions and Note XX)
     
     python3 ../../../epixtractor_structural_bepipred3.py --path /pathto/project/brewpitope
```

19. (Locally) Use the PDB at the [SEPPA3.0](https://services.healthtech.dtu.dk/services/BepiPred-3.0/) server
    Select the parameters accordingly to your system.
    
    Download the output .pdb and copy the file to B_structural_predictions/seppa/seppa_predict.pdb
    If Seppa3.0 does not predict any epitope in your target sequence see Note1.


20. (Locally) Extract the residues predicted by Seppa3.0, using a surface-based proximity approach

     It is important to have the script prot_surface_cluster.py on same folder
     as epixtractor_structural_seppa3.py
     You will need to have MSMS installed (see script instructions and Note XX)
     
      python3 ../../../epixtractor_structural_seppa3.py --path /pathto/project/brewpitope


21. (Locally) Use the FASTA file at the [SeRenDIP](https://www.ibi.vu.nl/programs/serendipwww/) server
    Select the RF model of "Epitopes".
    
     Download (copy-paste)  the output and copy the file to B_structural_predictions/serendipce/serendipce_predict.csv
     Copy the renumbered PDB at B_structural_predictions/serendipce/serendipce_renumpdb.pdb
     If SeRendip-CE does not predict any epitope in your target sequence see Note1.
     
22. (Locally) Extract the residues predicted by SeRenDIP-CE, using a surface-based proximity approach

     It is important to have the script prot_surface_cluster.py on same folder
     as epixtractor_structural_serendip.py
     You will need to have MSMS installed (see script instructions and Note XX)
     
     python3 ../../../epixtractor_structural_serendipce.py --path /pathto/project/brewpitope


23. Merge the epitopes extracted from different epitopes methods using epimerger.R
```
Rscript ../../../epimerger.R --path /pathto/project/brewpitope
```

Take steps 15, 16 and 17.1 if you want to predict protein topology using CCTOP. Otherwise, if the topology of your protein is already described and you want to add it manually go directly to step 17.2.

24. (Locally) Predict the protein topology using [CCTOP](http://cctop.enzim.ttk.mta.hu/?_=/jobs/submit) server.  
      Copy the .xml output from the main page into a text file (do not copy the downloadable file).    
      Save as E_topology/CCTOP/cctop.xml

25. Extract the topological domains using xml_cctop_parser.R. (See Note 4)
```
Rscript ../../../xml_cctop_parser.R --path /pathto/project/brewpitope
```

26. Label the epitopes based on their topology (intracellular, membrane or extracellular) using epitopology.R

26.1 Using CCTOP predictions --> use epitopology_cctop.R
```
Rscript ../../../epitopology_cctop.R --path /pathto/project/brewpitope
```

26.2 Manual topology annotation --> use epitopology_manual.R. Add manually the starting and ending positions of your extracellular domains at --start_pos and end_pos.
```
Rscript ../../../epitopology_manual.R --start_pos 1,12,22 --end_pos 8,18,28 --input_epitopes D_epimerger/merged.csv --outdir E_epitopology
```

27. (Locally) Predict the Post-translational modifications profile of the protein using the FASTA file.  
      
      N-GLYCOSYLATIONS at [NetNGlyc 1.0](https://services.healthtech.dtu.dk/service.php?NetNGlyc-1.0) server.       
      Copy manually into a text editor the output table headed: SeqName	Position	Potential	Jury_agreement	NGlyc_result      
      Do NOT include the header(error prone), only the data.            
      SAVE AS CSV as F_epitpm/netnglyc/nglyc_predict.csv       
      (If you get an error, please see Note 5)      

      O-GLYCOSYLATIONS AT [NetOGlyc 4.0](https://services.healthtech.dtu.dk/service.php?NetOGlyc-4.0) server.     
      Copy manually into a text editor the rows of the output table headed: seqName  	source	feature	start 	end	score strand      frame       comment.          
      Do NOT include the header(error prone), only the data.               
      SAVE AS CSV at F_epitpm/netoglyc/oglyc_predict.csv       
      (If you get an error, please see Note 5)

      C-mannosylation Glycosilation AT [NetCGlyc 1.0](https://services.healthtech.dtu.dk/services/NetCGlyc-1.0/) server.
      Use the GFF output format. Copy the output table manually to a text editor
      Copy manually into a text editor the rows of the output        
      Do NOT include the header(i.e rows starting with #), only the data.               
      Save file as F_epitpm/netcglyc/netcglyc_predict.csv       
      (If you get an error, please see Note 5)
      
      Phosphorylation sites using [NetPhos 3.1](https://services.healthtech.dtu.dk/services/NetPhos-3.1/) server
      Copy manually the output GFF table  into a text editor.
      Skip header rows (i.e ALL starting with #)
      Save file as: F_epitpm/netphos/phospho_predict.csv    
      (If you get an error, please see Note 5)  
      
      
28. Extract the PTM positions from N-glyc, O-glyc, C-glyc and Phospho outputs using epiptm_extractor.R
```
Rscript ../../../epiptm_extractor.R --path /pathto/project/brewpitope
```

29. Use epiptm.R to label the Post-translational modifications epitopes.  
```
Rscript ../../../epiptm.R --path /pathto/project/brewpitope
```

30. Label the epitopes based on their buried positions using episurf.py

      We will compute the SASA of each residue, and label the epitopes
      based on theri >80% of residues accesibility
      (If you get an error, please see Note 5)
      
      PDB file must be copied at "G_episurf/episurf.pdb"
      
```
python ../../../episurf.py --path /pathto/project/brewpitope
```

31. Select the epitopes that are extraviral, non-glycosilated, exposed and length >= 5 using epifilter.R  
```
Rscript ../../../epifilter.R --path /pathto/project/brewpitope
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
5. You can use any text editor but make sure to place an empty line at the end of the TSV file. Otherwise, you might get an error due to incomplete final line.
