# BREWPITOPES
Set of tools to manage epitope prediction results from linear and structural origin and to integrate a pipeline of prioritization filters.

## PIPELINE
1. Create and R project for the corresponding protein.
2. Use directories.R to create the folder environment.
3. Download the FASTA file of the target protein at Uniprot.
      Save at /Z_fasta
4. Use the FASTA to predict linear epitopes using Bepipred and export results as csv (default parameters).
      https://services.healthtech.dtu.dk/service.php?BepiPred-2.0
            Save at path/to/A_linear_predictions/bebipred
5. Extract epitopes from Bepipred results using epixtractor_linear_bebipred.py.
      python3 epixtractor_linear_bebipred.py
            path/to/bepipred_results.csv
                  Copy output to /C_epixtractor.
6. Use the FASTA to predict linear epitopes using ABCpred.
      https://webs.iiitd.edu.in/raghava/abcpred/ABC_submission.html
      Predict using all the epitope windows and overlapping filter ON.
      Copy results and save as .csv
      Save at path/to/A_linear_predictions/abcpred
7. Extract epitopes from ABCpred results using epixtractor_linear_abcpred.R
      Add paths to ABCpred result files.
      Follow the instructions in the R file.
      Copy output to /C_epixtractor
8. Download the PDB file of the target protein at PDB DB. 
       Save at /B_structural_predictions/pdb
10. Use PDBrenum to renumerate the PDB residues according to its corresponding FASTA file in Uniprot.
      http://dunbrack3.fccc.edu/PDBrenum/
      Save results as .pdb
      Save at /B_structural_predictions/pdbrenum
11. Use the renumbered PDB to predict structural epitopes using Discotope and export the results as csv.
      https://services.healthtech.dtu.dk/service.php?DiscoTope-2.0
      Default threshold.
      Select chain A.
      Save at /B_structural_predictions/discotope
12. Extract epitopes from Discotope results using epixtract_structural.py
      python3 epixtract_structural.py
      path/to/B_structural_predictions/discotopediscotope_results.csv
13. Use epimerger.py to merge the epitopes extracted from Bepipred, ABCpred and Discotope results.
      Add the paths to the corresponding files at C_epixtractor folder.
      Follow the instructions in the R file.
      Save at /D_epimerger
14. Predict the protein topology using CCTOP.
      http://cctop.enzim.ttk.mta.hu/?_=/jobs/submit
      Donwload results as .xml
      Save at /E_topology/CCTOP
15. Extract the extraviral domains using xml_cctop_parser.R
      Follow the instructions at the R file.
16. Use epitopology.R to label the epitopes based on the extraviral domains.
      Follow the instructions at the R file. 
      Save results at /E_topology
17. Predict the glycosilation profile of the protein using the FASTA file.
      N-GLYCOSILATIONS AT:
      https://services.healthtech.dtu.dk/service.php?NetNGlyc-1.0
      SAVE THE DATAFRAME HEADED SeqName	Position	Potential	Jury_agreement	NGlyc_result	Prediction
      AS CSV
      O-GLYCOSILATIONS AT:
      https://services.healthtech.dtu.dk/service.php?NetNGlyc-1.0
      SAVE THE DATAFRAME HEADED seqname	source	feature	start	end	score	strand	frame	comment
      AS CSV
18. Use epiglycan_extractor.R to extract the glycosilated positions.
      Follow the instructions at the R file.
      Save at /F_epiglycan
19. Use epiglycan.py to label the glycosilated epitopes.
      python3 epiglycan.py
      /E_epitopology/XXX_epitopology.csv
      /F_epiglycan/glycosilated_positions.csv
      Save as /F_epiglycan/XXX_epitopology_glycans.csv
20. Use ICM_browser (MOLSOFT) to extract the RSA values for accessibility calculation.
      Download ICM_browser from http://www.molsoft.com/icm_browser.html
      Open the PDB renumbered file of the corresponding protein.
      Execute in the command line the code in Compute_ASA.icm
      Save results at /G_episurf
21. Use icm_extractor.R to extract the buried positions.
      Follow the instructions at the R file. 
      Save at /G_episurf
22. Use episurf.py to label the epitopes with buried residues.
      python3 episurf.py
      path/to/F_epiglycan/XXX_epitoplogy_glycans.csv
      path/to/G_episurft/buried_positions.csv
      Save as /G_episurf/XXX_epitoplogy_glycans_surf.csv
23. Use epifilter.R to filter out the epitopes that are extraviral, non-glycosilated, exposed and length >= 5.
      Follow the R file instructions.
      Save at /I_final_candidates.
24. Use epicontig.ipynb (Jupiter Notebook) to extract the epitopic regions / contigs.
      Follow the instructions in the Notebook.
25. Use yield_plot.R to plot the results of the pipeline.
      Follow the instructions in the R file.
