## MAKE DIRECTORIES FOR BREWPITOPES

library(argparser)

### SCRIPT ARGUMENTS
# Create a parser
p <- arg_parser("PROJECT DIRECTORIES")

# Add command line arguments
p <- add_argument(p, "--path", help= "Path to project directory", type="character", default = ".")

# Parse the command line arguments
argv <- parse_args(p)

### CREATE PROJECT DIRECTORY
dir.create(argv$path)

### BREWPITOPES DIRECTORY
dir.create(paste0(argv$path, "/", "brewpitopes", sep = ""))

## INITIAL DATA - FASTA
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "Z_fasta", sep = ""))

## DIRECTORIES
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "A_linear_predictions", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "A_linear_predictions/bebipred", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "A_linear_predictions/abcpred", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "B_structural_predictions", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "B_structural_predictions/discotope", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "ZZ_pdb/pdb", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "ZZ_pdb/alphafold", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "ZZ_pdb/pdbrenum", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "C_epixtractor", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "D_epimerger", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "E_epitopology", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "E_epitopology/CCTOP", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "F_epiglycan", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "F_epiglycan/netnglyc", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "F_epiglycan/netoglyc", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "G_episurf", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "G_episurf/icm", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "H_epifilter", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "I_final_candidates", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "J_plots", sep = ""))
dir.create(paste0(argv$path, "/", "brewpitopes", "/", "K_epitope_regions", sep = ""))
