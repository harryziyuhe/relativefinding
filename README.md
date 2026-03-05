# CSE 284 Project --- Comparing IBD Detection Methods for Relative Finding

## Project Summary
This project compares two identity-by-descent (IBD) detection methods used for identifying genetic relatives: PLINK, a genotype-based approach, and GERMLINE, a haplotype-based method. Using subsets of the 1000 Genomes Project dataset, I evaluate the performance of these two methods in terms of both accuracy and computational efficiency. PLINK estimates pairwise relatedness using genotype similarity statistics, whereas GERMLINE identifies shared haplotype segments that indicate inheritance from a common ancestor. By comparing PLINK’s IBD proportion statistics with the IBD segment lengths detected by GERMLINE, and validating results against the ground-truth pedigree information available for known relatives, this project examines how effectively the two approaches detect related individuals across different population groups.

## Dataset
The dataset used in this project comes from the 1000 Genomes Project, specifically the expanded collection containing 3,202 individuals from 26 populations. While the original dataset includes 2,504 unrelated individuals, the expanded release adds 698 individuals with known pedigree relationships, including parent–child trios. Because related individuals are relatively rare in the dataset and most populations contain only a small number of related pairs, the analysis focuses on four specific populations in order to reduce computational cost and enable clearer comparisons:
- Luhya in Kenya (LWK)
- Yoruba in Nigeria (YRI)
- Finnish in Finland (FIN)
- Southern Han Chinese (CHS)

The two African populations are included because they are frequently used as reference populations with high genetic diversity. The Finnish population represents a recently bottlenecked population, while the Southern Han Chinese population represents a large and moderately diverse East Asian population. The analysis is restricted to chromosome 10. Before running either IBD detection method, the raw data are processed to retain only biallelic SNP sites. The data are then subset by population to generate separate input files for downstream analysis.

The processing step is completed by calling
```bash
bash process.sh
```

The repository has the following file structure
```
.
├── data/                     # 1000G VCF, panel, and pedigree files
├── germline/                 # Processed data and results file for GERMLINE analysis
├── plink/                    # Processed data and results file for PLINK analysis
├── scripts/
│   ├── process.sh                     # Process raw phased VCF file 
│   ├── plink.sh                       # PLINK data preparation and analysis
│   ├── germline.sh                    # GERMLINE data preparation and analysis
│   ├── plot_analysis.ipynb            # Visualization
│   └── VCFtoHAP.py                    # Convert phased VCF data to haploid for GERMLINE analysis
├── figures/                  # Figures from analysis
└── README.md
```

## Analysis
The analysis can be executed using the following scripts:
```bash
# Run PLINK analysis
bash plink.sh

# Run GERMLINE analysis
bash germline.sh
```
The PLINK pipeline produces a ```.genome``` file containing pairwise relatedness statistics, while the GERMLINE pipeline produces a ```.match``` file containing detected IBD segments.

To compare the two methods, the notebook ```plot_analysis.ipynb``` aggregates the total shared IBD segment length for each pair of individuals detected by GERMLINE. These values are then compared with the IBD proportion estimates from PLINK, allowing a direct visualization of how the two approaches measure genetic relatedness.

## Remaining Work
- Run the GERMLINE command with different ```min_m``` size to see how sensitive the results are to parameter choices.
- Add runtime and memory comparison between the two methods.
- Make full use of the pedigree information in analysis.

## Challenges Encountered
- One challenge arises from the data format requirements of GERMLINE. The current pipeline uses the PLINK ```recode``` command to generate the ```.ped``` and ```.map``` files required for GERMLINE input. However, the GERMLINE documentation recommends first phasing genotype data using tools such as SHAPEIT, while providing limited guidance on converting phased data into the ```.ped``` and ```.map``` formats required by GERMLINE. Although the 1000 Genomes data used in this project are already phased, and the PLINK command was run with the ```--keep-allele-order``` option to preserve haplotype information, it remains possible that PLINK commands introduce inaccuracies in the GERMLINE results. To address this, I implemented a Python script (```VCFtoHAP.py```) to convert phased VCF data directly into haploid input for GERMLINE. However, due to the size of chromosome 10, this approach has encountered out-of-memory (OOM) errors when executed locally and on the DataHub environment.

- Another challenge relates to the availability of close relatives within bottlenecked populations. The original proposal aimed to compare IBD detection performance between populations with different demographic histories, particularly recently bottlenecked populations where reduced genetic diversity might increase noise in relatedness estimation. However, most bottlenecked populations in the 1000 Genomes dataset contain few or no documented close relatives, making it difficult to evaluate method accuracy in those cases. Determining how best to study the impact of population structure on relatedness inference remains an open question.

