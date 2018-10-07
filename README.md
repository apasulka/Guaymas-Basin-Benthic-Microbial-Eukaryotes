## Bioinformatic workflow - Pasulka *et al.*
Methods used for quality control, clustering, and analysis of tag-sequencing (18S rRNA gene) survey from Guaymas Basin benthic samples.

#### Alexis Pasulka, Sarah K. Hu, Lisa Mesrop, Craig Cary, Kathryn Coyne, Karla Heidelberg, Peter Countway, & David A. Caron. SSU-rRNA sequencing survey of sediment-hosted microbial eukaryotes from Guaymas Basin hydrothermal vent. _In prep._

## Sequence quality control
Benthic samples were collected from Guaymas Basin, MX, from a hydrothermally active 
Raw sequences are from NCBI (SRA) under project ID SRP110312 or accession numbers SAMN07274333 - SAMN07274355.

Script for quality filtering sequences and generating an OTU table: 'SeqQC_OTUclustering_Pasulka_et_al.pl'

## Data analysis

Import OTU table into R for all downstream analyses.

* Initial OTU QC *
1. Import OTU table
2. Calculate total number of sequences per OTU, per sample
3. Import 'NameSchematic.txt' to re-name sample names so they are more informative. Join with data.
4. Remove global singletons (OTUs that are found in only 1 sample with 1 sequence)
5. Plot supplemental figures showing total number of sequenes and the distribution of OTUs in each sample.
6. Run through 'pr2_rename_taxa' function to manually curate taxonomic group names (for summing sequences)
7. Remove unwanted samples due to low sequence number or high abundance of metazoa.

* Whole community plots *
1. Summarize the number of sequences in each sample by the manually designated "Taxa" column name
2. Pool replicate samples
3. Plot (ggplot2) relative abundance of each taxonomic group in the community
4. Plot OTU richness

* Composition of ciliates *
1. Aggregate data to major taxonomic group and Level 4 (approximately Class level).
2. Plot relative abundance of ciliate reads at the class level.

* Bubble plots - 3 groups *
1. Calculate relative abundance
2. Subset to three taxonomic groups of interest: Rhizaria, Ciliates, and Apicomplexa
3. Generate bubble plot

* Presence-absence UpsetR *
1. Aggregate the count of OTUs by habitat/sample type
2. Change to binary
3. Repeat with only ciliate data
4. Plot using UpSetR

* OTU richness - Ciliates only *
1. Subset ciliate reads from main data frame
2. If value does not equal 1, change to 1 (change to binary)
3. Aggregate by total number of OTUs in each ciliate class
4. Generate plot bubbles
5. Import distribution of ciliate OTUs (what samples were each OTUs found at)
6. Generate shaded grey area for each bubble plot

* Beta diversity metrics - MDS and ANOSIMS *
1. Transpose data and convert to numeric
2. Calculate relative abundance
3. Transform data for test for best fit: including 4th root, square root, and presence absence.
4. Calculate NMDS for each transformed data set
5. look at stress value to identify which transformation results in least stress
6. Import "meta_Vent.csv" and merge with data
7. This analysis uses 4th root transformed data for MDS plots.

* NMDS Figure *
1. Input points calculated from above section
2. Factor appropriate colors and shapes
3. Use ggplot2 to plot NMDS figure

* ANOSIM, SIMPER, & Alpha diversity *
ANOSIM
1. Import 4th root transformed data
2. Test various factors to run ANOSIM analyses - habitat, sediment depth, mat color
3. Save output results at .txt files.
4. Repeat, but remove control samples to run ANOSIM
SIMPER
1. Import 4th root transformed data
2. Run simper, again with various factors: habitat, mat color, & sediment horizon/depth
3. Save output as text file
Alpha diversity
1. Import R objects from before, need to use subsampled data.
2. Randomly sub-sample
3. Use 'diversity()' function on subsampled data to calculate Shannon and Inverse Simpson diversity metrics
4. Plot box plots to show distribution by sample type
Rarefaction curve
1. Use subsampled data
2. Factor sample types with desired colors to plot the rarefaction curve.
3. Use 'rarecurve()' function to generate rarefaction curves


## Contributors:
Alexis Pasulka & Sarah Hu - last updated October 2018
