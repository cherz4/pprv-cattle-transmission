# pprv-cattle-transmission
This repository accompanies the following article and it provides code for reproducing the analyses and modeling as well as associated figures.

Empirical and model-based evidence for a negligible role of cattle in peste des petits ruminants virus transmission and eradication

Catherine M. Herzog*, Fasil Aklilu*, Demeke Sibhatu, Dereje Shegu, Redeat Belaineh, Abde Aliy Mohammed, Menbere Kidane,  Claudia Schulz,
Brian J. Willett, Sarah Cleaveland, Dalan Bailey, Andrew R. Peters, Isabella M. Cattadori, Peter J. Hudson, Hagos Asgedom, Joram Buza, 
Mesfin Sahle Forza, Tesfaye Rufael Chibssa, Solomon Gebredufe, Nick Juleff, Ottar N. Bjørnstad, Michael D. Baron, Vivek Kapur (2024). Commmunications Biology.

\* These authors contributed equally


# License and citation information
If you use the code or data provided here, please make sure to cite the above manuscript.

# Directories
- `data`: all clinical, molecular, tissue culture, and serological data
- `scripts`: all code

# How to reproduce the analyses
A guide to reproducing the analysis from the paper follows.

## Getting the code
First download the code. The recommended way is to git clone our Github repository from the command line:

`git clone https://github.com/cherz4/pprv-cattle-transmission.git`

 or by downloading it manually via Github's download button (see green Code button on main page of repository for options).

## Running the code
All code provided are as R scripts. At the top of each script the packages (libraries) needed will be indicated and you will want to install these ahead of time with
`install.packages(package_name)`

All scripts should run as long as the relative paths are not disrupted. Scripts corresponding to clinical, molecular, and serological data analysis and plotting for all five experimental trials are present, as well as calculation of the transmission probabilities and credible intervals that were used in the comparatmental modeling next-generation matrix framework for calculating R0. Pay careful attention to commenting in the code, some blocks may need to be commented in/out (with #) and run again to reproduce all analyses or figures.


# Data Dictionary

## Clinical Data
Scoring was based on the scoring system used in the following publication:
Pope, R. A. et al. Early events following experimental infection with Peste-Des-Petits ruminants virus suggest immune cell targeting. PLoS One 8, e55830 (2013).
In general, a higher score indicates a more noticeable clinical sign.

Column name and definition
- eartag: unique animal ID on eartag
- date: date sample was taken
- dpi: day post innoculation
- barn: which part this animal was located in during the experiment, among 6 barns
- innoc: binary variable, 1 if animal was innoculated, 0 if not innoculated
- barni: combination of barn number and innoculation status (a = inoculated, b = not inoculated)
- general: general behavior score based on Pope et al 2013
- temp: temperature score based on Pope et al 2013
- feces: based on Pope et al 2013
- odnd: ocular and nasal discharge score based on Pope et al 2013
- headmucos: mucosal membrane score based on Pope et al 2013
- resp: based on Pope et al 2013
- totscore: sum of the six scoring columns, calculated by hand by experimenters
- autosum: sum of the six scoring columns, calculated by Excel formula
- tempc: rectal temperature reading in Celsius

## Molecular and Cell Culture Data
Not every study animal had every sample type taken or tested with each type of test in each experiment. Isolation in DMEM was only carried out in trial 1. See published paper and supplementary information for more details.

Column name and definition
- eartag: unique animal ID on eartag
- species: species type of animal with this ID (bovine/cattle, caprine/goat, ovine/sheep) 
- date: date sample was taken
- dpi: day post innoculation
- barn: which part this animal was located in during the experiment, among 6 barns
- innoc: binary variable, 1 if animal was innoculated, 0 if not innoculated
- barni: combination of barn number and innoculation status (a = inoculated, b = not inoculated)
- od: IDVET AgELISA test optical density reading on ocular swabs
- sp: IDVET AgELISA sample to positive control ratio (s/p %) calculation on ocular swabs
- status: IDVET AgELISA test result (Neg = Negative, Pos = Positive) on ocular swabs
- od_f: IDVET AgELISA test optical density reading on rectal (f = fecal) swabs
- sp_f: IDVET AgELISA sample to positive control ratio (s/p %) calculation on rectal (f = fecal) swabs
- status_f: IDVET AgELISA test result (Neg = Negative, Pos = Positive) on rectal (f = fecal) swabs
- notes1: Any relevant notes about animal or test outcomes
- notes2: Any additional relevant notes about animal or test outcomes
- expt: Experiment type: Experimental Barn, Positive Control, Negative Control, or left blank
- o1: Cycle threshold (ct) value from rQT-PCR measurement of ocular swab, technical replicate 1, U = tested but result undetermined, blank = not tested
- o2: Cycle threshold (ct) value from rQT-PCR measurement of ocular swab, technical replicate 2, U = tested but result undetermined, blank = not tested
- n1: Cycle threshold (ct) value from rQT-PCR measurement of nasal swab, technical replicate 1, U = tested but result undetermined, blank = not tested
- n2: Cycle threshold (ct) value from rQT-PCR measurement of nasal swab, technical replicate 2, U = tested but result undetermined, blank = not tested
- r1: Cycle threshold (ct) value from rQT-PCR measurement of rectal swab, technical replicate 1, U = tested but result undetermined, blank = not tested
- r2: Cycle threshold (ct) value from rQT-PCR measurement of rectal swab, technical replicate 2, U = tested but result undetermined, blank = not tested
- wb1: Cycle threshold (ct) value from rQT-PCR measurement of whole blood sample, technical replicate 1, U = tested but result undetermined, blank = not tested
- wb2: Cycle threshold (ct) value from rQT-PCR measurement of whole blood sample swab, technical replicate 2, U = tested but result undetermined, blank = not tested
- omean: Mean of o1 and o2 ct values
- nmean: Mean of n1 and n2 ct values
- rmean: Mean of r1 and r2 ct values
- wbmean: Mean of wb1 and wb2 ct values
- retest:  (often blank)
- isolation1_dbps: binary variable, 1 if any tissue culture score indicated successful virus isolation from any swab type in the next few columns, in DBPS media
- i1o: binary variable, 1 if any tissue culture score indicated successful virus isolation from an ocular swab, 0 if no isolation from an ocular swab, in DBPS media
- i1n: binary variable, 1 if any tissue culture score indicated successful virus isolation from an ocular swab, 0 if no isolation from an nasal swab, in DBPS media
- i1r: binary variable, 1 if any tissue culture score indicated successful virus isolation from an ocular swab, 0 if no isolation from an rectal swab, in DBPS media
- isolation2_dmem: binary variable, 1 if any tissue culture score indicated successful virus isolation from any swab type in the next few columns, in DMEM media
- i2o: binary variable, 1 if any tissue culture score indicated successful virus isolation from an ocular swab, 0 if no isolation from an ocular swab, in DMEM media
- i2n: binary variable, 1 if any tissue culture score indicated successful virus isolation from an ocular swab, 0 if no isolation from an nasal swab, in DMEM media
- i2r: binary variable, 1 if any tissue culture score indicated successful virus isolation from an ocular swab, 0 if no isolation from an rectal swab, in DMEM media

## Serology Data

Column name and definition
- eartag: unique animal ID on eartag
- species: species type of animal with this ID (bovine/cattle, caprine/goat, ovine/sheep) 
- barn: which part this animal was located in during the experiment, among 6 barns
- dpi: day post innoculation
- od1: IDVET cELISA test optical density reading on serum sample, technical replicate 1
- od2: IDVET cELISA test optical density reading on serum sample, technical replicate 2
- od_avg: IDVET cELISA test optical density reading on serum sample, average of technical replicates
- sn: IDVET cELISA measure of serum neutralization (SN)
- status: IDVET cELISA test result (Neg = Negative, Pos = Positive) on serum samples based on SN.

