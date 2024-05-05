

## h-BiPs
This is a modified version of h-BiP that will allow prediction on short reads. This new functionality is underconstruction. Stay tunned!

h-BiP is a python package that predicts the human Binding Potential from the S protein sequences of  alpha and beta coronaviruses as described in the following paper:

*Gonzalez-Isunza, G., Jawaid, M.Z., Liu, P., Cox, D.L., Vazquez, M. and Arsuaga, J., 2023. Using machine learning to detect coronaviruses potentially infectious to humans. Scientific Reports, 13(1), p.9319.*

Results presented at Gonzalez-Isunza et al. are available at: 
+ Alpha and beta coronaviruses:   
 ./hbip_scores/alpha_beta_scores.csv
+ Alpha and beta coronaviruses excluding SARS2:   
 ./hbip_scores/alpha_beta_before_2020_scores.csv  
 
Any publication that discloses findings arising from using this data or the source code should [cite Gonzalez-Isunza et al.](#citing-this-work)

### Basic usage
**Setup**
+ Clone or download ZIP the repository.
+ Create and activate an environment using `conda` or a virtual environment with Pyhon=3.9.7.
+ CD into the repository directory.
+ Pip install the hbip package
```
myterminal ~/h-BiP $  pip install .
```

**Reproducing results from Gonzalez-Isunza et al.**   
+ Open the terminal, navigate to /h-BiP and type:    
```
myterminal ~/h-BiP $ python3 hbip_reproduce.py
``` 
+ To train any other dataset, create a config file using the template provided at the data folder and add the file path at the end.    

```
myterminal ~/h-BiP $ python3 hbip_reproduce.py ./data/my_own_data_config.yml
```

**Predicting h-BiP scores for spike amino acid sequences**  
By default, h-BiP will compute the score from the alpha_beta model (full dataset).    
+ If no fasta file is provided, it will compute de score for SARS-CoV-2.
After navigating in the terminal to /h-BiP type:
```
myterminal ~/h-BiP $  python3 hbip_predict.py
```
+ From a fasta file:
```
myterminal ~/h-BiP $  python3 hbip_predict.py path_to_fasta_file
```
+ To use a different model for prediction, add the model name at the end. 
```
myterminal ~/h-BiP $  python3 hbip_predict.py path_to_fasta_file model_name
```

### Data

#### Available datasets
+ Alpha and beta coronaviruses (full dataset):
    + 2,534 unique spike protein sequences annotated for binding condition to human receptor.
    + alpha_beta_hip_train.csv and alpha_beta_hip_test.csv  
+ Alpha and beta coronaviruses excluding SARS2 viruses and all viruses uploaded to NCBI after Dec. 31, 2019:
    + alpha_beta_before_2020_hip_train.csv and alpha_beta_before_2020_hip_test.csv   
 
#### Dataset structure
Fields marked with asterix are used in the code (position is irrelevant).
+ *Accession: string   
Identifier for the sequence.
+ *Sequence: string   
Amino acid sequence for the spike protein (S).
+ Species: string   
Virus' species.
+ *Species_agg: string   
Simplified (aggregated) label for Species. Current labels are:   
    + hCoV-OC43, hCoV-HKU1, Beta 1, MERS, MERS-related, other, SARS-CoV-2, SARS-CoV-1, Other Sarbecovirus, Beta other, PorcineEp, hCoV-NL63, hCoV-229E, Other Alpha.       
    + The only label needed in the code is 'SARS-CoV-2'   
+ *Virus: string   
Virus name.
+ Host: string   
Host's name.
+ Host_agg: string   
Simplified (aggregated) label for Host.
+ Country: string   
Virus' origin.
+ Human: integer (0/1)  
Human host will have a value of 1 and 0 otherwise.
+ Binds: integer (0/1)  
Experimental evidence of binding to human receptor.

### Citing this work
If you use the code or data in this package, please cite:
```
@article {Gonzalez-Isunza2023using,
	author = {Gonzalez-Isunza, Georgina and Jawaid, M. Zaki and Liu, Pengyu and Cox, Daniel L. and Vazquez, Mariel and Arsuaga, Javier},
	title = {Using machine learning to detect coronaviruses potentially infectious to humans},
	journal = {Scientific reports},
	year = {2023},
	publisher = {Nature Publishing Group UK London}
}
```


