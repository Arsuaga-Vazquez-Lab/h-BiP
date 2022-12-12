

## h-BiP
h-BiP is a python package that predicts the Human-Infection Potential from the S protein sequences of  alpha and beta coronaviruses as described in the following paper:

*Gonzalez-Isunza, G., Jawaid M., Cox, D., Vazquez, M., Arsuaga, J. 2022. A machine learning method to predict the human infection potential of animal coronaviruses. Manuscript in preparation.*

Results presented at Gonzalez-Isunza et al. are available at: 
+ Alpha and beta coronaviruses:   
 ./hip_scores/alpha_beta_scores.csv
+ Alpha and beta coronaviruses excluding SARS2:   
 ./hip_scores/ab_no_sars2_scores.csv  

### Basic usage
**Setup**
+ Clone or download ZIP the repository.
+ Setup up the environment by navigating in the terminal to /HIP and then typing:   
```
myterminal$ pip install .
```  
**Reproducing results from Gonzalez-Isunza et al.**   
+ After navigating in the terminal to /HIP type:    
```
myterminal$ python3 hip_reproduce.py
``` 
+ To train any other dataset, create a config file using the template provided at the data folder and add the file path at the end.    
  ```
  myterminal$ python3 hip_reproduce.py ./data/my_own_data_config.yml
  ```
**Predicting h-BiP scores for spike amino acid sequences**  
By default, h-BiP will compute the score from the alpha_beta model (full dataset).    
+ If no fasta file is provided, it will compute de score for SARS-CoV-2.
After navigating in the terminal to /HIP type:
```
myterminal$ python3 hbip_predict.py
```  
+ From a fasta file:
```
myterminal$ python3 hbip_predict.py path_to_fasta_file
```
+ To use a different model for prediction, add the model name at the end. 
```
myterminal$ python3 hbip_predict.py path_to_fasta_file model_name
```

#### Available datasets
+ Alpha and beta coronaviruses (full dataset):
    + 2,534 unique spike protein sequences annotated for binding condition to human receptor.
    + alpha_beta.csv  
+ Alpha and beta coronaviruses excluding SARS2 viruses:
    + After removing all SARS2 viruses, this dataset is identical to alpha_beta.csv.
    + ab_no_sars2.csv   
 
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




