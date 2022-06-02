

## HIP
HIP is a python package that predicts the Human-Infection Potential from the S protein sequences of  alpha and beta coronaviruses as described in the following paper:

*Gonzalez-Isunza, G., Jawaid M., Cox, D., Vazquez, M., Arsuaga, J. 2022. A machine learning method to predict the human infection potential of animal coronaviruses. Manuscript in preparation.*

Results presented at Gonzalez-Isunza et al. are available at: 
+ Alpha and beta coronaviruses:   
 ./hip_scores/alpha_beta_scores.csv
+ Alpha and beta coronaviruses excluding SARS2:   
 ./hip_scores/ab_no_sars2_scores.csv  

#### Basic usage
**Setup**
+ Clone or download ZIP the repository
+ Setup up the environment by navigating in the terminal to /HIP and then typing:   
```
pip install .
```  
**Reproducing results from Gonzalez-Isunza et al.**   
+ After navigating in the terminal to /HIP type:    
```
python3 hip_reproduce.py
``` 
+ To train any other dataset, create a config file using the template provided at the data folder and add the file path at the end.    
  ```
  python3 hip_reproduce.py ./data/my_own_data_config.yml
  ```
**Predicting HIP scores for spike amino acid sequences**  
By default, HIP will compute the score from the alpha_beta model (full dataset).    
+ If no fasta file is provided, it will compute de score for SARS-CoV-2.
After navigating in the terminal to /HIP type:
```
python3 hip_predict.py
```  
+ From a fasta file:
```
python3 hip_predict.py path_to_fasta_file
```
+ To use a different model for prediction, add the model name at the end 
```
python3 hip_predict.py path_to_fasta_file model_name
```

#### Available datasets
+ Alpha and beta coronaviruses (full dataset)
    + 2,534 unique spike protein sequences annotated for binding condition to human receptor
    + alpha_beta.csv  
+ Alpha and beta coronaviruses excluding SARS2 viruses
    + After removing all SARS2 viruses, this dataset is identical to alpha_beta.csv
    + ab_no_sars2.csv   
 




