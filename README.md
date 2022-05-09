

## HIP
HIP is a python package that predicts the Human-Infection Potential from the S protein sequences of  alpha and beta coronaviruses as described in the following paper:

*Gonzalez-Isunza, G., Jawaid M., Cox, D., Vazquez, M., Arsuaga, J. 2022. A machine learning method to predict the human infection potential of animal coronaviruses. Manuscript in preparation.*

Results presented at Gonzalez-Isunza et al. are available at: 
+ Alpha and beta coronaviruses:   
 ./outputs/alpha_beta_scores.csv
+ Alpha and beta coronaviruses excluding SARS2:   
 ./outputs/ab_no_sars2_scores.csv  

#### Basic usage
1. Clone or download ZIP the repository
2. Setup up the environment by navigating in the terminal to /HIP and then typing:   
```pip install .```
3. To reproduce results for all alpha and beta coronaviruses go to the terminal and type:    
```python3 main.py``` 
4. To train any other dataset, create a config file using the template provided at the data folder and add the file path at the end.    
  ```python3 main.py ./data/my_own_data_config.yml``` 

#### Available datasets
+ alpha_beta.csv   
Alpha and beta coronaviruses annotated for binding condition to human receptor. It consists of
2,534 unique spike protein sequences. 
+ ab_no_sars2.csv   
After removing all SARS2 viruses, this dataset is identical to alpha_beta.csv. 




