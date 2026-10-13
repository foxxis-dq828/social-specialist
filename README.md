# Social Specialist

MATLAB and R code for analysing social inference across three behavioural experiments. The analyses compare mentalizing, vote-based, diversity-weighted, and knowledge-weighted models using hierarchical Bayesian inference (HBI).

## Repository structure

```text
social-specialist/
├── README.md
├── code/
│   ├── exp1/
│   ├── exp2/
│   ├── exp3/
│   ├── model recovery/
│   └── posterior/
└── data/
    ├── experiment1.xlsx
    ├── experiment2.xlsx
    └── experiment3.xlsx
```

Each experiment directory contains:

- `BMI_run_exp*.mlx`: model fitting, HBI model comparison, predictions, and prior robustness analyses.
- `BMI_compare_exp*.mlx`: behavioural comparisons across strategy groups and robustness analyses using alternative classification thresholds.
- `src/models/`: model likelihood functions.
- `src/plotting/`: prediction and plotting functions.
- `src/r2_assignment/`: R²-based assignment analyses.
- `results/fitting/`: saved fitting results and logs.

## Models

| Model       | Description                                                   | Source file           |
| ----------- | ------------------------------------------------------------- | --------------------- |
| Mentalizing | Infers the source jar from agents’ judgments and knowledge    | `mentalizing_model.m` |
| VT          | Combines agents’ votes with a common vote weight              | `vt_model.m`          |
| DW          | Weights evidence according to knowledge diversity             | `dw_model.m`          |
| KW          | Weights evidence according to agents’ knowledge diagnosticity | `kw_model.m`          |
| Baseline    | Assigns equal probability to all six response categories      | `baseline_model.m`    |

The substantive models use beta distributions to generate probabilities over six response categories. Positive parameters are represented on the logarithmic scale during fitting. Model implementations and experimental conditions are defined separately for each experiment.

## Requirements

For the MATLAB analyses:

- MATLAB with Live Script support.
- The Computational and Behavioral Modeling (CBM) toolbox ([CBM toolbox](https://github.com/payampiray/cbm)).
- Optimization Toolbox.
- Statistics and Machine Learning Toolbox.
- Parallel Computing Toolbox for the parallel model recovery workflow.

For the R analysis:

- R.
- The `tidyverse` package.
- Quarto to render the `.qmd` document.

## Model recovery

`code/model recovery/model_recovery.mlx` simulates responses from the candidate models, refits the models, and summarises model recovery.

The current configuration uses 50 repetitions per generating model for each experiment and an eight-worker parallel pool. Update the experiment directories, CBM location, and result directory in the settings section before running it.

Saved recovery outputs are stored under `code/model recovery/result/`.

## Citation

Please cite the following paper when using the CBM toolbox:

Piray et al., "Hierarchical Bayesian inference for concurrent model fitting and comparison for group studies", *PLoS Computational Biology*, 2019.
