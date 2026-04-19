# FES-MADM III Decision Studio

## User Guide and Operational Instructions

This directory contains the executable R-based analytical platform developed for implementation of the FES-MADM III framework. The platform provides a transparent, structured, and user-friendly environment for multicriteria evaluation, scenario comparison, and reproducible decision analytics.

The application supports both single-scenario execution and simultaneous multi-scenario processing, allowing users to evaluate alternatives under multiple external conditions within one integrated workflow.

---

## Main File

- `FESMADMIII_Decision_Studio.R`

This script launches the analytical environment and performs the full computational workflow of the framework.

---

## Core Functionalities

The platform supports:

- Import of structured Excel scenario datasets  
- Automatic preprocessing of input matrices  
- Execution of the FES-MADM III methodology  
- Entropy-based weighting procedures  
- Multicriteria aggregation and ranking  
- Scenario-by-scenario comparison  
- Consolidated export of outputs  
- Reproducible re-execution with updated data

---

## Input Modes

### 1. Single Scenario Execution

Users may upload one scenario dataset only (for example S0 Baseline) in order to perform a focused standalone evaluation.

Typical use cases include:

- Baseline ranking generation  
- Sensitivity testing of one environment  
- Preliminary model inspection  
- Isolated scenario analysis

### 2. Multi-Scenario Execution (Recommended)

The platform also supports simultaneous upload of all structured scenario files.

Example scenarios:

- S0 Baseline  
- S1 Inflation / Market Tightening  
- S2 Cyber Procurement Disruption  
- S3 Market Consolidation  
- S4 Budget Tightening

When all scenarios are imported together, the platform executes the full comparative analytical workflow and generates integrated outputs across all operating environments.

This mode is recommended for robust strategic assessment and comparative decision support.

---

## Recommended Workflow

### Step 1 — Launch the Application

Run:

`source("FESMADMIII_Decision_Studio.R")`

or open the file directly in RStudio and execute the script.

### Step 2 — Upload Scenario Files

Import one or more Excel scenario files from the `data/` directory.

All datasets follow a harmonized structure and can be processed directly without additional formatting.

### Step 3 — Automatic Processing

After import, the platform performs:

- Data recognition  
- Matrix validation  
- Normalization procedures  
- Entropy calculations  
- Weight estimation  
- Ranking synthesis  
- Comparative scenario calculations

### Step 4 — Results Review

Users may inspect generated:

- Rankings  
- Scenario comparisons  
- Aggregated outputs  
- Robustness signals  
- Structured result tables

### Step 5 — Export Outputs

Generated outputs may be exported for:

- Reporting  
- Manuscript support  
- Supplementary materials  
- Decision documentation  
- Archival storage

---

## Why Multi-Scenario Upload Matters

Uploading all scenarios simultaneously enables the platform to evaluate how rankings and priorities evolve under changing environments.

This allows users to identify:

- Stable high-performing alternatives  
- Scenario-sensitive alternatives  
- Robust decision options  
- Vulnerability to external shocks  
- Consistency of rankings across regimes

Simultaneous scenario execution therefore provides substantially richer decision intelligence than isolated single-case analysis.

---

## Input File Requirements

All uploaded Excel files should preserve the provided repository structure and formatting.

Users are encouraged to use the supplied scenario datasets directly or modify copies while maintaining the same schema.

---

## Reproducibility Note

Because all scenarios are standardized, repeated execution of the same files should reproduce identical outputs under the same computational environment.

This supports transparency, auditability, and independent verification.

---

## Technical Environment

Recommended:

- Recent version of R  
- RStudio  
- Required packages referenced inside the script

---

## Supportive Research Use

The platform may be used for:

- Procurement analytics  
- Supplier evaluation  
- Strategic sourcing  
- Comparative benchmarking  
- Resilience-oriented ranking  
- Multicriteria scenario planning

---

## Final Remark

The Decision Studio was developed not merely as a calculator, but as an integrated analytical environment enabling structured and auditable scenario-based decision support.
