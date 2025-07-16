# Discovery Answers

## Q1: Will you test each foundation model environment on the same dataset for comparison?
**Answer:** Yes, let's actually use that

## Q2: Do you need to validate that all model downloads completed successfully before testing?
**Answer:** Definitely and each model should have its own environment setup, demo script and work with a .h5ad file

## Q3: Should the testing include both environment setup verification and actual embedding generation?
**Answer:** Yes, of course, keep in mind that I'm prefixing all environments to my scratch

## Q4: Will you run these tests on the BigPurple HPC cluster with CUDA support?
**Answer:** Yes, on v100 or a100 with cuda 12

## Q5: Do you want automated scripts that can test all environments sequentially without manual intervention?
**Answer:** Lets do that yes