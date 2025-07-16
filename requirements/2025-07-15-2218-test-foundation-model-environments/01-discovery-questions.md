# Discovery Questions

## Q1: Will you test each foundation model environment on the same dataset for comparison?
**Default if unknown:** Yes (enables meaningful comparison between models)

## Q2: Do you need to validate that all model downloads completed successfully before testing?
**Default if unknown:** Yes (prevents runtime failures due to missing models)

## Q3: Should the testing include both environment setup verification and actual embedding generation?
**Default if unknown:** Yes (comprehensive testing ensures full functionality)

## Q4: Will you run these tests on the BigPurple HPC cluster with CUDA support?
**Default if unknown:** Yes (matches the project's target deployment environment)

## Q5: Do you want automated scripts that can test all environments sequentially without manual intervention?
**Default if unknown:** Yes (enables efficient batch testing and reduces manual errors)