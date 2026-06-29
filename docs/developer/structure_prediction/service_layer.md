# Structural prediction service layer

The service layer provides a GUI-safe entry point for preparing structural prediction batches.

Main file:

    ppigfinder/services/structure_prediction_service.py

Main API:

    StructuralPredictionRequest
    prepare_structural_prediction_batch

This layer keeps Qt code thin. A future button or workflow should collect user input, build a StructuralPredictionRequest and call the service.

The service performs:

- multibackend job creation;
- token/resource planning;
- artifact generation;
- manifest writing and validation;
- backend input writing;
- Slurm script rendering.

It does not submit jobs.
