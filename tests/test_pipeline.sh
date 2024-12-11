#!/bin/bash

# Create test data directories
mkdir -p data/test_data

# Generate test data
python scripts/generate_test_data.py

# Run pipeline steps
echo "Running diffusion network model..."
python src/drivers/driver_diffusion_network.py --config=config/test_config.json

echo 'Gather vulnerability covariates...'
Rscript src/drivers/driver_gather_vulnerability_covariates.R --config=config/test_config.json

echo "Running LGCM model..."
Rscript src/drivers/driver_LGCM.R --config=config/test_config.json

echo "Running diffusion network model with weights..."
python src/drivers/driver_diffusion_network.py --config=config/test_config.json --get_weights

echo 'Gather receptivity covariates...'
Rscript src/drivers/driver_gather_receptivity_covariates.R --config=config/test_config.json

echo "Running geospatial model..."
Rscript src/drivers/driver_geospatial_model_Re.R --config=config/test_config.json

# Verify outputs exist
echo "Verifying outputs..."
for file in $(cat <<EOF
data/test_data/case_data_withRe.csv
data/test_data/case_data_withRe_Covs.csv
data/test_data/case_data_withVuln_Covs.csv
data/test_data/imported_Pf_allYears.tif
data/test_data/imported_Pv_allYears.tif
EOF
); do
    if [ -f "$file" ]; then
        echo "✓ $file exists"
    else
        echo "✗ $file missing"
    fi
done