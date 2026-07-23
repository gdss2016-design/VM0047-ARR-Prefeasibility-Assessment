# VM0047-ARR-Prefeasibility-Assessment

Google Earth Engine workflow for VM0047 ARR prefeasibility assessment, including eligibility screening, historical forest-loss analysis, donor/control matching, balance diagnostics, vegetation-condition time series analysis, and project candidate prioritisation.

This repository contains a modular Google Earth Engine (GEE) workflow for VM0047-style ARR (Afforestation, Reforestation and Revegetation) prefeasibility assessment. The workflow identifies and evaluates potential project areas through land eligibility screening, historical forest-loss assessment, true eligible-area calculation, environmental covariate extraction, donor/control matching, balance diagnostics, and annual vegetation-condition analysis.

The workflow is designed to support early-stage project screening and prioritisation before detailed field surveys, baseline development, and carbon accounting. Outputs include eligible candidate locations, donor-control comparisons, balance metrics, annual vegetation-condition indicators (NDVI, EVI, NDMI and NBR), and true pixel-based eligible area estimates.

The workflow consists of six stages:

- Script 1 – Eligibility Screening and True Area Calculation
- Script 2 – Environmental Covariate Extraction
- Script 3 – Project-Control Matching
- Script 4 – Balance Diagnostics
- Script 5 – Annual Vegetation Condition Time Series
- Script 6 – Final Candidate Location Summary

This repository supports prefeasibility assessment only and does not calculate carbon stocks, carbon removals, carbon credits, or verified carbon units (VCUs). Additional field validation and methodology-specific carbon accounting are required before project implementation.
