/***************************************************************
FINAL SCRIPT 2A:
MEMORY-SAFE COVARIATE EXTRACTION FOR MATCHING

Run twice:
1. PROCESS_PROJECT = true
** PROCESS_PROJECT = false

Inputs:** 01_Project_EligibleScored_1000m
- 01_Donor_EligibleScored_1000m

Outputs:
- 02A_Project_**variates_1000m
- 02A_Donor**ovariates_1000m


****************************************************************/


/***************************************************************
0. CONFIGURATION
****************************************************************/

var ASSET_ROOT = 'projects/ee-gdss2016/assets/Pre-feasibility_assessment/';


// Run first with true, then run again with false.
var PROCESS_PROJECT = false;

var PROJECT_INPUT_ASSET =
  ASSET_ROOT + '01_Project_EligibleScored_1000m_';

var DONOR_INPUT_ASSET =
  ASSET_ROOT + '01_Donor_EligibleScored_1000m_';

var GRID_SIZE_M = 1000;
var CRS = 'EPSG:32644';

var PRE_START_YEAR = 2021;
var PRE_END_YEAR = 2025;

var SEASON_START_MONTH = 10;
var SEASON_START_DAY = 1;
var SEASON_END_MONTH = 12;
var SEASON_END_DAY = 31;

// Memory-safe limits.
var MAX_PROJECT_POINTS = 300;
var MAX_DONOR_POINTS = 1500;


/***************************************************************
1. SELECT INPUT AND OUTPUT
****************************************************************/

var inputAsset;
var outputAsset;
var outputDescription;
var processingType;
var maxPoints;

if (PROCESS_PROJECT) {
  inputAsset = PROJECT_INPUT_ASSET;
  outputAsset = ASSET_ROOT + '02A_Project_Covariates_1000m_';
  outputDescription = '02A_Project_Covariates_1000m_';
  processingType = 'project';
  maxPoints = MAX_PROJECT_POINTS;
} else {
  inputAsset = DONOR_INPUT_ASSET;
  outputAsset = ASSET_ROOT + '02A_Donor_Covariates_1000m_';
  outputDescription = '02A_Donor_Covariates_1000m_';
  processingType = 'donor';
  maxPoints = MAX_DONOR_POINTS;
}

print('Running Script 2A for', processingType);
print('Input asset', inputAsset);
print('Output asset', outputAsset);
print('Max points used', maxPoints);


/***************************************************************
2. LOAD AND LIMIT SELECTED CELLS
****************************************************************/

var allCells = ee.FeatureCollection(inputAsset);

// selected_for_matching is numeric 1/0 from Script 1.
var selectedCells = allCells
  .filter(ee.Filter.gt('selected_for_matching', 0))
  .filter(ee.Filter.gt('eligible_area_ha', 0));

// Sort by non-carbon feasibility score and eligible area.
// Limit BEFORE Landsat sampling to avoid memory error.
var points = selectedCells
  .sort('feasibility_score_non_carbon', false)
  .limit(maxPoints);

print('Selected limited cells sample', points.limit(10));

// Do not print selectedCells.size() for very large donor asset.
// Do not Map.addLayer points.


/***************************************************************
3. ANALYSIS GEOMETRY
****************************************************************/

// Use limited points geometry only.
// This keeps Landsat filtering spatially small.
var analysisGeom = points.geometry()
  .buffer(5000)
  .bounds(1000);


/***************************************************************
4. PROJECTION
****************************************************************/

var PROJ_30M = ee.Projection(CRS).atScale(30);


/***************************************************************
5. LANDSAT CLOUD MASK AND SCALE
****************************************************************/

function maskAndScaleLandsat(img) {
  var qa = img.select('QA_PIXEL');

  var cloud = qa.bitwiseAnd(1 << 3).eq(0);
  var shadow = qa.bitwiseAnd(1 << 4).eq(0);
  var snow = qa.bitwiseAnd(1 << 5).eq(0);

  var mask = cloud.and(shadow).and(snow);

  var optical = img.select([
      'SR_B2',
      'SR_B3',
      'SR_B4',
      'SR_B5',
      'SR_B6',
      'SR_B7'
    ])
    .multiply(0.0000275)
    .add(-0.2)
    .rename([
      'BLUE',
      'GREEN',
      'RED',
      'NIR',
      'SWIR1',
      'SWIR2'
    ]);

  return optical
    .updateMask(mask)
    .copyProperties(img, ['system:time_start']);
}


/***************************************************************
6. SEASONAL LANDSAT COMPOSITE
****************************************************************/

function seasonalComposite(year) {
  year = ee.Number(year);

  var start = ee.Date.fromYMD(
    year,
    SEASON_START_MONTH,
    SEASON_START_DAY
  );

  var end = ee.Date.fromYMD(
    year,
    SEASON_END_MONTH,
    SEASON_END_DAY
  );

  var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate(start, end)
    .filterBounds(analysisGeom)
    .map(maskAndScaleLandsat);

  var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
    .filterDate(start, end)
    .filterBounds(analysisGeom)
    .map(maskAndScaleLandsat);

  var composite = l8.merge(l9).median();

  var ndvi = composite.normalizedDifference(['NIR', 'RED'])
    .rename('NDVI');

  var evi = composite.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
    {
      NIR: composite.select('NIR'),
      RED: composite.select('RED'),
      BLUE: composite.select('BLUE')
    }
  ).rename('EVI');

  var ndmi = composite.normalizedDifference(['NIR', 'SWIR1'])
    .rename('NDMI');

  var nbr = composite.normalizedDifference(['NIR', 'SWIR2'])
    .rename('NBR');

  return ee.Image.cat([
    ndvi,
    evi,
    ndmi,
    nbr
  ])
  .setDefaultProjection(PROJ_30M)
  .set('year', year);
}


/***************************************************************
7. BUILD PRE-PROJECT SI IMAGE
****************************************************************/

var years = ee.List.sequence(PRE_START_YEAR, PRE_END_YEAR);

var annualSI = ee.ImageCollection(
  years.map(function(y) {
    return seasonalComposite(y);
  })
);

var preMeanSI = annualSI.mean()
  .setDefaultProjection(PROJ_30M)
  .rename([
    'pre_NDVI_mean',
    'pre_EVI_mean',
    'pre_NDMI_mean',
    'pre_NBR_mean'
  ]);

var preSdSI = annualSI.reduce(ee.Reducer.stdDev())
  .setDefaultProjection(PROJ_30M)
  .rename([
    'pre_NDVI_sd',
    'pre_EVI_sd',
    'pre_NDMI_sd',
    'pre_NBR_sd'
  ]);


/***************************************************************
8. TERRAIN COVARIATES
****************************************************************/

var terrain = ee.Algorithms.Terrain(
  ee.Image('USGS/SRTMGL1_003')
);

var elevation = terrain.select('elevation')
  .setDefaultProjection(PROJ_30M)
  .rename('elevation');

var slope = terrain.select('slope')
  .setDefaultProjection(PROJ_30M)
  .rename('slope');


/***************************************************************
9. ESA CCI AGB REFERENCE

Reference covariate only.
Not carbon sequestration.
****************************************************************/

var agbReference = ee.ImageCollection('ESA/CCI/Above_Ground_Biomass/V6_0')
  .filterDate(PRE_START_YEAR + '-01-01', (PRE_END_YEAR + 1) + '-01-01')
  .filterBounds(analysisGeom)
  .mean()
  .select(0)
  .setDefaultProjection(PROJ_30M)
  .rename('pre_AGB_reference');


/***************************************************************
10. COMBINE COVARIATE IMAGE
****************************************************************/

var covariateImage = ee.Image.cat([
  preMeanSI,
  preSdSI,
  elevation,
  slope,
  agbReference
]).unmask(-9999);


/***************************************************************
11. PROPERTIES TO KEEP FROM SCRIPT 1
****************************************************************/

var KEEP_PROPS = [
  'cell_id',
  'grid_source',
  'grid_size_m',
  'nominal_area_ha',
  'x_index',
  'y_index',

  'processing_type',
  'project_start_year',

  'recent_loss_check_done_in_script',
  'recent_loss_period',
  'recent_lossyear_min_code',
  'recent_lossyear_max_code',

  'min_priority_eligible_area_ha',
  'min_priority_candidate_fraction',

  'true_cell_area_inside_boundary_ha',

  'broad_candidate_area_ha',
  'priority_candidate_area_ha',
  'cropland_candidate_area_ha',
  'hard_excluded_area_ha',

  'wc_broad_candidate_fraction',
  'wc_priority_candidate_fraction',
  'wc_cropland_candidate_fraction',
  'hard_excluded_fraction',

  'pass_min_priority_area',
  'pass_min_priority_fraction',
  'pass_forest2000_30pct',
  'pass_recent_loss_threshold',

  'wc_tree_fraction',
  'wc_built_fraction',
  'wc_water_fraction',
  'wc_wetland_fraction',
  'wc_mangrove_fraction',

  'forest2000_10pct_fraction',
  'forest2000_30pct_fraction',

  'recent_forest_loss_fraction',
  'recent_forest_loss_pct',
  'recent_forest_loss_area_ha',

  'strict_eligible_area_ha',
  'priority_eligible_area_ha',
  'broad_eligible_area_ha',
  'eligible_area_ha',

  'strict_eligibility_pass',
  'priority_eligibility_pass',
  'selected_for_matching',

  'feasibility_score_non_carbon',
  'area_method',
  'note_script1'
];


/***************************************************************
12. SAMPLE COVARIATES
****************************************************************/

var sampled = covariateImage.sampleRegions({
  collection: points,
  properties: KEEP_PROPS,
  scale: 30,
  projection: PROJ_30M,
  geometries: true,
  tileScale: 16
});


/***************************************************************
13. ADD READY FLAGS
****************************************************************/

var withFlags = sampled.map(function(f) {
  var hasValidSI = ee.Number(f.get('pre_NDVI_mean')).neq(-9999)
    .and(ee.Number(f.get('pre_EVI_mean')).neq(-9999))
    .and(ee.Number(f.get('pre_NDMI_mean')).neq(-9999))
    .and(ee.Number(f.get('pre_NBR_mean')).neq(-9999));

  var ready = hasValidSI
    .and(ee.Number(f.get('selected_for_matching')).gt(0))
    .and(ee.Number(f.get('eligible_area_ha')).gt(0));

  return f.set({
    script2a_processing_type: processingType,
    has_valid_pre_project_si: ee.Number(
      ee.Algorithms.If(hasValidSI, 1, 0)
    ),
    ready_for_matching: ee.Number(
      ee.Algorithms.If(ready, 1, 0)
    ),
    covariate_years: PRE_START_YEAR + '-' + PRE_END_YEAR,
    seasonal_window: SEASON_START_MONTH + '-' + SEASON_END_MONTH,
    max_points_used_script2a: maxPoints,
    note_script2a: 'Memory-safe covariate extraction only.'
  });
});


/***************************************************************
14. DIAGNOSTICS
****************************************************************/

var readyRows = withFlags.filter(
  ee.Filter.gt('ready_for_matching', 0)
);

print('Script 2A output sample', withFlags.limit(10));
print('Script 2A ready rows', readyRows.size());

print(
  'Script 2A ready eligible area ha',
  readyRows.aggregate_sum('eligible_area_ha')
);


/***************************************************************
15. EXPORT
****************************************************************/

Export.table.toAsset({
  collection: withFlags,
  description: outputDescription,
  assetId: outputAsset
});

Export.table.toDrive({
  collection: withFlags,
  description: outputDescription + '_CSV',
  fileNamePrefix: outputDescription,
  fileFormat: 'CSV'
});
