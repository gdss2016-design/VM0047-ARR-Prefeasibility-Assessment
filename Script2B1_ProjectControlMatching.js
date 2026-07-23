/***************************************************************
FINAL SCRIPT 3:
MEMORY-SAFE MATCHING ONLY

Inputs:
- 02_Project_Covariates_1000m
- 02_Donor_Covariates_1000m

Output:
- 03_Matched_Project_Control_Pairs_1000m

Purpose:
- Loads Script 2 covariate exports.
- Filters ready_for_matching = 1.
- Uses a smaller project and donor subset.
- Matches each project cell to one best donor cell.
***********************************/


/***************************************************************
0. CONFIGURATION
****************************************************************/

var ASSET_ROOT = 'projects/ee-gdss2016/assets/Pre-feasibility_assessment/';

var PROJECT_COVARIATE_ASSET =
  ASSET_ROOT + '02_Project_Covariates_1000m_';

var DONOR_COVARIATE_ASSET =
  ASSET_ROOT + '02_Donor_Covariates_1000m_';

var OUTPUT_ASSET =
  ASSET_ROOT + '03_Matched_Project_Control_Pairs_1000m_';

var GRID_SIZE_M = 1000;


var MAX_PROJECT_POINTS = 20;
var MAX_DONOR_POINTS = 300;

var SCENARIO_NAME =
  'S1_ConservativeARR_MatchingOnly_1000m';

// Reduced covariates for memory-safe first matching.
// Keep only strongest matching variables first.
var COVARIATES = [
  'pre_NDVI_mean',
  'pre_EVI_mean',
  'pre_NDMI_mean',
  'pre_NBR_mean',
  'elevation',
  'slope',
  'wc_priority_candidate_fraction',
  'forest2000_30pct_fraction',
  'recent_forest_loss_fraction'
];


/***************************************************************
1. LOAD INPUTS
****************************************************************/

var projectAll = ee.FeatureCollection(PROJECT_COVARIATE_ASSET);
var donorAll = ee.FeatureCollection(DONOR_COVARIATE_ASSET);

print('Running Script 2B1 memory-safe matching only');
print('Project covariate asset', PROJECT_COVARIATE_ASSET);
print('Donor covariate asset', DONOR_COVARIATE_ASSET);
print('Output asset', OUTPUT_ASSET);
print('MAX_PROJECT_POINTS', MAX_PROJECT_POINTS);
print('MAX_DONOR_POINTS', MAX_DONOR_POINTS);


/***************************************************************
2. FILTER VALID ROWS
****************************************************************/

function filterValidCovariates(fc) {
  var out = fc
    .filter(ee.Filter.gt('ready_for_matching', 0))
    .filter(ee.Filter.gt('selected_for_matching', 0))
    .filter(ee.Filter.gt('eligible_area_ha', 0));

  COVARIATES.forEach(function(c) {
    out = out
      .filter(ee.Filter.notNull([c]))
      .filter(ee.Filter.neq(c, -9999));
  });

  return out;
}

// Keep best-scored cells only.
var projectReady = filterValidCovariates(projectAll)
  .sort('feasibility_score_non_carbon', false)
  .limit(MAX_PROJECT_POINTS);

var donorReady = filterValidCovariates(donorAll)
  .sort('feasibility_score_non_carbon', false)
  .limit(MAX_DONOR_POINTS);

// Keep diagnostics light.
print('Project ready count', projectReady.size());
print('Donor ready count', donorReady.size());


/***************************************************************
3. SAFE NUMBER HELPERS
****************************************************************/

function safeNumber(value, fallback) {
  return ee.Number(
    ee.Algorithms.If(
      value,
      value,
      fallback
    )
  );
}

function safeSd(value) {
  return ee.Number(
    ee.Algorithms.If(
      value,
      ee.Algorithms.If(
        ee.Number(value).gt(0),
        value,
        1
      ),
      1
    )
  );
}

function donorStats(fc, propertyName) {
  var meanRaw = fc.aggregate_mean(propertyName);
  var sdRaw = fc.aggregate_total_sd(propertyName);

  return ee.Dictionary({
    mean: safeNumber(meanRaw, 0),
    sd: safeSd(sdRaw)
  });
}


/***************************************************************
4. DONOR STANDARDIZATION STATS
****************************************************************/

var statsDict = ee.Dictionary({});

COVARIATES.forEach(function(c) {
  statsDict = statsDict.set(c, donorStats(donorReady, c));
});


/***************************************************************
5. MATCH ONE PROJECT CELL TO BEST DONOR CELL
****************************************************************/

function matchOneProject(projectFeature) {
  projectFeature = ee.Feature(projectFeature);

  var scoredDonors = donorReady.map(function(donorFeature) {
    donorFeature = ee.Feature(donorFeature);

    var d2 = ee.Number(0);

    COVARIATES.forEach(function(c) {
      var stat = ee.Dictionary(statsDict.get(c));
      var sd = safeSd(stat.get('sd'));

      var pVal = safeNumber(projectFeature.get(c), 0);
      var dVal = safeNumber(donorFeature.get(c), 0);

      var z = pVal.subtract(dVal).divide(sd);
      d2 = d2.add(z.pow(2));
    });

    return donorFeature.set({
      temp_project_id: projectFeature.get('cell_id'),
      temp_match_d2: d2,
      temp_match_distance: d2.sqrt()
    });
  });

  var bestDonor = ee.Feature(
    scoredDonors.sort('temp_match_distance', true).first()
  );

  var projectCoords = projectFeature.geometry().coordinates();
  var donorCoords = bestDonor.geometry().coordinates();

  var projectLon = ee.Number(projectCoords.get(0));
  var projectLat = ee.Number(projectCoords.get(1));
  var donorLon = ee.Number(donorCoords.get(0));
  var donorLat = ee.Number(donorCoords.get(1));

  return ee.Feature(projectFeature.geometry()).set({
    scenario_name: SCENARIO_NAME,

    project_id: projectFeature.get('cell_id'),
    donor_id: bestDonor.get('cell_id'),

    project_lon: projectLon,
    project_lat: projectLat,
    donor_lon: donorLon,
    donor_lat: donorLat,

    match_d2: bestDonor.get('temp_match_d2'),
    match_distance: bestDonor.get('temp_match_distance'),

    project_eligible_area_ha: projectFeature.get('eligible_area_ha'),
    donor_eligible_area_ha: bestDonor.get('eligible_area_ha'),

    project_priority_eligible_area_ha: projectFeature.get('priority_eligible_area_ha'),
    donor_priority_eligible_area_ha: bestDonor.get('priority_eligible_area_ha'),

    project_broad_eligible_area_ha: projectFeature.get('broad_eligible_area_ha'),
    donor_broad_eligible_area_ha: bestDonor.get('broad_eligible_area_ha'),

    project_cropland_candidate_area_ha: projectFeature.get('cropland_candidate_area_ha'),
    donor_cropland_candidate_area_ha: bestDonor.get('cropland_candidate_area_ha'),

    project_feasibility_score_non_carbon: projectFeature.get('feasibility_score_non_carbon'),
    donor_feasibility_score_non_carbon: bestDonor.get('feasibility_score_non_carbon'),

    project_wc_priority_candidate_fraction: projectFeature.get('wc_priority_candidate_fraction'),
    donor_wc_priority_candidate_fraction: bestDonor.get('wc_priority_candidate_fraction'),

    project_wc_cropland_candidate_fraction: projectFeature.get('wc_cropland_candidate_fraction'),
    donor_wc_cropland_candidate_fraction: bestDonor.get('wc_cropland_candidate_fraction'),

    project_hard_excluded_fraction: projectFeature.get('hard_excluded_fraction'),
    donor_hard_excluded_fraction: bestDonor.get('hard_excluded_fraction'),

    project_forest2000_10pct_fraction: projectFeature.get('forest2000_10pct_fraction'),
    donor_forest2000_10pct_fraction: bestDonor.get('forest2000_10pct_fraction'),

    project_forest2000_30pct_fraction: projectFeature.get('forest2000_30pct_fraction'),
    donor_forest2000_30pct_fraction: bestDonor.get('forest2000_30pct_fraction'),

    project_recent_forest_loss_fraction: projectFeature.get('recent_forest_loss_fraction'),
    donor_recent_forest_loss_fraction: bestDonor.get('recent_forest_loss_fraction'),

    project_recent_forest_loss_pct: projectFeature.get('recent_forest_loss_pct'),
    donor_recent_forest_loss_pct: bestDonor.get('recent_forest_loss_pct'),

    project_pre_NDVI_mean: projectFeature.get('pre_NDVI_mean'),
    donor_pre_NDVI_mean: bestDonor.get('pre_NDVI_mean'),

    project_pre_EVI_mean: projectFeature.get('pre_EVI_mean'),
    donor_pre_EVI_mean: bestDonor.get('pre_EVI_mean'),

    project_pre_NDMI_mean: projectFeature.get('pre_NDMI_mean'),
    donor_pre_NDMI_mean: bestDonor.get('pre_NDMI_mean'),

    project_pre_NBR_mean: projectFeature.get('pre_NBR_mean'),
    donor_pre_NBR_mean: bestDonor.get('pre_NBR_mean'),

    project_pre_AGB_reference: projectFeature.get('pre_AGB_reference'),
    donor_pre_AGB_reference: bestDonor.get('pre_AGB_reference'),

    project_elevation: projectFeature.get('elevation'),
    donor_elevation: bestDonor.get('elevation'),

    project_slope: projectFeature.get('slope'),
    donor_slope: bestDonor.get('slope'),

    recent_loss_check_done_in_script: projectFeature.get('recent_loss_check_done_in_script'),
    recent_loss_period: projectFeature.get('recent_loss_period'),

    matching_grid_m: GRID_SIZE_M,
    max_project_points_used: MAX_PROJECT_POINTS,
    max_donor_points_used: MAX_DONOR_POINTS,
    covariates_used: COVARIATES.join(','),

    note_script2b1: 'Memory-safe matching only. Balance check will be done separately.'
  });
}


/***************************************************************
6. RUN MATCHING
****************************************************************/

var matchedPairs = projectReady.map(matchOneProject);


/***************************************************************
7. LIGHT DIAGNOSTICS
****************************************************************/

print('Matched pairs count', matchedPairs.size());
print('Matched pairs sample', matchedPairs.limit(5));


/***************************************************************
8. EXPORT MATCHED PAIRS
****************************************************************/

Export.table.toAsset({
  collection: matchedPairs,
  description: '03_Matched_Project_Control_Pairs_1000m_',
  assetId: OUTPUT_ASSET
});

Export.table.toDrive({
  collection: matchedPairs,
  description: '03_Matched_Project_Control_Pairs_1000m_' + '_CSV',
  fileNamePrefix: '03_Matched_Project_Control_Pairs_1000m_',
  fileFormat: 'CSV'
});
