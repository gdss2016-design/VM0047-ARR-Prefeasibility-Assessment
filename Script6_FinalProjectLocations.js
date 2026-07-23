/***************************************************************
FINAL SCRIPT 6:
FINAL PROJECT LOCATIONS + TRUE FINAL AREA

Input:
- 05_Annual_SI_TimeSeries_MatchedPairs_1000m_

Output:
- 06_Final_Project_Locations_TrueArea_1000m_

Purpose:
- Loads annual SI time-series rows from Script 5.
- Keeps one row per final project location.
****************************************************************/


/***************************************************************
0. CONFIGURATION
****************************************************************/

var ASSET_ROOT = 'projects/ee-gdss2016/assets/Pre-feasibility_assessment/';

var INPUT_ASSET =
  ASSET_ROOT + '05_Annual_SI_TimeSeries_MatchedPairs_1000m_';

var OUTPUT_ASSET =
  ASSET_ROOT + '06_Final_Project_Locations_TrueArea_1000m_';


/***************************************************************
1. LOAD SCRIPT 3 OUTPUT
****************************************************************/

var annualRows = ee.FeatureCollection(INPUT_ASSET);

print('Running Final Script 6');
print('Input annual SI asset', INPUT_ASSET);
print('Output final locations asset', OUTPUT_ASSET);

print('Annual SI rows sample', annualRows.limit(10));
print('Annual SI row count', annualRows.size());


/***************************************************************
2. ADD LOCATION KEY
****************************************************************/

function addLocationKey(feature) {
  var centroid = feature.geometry().centroid(1);
  var coords = centroid.coordinates();

  var lon = ee.Number(coords.get(0));
  var lat = ee.Number(coords.get(1));

  var lonKey = lon.multiply(1000000).round().format();
  var latKey = lat.multiply(1000000).round().format();

  var locationKey = lonKey.cat('_').cat(latKey);

  return feature.set({
    final_project_lon: lon,
    final_project_lat: lat,
    location_key: locationKey
  });
}

var rowsWithKeys = annualRows.map(addLocationKey);


/***************************************************************
3. KEEP ONE ROW PER PROJECT LOCATION

Script 5 has repeated rows by year.
This keeps one project location only once.
****************************************************************/

var uniqueLocations = rowsWithKeys.distinct([
  'location_key'
]);

print('Unique final project locations', uniqueLocations.size());
print('Unique final locations sample', uniqueLocations.limit(10));


/***************************************************************
4. CREATE FINAL COMPACT OUTPUT
****************************************************************/

var finalLocations = uniqueLocations.map(function(f) {
  var areaHa = ee.Number(f.get('project_eligible_area_ha'));

  return ee.Feature(f.geometry(), {
    location_key: f.get('location_key'),

    project_id: f.get('project_id'),
    donor_id: f.get('donor_id'),

    project_lon: f.get('final_project_lon'),
    project_lat: f.get('final_project_lat'),
    donor_lon: f.get('donor_lon'),
    donor_lat: f.get('donor_lat'),

    final_eligible_area_ha: areaHa,
    final_eligible_area_sqkm: areaHa.divide(100),

    project_priority_eligible_area_ha: f.get('project_priority_eligible_area_ha'),
    project_broad_eligible_area_ha: f.get('project_broad_eligible_area_ha'),
    project_cropland_candidate_area_ha: f.get('project_cropland_candidate_area_ha'),

    donor_eligible_area_ha: f.get('donor_eligible_area_ha'),
    donor_priority_eligible_area_ha: f.get('donor_priority_eligible_area_ha'),
    donor_broad_eligible_area_ha: f.get('donor_broad_eligible_area_ha'),
    donor_cropland_candidate_area_ha: f.get('donor_cropland_candidate_area_ha'),

    project_recent_forest_loss_pct: f.get('project_recent_forest_loss_pct'),
    donor_recent_forest_loss_pct: f.get('donor_recent_forest_loss_pct'),

    recent_loss_check_done_in_script: f.get('recent_loss_check_done_in_script'),
    recent_loss_period: f.get('recent_loss_period'),

    match_distance: f.get('match_distance'),
    match_d2: f.get('match_d2'),

    overall_balance_pass: f.get('overall_balance_pass'),
    matched_sample_size_used: f.get('matched_sample_size_used'),
    standardised_diff_threshold: f.get('standardised_diff_threshold'),

    grid_size_m: f.get('grid_size_m'),

    area_method: 'Final area uses project_eligible_area_ha from Script 1 conservative ARR priority pixels, not nominal 100 ha grid area.',

    note_script4: 'Final project locations and true pixel-based eligible area. No carbon calculated. Recent loss check was already completed in Script 1.'
  });
});


/***************************************************************
5. SUMMARY DIAGNOSTICS
****************************************************************/

var totalFinalAreaHa = finalLocations.aggregate_sum(
  'final_eligible_area_ha'
);

var totalFinalAreaSqKm = ee.Number(totalFinalAreaHa).divide(100);

print('Final unique project locations', finalLocations.size());

print(
  'Total final eligible area ha',
  totalFinalAreaHa
);

print(
  'Total final eligible area sq km',
  totalFinalAreaSqKm
);

print(
  'Final locations sample',
  finalLocations.limit(20)
);

print(
  'Final area values sample',
  finalLocations
    .aggregate_array('final_eligible_area_ha')
    .slice(0, 30)
);

print(
  'Rows with recent loss check done in Script 1',
  finalLocations
    .filter(ee.Filter.eq('recent_loss_check_done_in_script', 'Script 1'))
    .size()
);

print(
  'Rows with final eligible area greater than 0',
  finalLocations
    .filter(ee.Filter.gt('final_eligible_area_ha', 0))
    .size()
);


/***************************************************************
6. MAP DISPLAY
****************************************************************/

Map.centerObject(finalLocations, 10);

Map.addLayer(
  finalLocations,
  {color: 'green'},
  '04 final project locations true area'
);


/***************************************************************
7. EXPORT OUTPUT
****************************************************************/

Export.table.toAsset({
  collection: finalLocations,
  description: '06_Final_Project_Locations_TrueArea_1000m_',
  assetId: OUTPUT_ASSET
});

Export.table.toDrive({
  collection: finalLocations,
  description: '06_Final_Project_Locations_TrueArea_1000m_' + '_CSV',
  fileNamePrefix: '06_Final_Project_Locations_TrueArea_1000m_',
  fileFormat: 'CSV'
});
