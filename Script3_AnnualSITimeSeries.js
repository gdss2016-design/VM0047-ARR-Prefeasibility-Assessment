/***************************************************************
FINAL SCRIPT 3:
ANNUAL SI TIME SERIES FOR MATCHED PROJECT-CONTROL PAIRS

Input:
- 02C_Matched_Project_Control_Balance_1000m_

Output:
- 03_Annual_SI_TimeSeries_MatchedPairs_1000m_

Purpose:
- Loads balanced matched project-control pairs.
- Extracts annual Landsat SI proxies:
  NDVI, EVI, NDMI, NBR.
- Calculates project-control annual differences.
- Uses project and donor coordinates exported from Script 2B1.
****************************************************************/


/***************************************************************
0. CONFIGURATION
****************************************************************/

var ASSET_ROOT = 'projects/ee-gdss2016/assets/Pre-feasibility_assessment/';

var MATCHED_INPUT_ASSET =
  ASSET_ROOT + '02C_Matched_Project_Control_Balance_1000m_';

var SI_OUTPUT_ASSET =
  ASSET_ROOT + '03_Annual_SI_TimeSeries_MatchedPairs_1000m_';

var START_YEAR = 2021;
var END_YEAR = 2025;

var GRID_SIZE_M = 1000;
var CELL_BUFFER_M = GRID_SIZE_M / 2;

var CRS = 'EPSG:32644';

var SEASON_START_MONTH = 10;
var SEASON_START_DAY = 1;
var SEASON_END_MONTH = 12;
var SEASON_END_DAY = 31;


var MAX_MATCHED_PAIRS = 50;


/***************************************************************
1. LOAD MATCHED PAIRS
****************************************************************/

var matchedAll = ee.FeatureCollection(MATCHED_INPUT_ASSET);

var matchedPairs = matchedAll
  .sort('match_distance', true)
  .limit(MAX_MATCHED_PAIRS);

print('Running Final Script 3: annual SI time series');
print('Matched input asset', MATCHED_INPUT_ASSET);
print('Output SI asset', SI_OUTPUT_ASSET);
print('Matched pairs used', matchedPairs.size());
print('Matched pair sample', matchedPairs.limit(10));


/***************************************************************
2. LANDSAT CLOUD MASK AND SCALE
****************************************************************/

function maskAndScaleLandsat(img) {
  var qa = img.select('QA_PIXEL');

  // Landsat Collection 2 QA bits:
  // 3 = cloud
  // 4 = cloud shadow
  // 5 = snow
  var cloudMask = qa.bitwiseAnd(1 << 3).eq(0);
  var shadowMask = qa.bitwiseAnd(1 << 4).eq(0);
  var snowMask = qa.bitwiseAnd(1 << 5).eq(0);

  var mask = cloudMask.and(shadowMask).and(snowMask);

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
3. SEASONAL SI IMAGE
****************************************************************/

function seasonalSIImage(year, analysisGeom) {
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
  ]).unmask(-9999);
}


/***************************************************************
4. SAFE NUMBER HELPER
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


/***************************************************************
5. EXTRACT SI FOR ONE PAIR AND ONE YEAR
****************************************************************/

function extractPairYear(pair, year) {
  pair = ee.Feature(pair);
  year = ee.Number(year);

  var projectPoint = ee.Geometry.Point([
    ee.Number(pair.get('project_lon')),
    ee.Number(pair.get('project_lat'))
  ]);

  var donorPoint = ee.Geometry.Point([
    ee.Number(pair.get('donor_lon')),
    ee.Number(pair.get('donor_lat'))
  ]);

  var projectGeom = projectPoint
    .transform(CRS, 1)
    .buffer(CELL_BUFFER_M, 1)
    .bounds(1, ee.Projection(CRS));

  var donorGeom = donorPoint
    .transform(CRS, 1)
    .buffer(CELL_BUFFER_M, 1)
    .bounds(1, ee.Projection(CRS));

  var combinedGeom = projectGeom
    .union(donorGeom, 1)
    .bounds(1, ee.Projection(CRS));

  var img = seasonalSIImage(year, combinedGeom);

  var projectStats = img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: projectGeom,
    scale: 30,
    maxPixels: 1e8,
    bestEffort: true,
    tileScale: 8
  });

  var donorStats = img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: donorGeom,
    scale: 30,
    maxPixels: 1e8,
    bestEffort: true,
    tileScale: 8
  });

  var pNDVI = safeNumber(projectStats.get('NDVI'), -9999);
  var pEVI = safeNumber(projectStats.get('EVI'), -9999);
  var pNDMI = safeNumber(projectStats.get('NDMI'), -9999);
  var pNBR = safeNumber(projectStats.get('NBR'), -9999);

  var dNDVI = safeNumber(donorStats.get('NDVI'), -9999);
  var dEVI = safeNumber(donorStats.get('EVI'), -9999);
  var dNDMI = safeNumber(donorStats.get('NDMI'), -9999);
  var dNBR = safeNumber(donorStats.get('NBR'), -9999);

  return ee.Feature(projectPoint).set({
    scenario_name: pair.get('scenario_name'),

    project_id: pair.get('project_id'),
    donor_id: pair.get('donor_id'),

    project_lon: pair.get('project_lon'),
    project_lat: pair.get('project_lat'),
    donor_lon: pair.get('donor_lon'),
    donor_lat: pair.get('donor_lat'),

    year: year,

    project_NDVI: pNDVI,
    donor_NDVI: dNDVI,
    delta_NDVI_project_minus_donor: pNDVI.subtract(dNDVI),

    project_EVI: pEVI,
    donor_EVI: dEVI,
    delta_EVI_project_minus_donor: pEVI.subtract(dEVI),

    project_NDMI: pNDMI,
    donor_NDMI: dNDMI,
    delta_NDMI_project_minus_donor: pNDMI.subtract(dNDMI),

    project_NBR: pNBR,
    donor_NBR: dNBR,
    delta_NBR_project_minus_donor: pNBR.subtract(dNBR),

    match_distance: pair.get('match_distance'),
    match_d2: pair.get('match_d2'),

    overall_balance_pass: pair.get('overall_balance_pass'),
    matched_sample_size_used: pair.get('matched_sample_size_used'),
    standardised_diff_threshold: pair.get('standardised_diff_threshold'),

    project_eligible_area_ha: pair.get('project_eligible_area_ha'),
    donor_eligible_area_ha: pair.get('donor_eligible_area_ha'),

    project_priority_eligible_area_ha: pair.get('project_priority_eligible_area_ha'),
    donor_priority_eligible_area_ha: pair.get('donor_priority_eligible_area_ha'),

    project_broad_eligible_area_ha: pair.get('project_broad_eligible_area_ha'),
    donor_broad_eligible_area_ha: pair.get('donor_broad_eligible_area_ha'),

    project_cropland_candidate_area_ha: pair.get('project_cropland_candidate_area_ha'),
    donor_cropland_candidate_area_ha: pair.get('donor_cropland_candidate_area_ha'),

    project_recent_forest_loss_pct: pair.get('project_recent_forest_loss_pct'),
    donor_recent_forest_loss_pct: pair.get('donor_recent_forest_loss_pct'),

    recent_loss_check_done_in_script: pair.get('recent_loss_check_done_in_script'),
    recent_loss_period: pair.get('recent_loss_period'),

    grid_size_m: GRID_SIZE_M,
    cell_buffer_m: CELL_BUFFER_M,
    season_start_month: SEASON_START_MONTH,
    season_end_month: SEASON_END_MONTH,

    note_script3: 'Annual SI proxy time series for matched project-control pairs.'
  });
}


/***************************************************************
6. BUILD ANNUAL TIME SERIES TABLE
****************************************************************/

var years = ee.List.sequence(START_YEAR, END_YEAR);

var pairCount = matchedPairs.size();
var pairList = matchedPairs.toList(pairCount);

var allRows = ee.FeatureCollection(
  ee.List.sequence(0, pairCount.subtract(1)).map(function(i) {
    var pair = ee.Feature(pairList.get(i));

    var rowsForPair = ee.FeatureCollection(
      years.map(function(y) {
        return extractPairYear(pair, y);
      })
    );

    return rowsForPair;
  })
).flatten();


/***************************************************************
7. DIAGNOSTICS
****************************************************************/

print('Annual SI rows count', allRows.size());
print('Annual SI sample', allRows.limit(10));

print(
  'Rows with invalid project NDVI',
  allRows.filter(ee.Filter.eq('project_NDVI', -9999)).size()
);

print(
  'Rows with invalid donor NDVI',
  allRows.filter(ee.Filter.eq('donor_NDVI', -9999)).size()
);


/***************************************************************
8. EXPORT OUTPUT
****************************************************************/

Export.table.toAsset({
  collection: allRows,
  description: '03_Annual_SI_TimeSeries_MatchedPairs_1000m_',
  assetId: SI_OUTPUT_ASSET
});

Export.table.toDrive({
  collection: allRows,
  description: '03_Annual_SI_TimeSeries_MatchedPairs_1000m_'+ '_CSV',
  fileNamePrefix: '03_Annual_SI_TimeSeries_MatchedPairs_1000m_',
  fileFormat: 'CSV'
});
