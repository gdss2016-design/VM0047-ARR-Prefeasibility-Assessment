/***************************************************************
FINAL SCRIPT 1:
CONSERVATIVE ARR PRIORITY AREA VERSION

GRID + ELIGIBILITY + 10-YEAR RECENT FOREST-LOSS CHECK
+ TRUE ELIGIBLE AREA

Outputs:
- 01_Project_EligibleScored_1000m
- 01_Donor_EligibleScored_1000m

Main logic:
- Conservative ARR priority area = shrubland + grassland + bare/sparse.
- Cropland is excluded from main eligible_area_ha.
- Cropland is exported separately as diagnostic area.
- Hansen lossyear and treecover2000 are unmasked to 0.
- 10-year recent forest-loss check is done here, before matching.
- selected_for_matching uses numeric 1/0 flag.
****************************************************************/


/***************************************************************
0. CONFIGURATION
****************************************************************/

var ASSET_ROOT = 'projects/ee-gdss2016/assets/Pre-feasibility_assessment/';

var PROJECT_BLOCK_ASSET = 'projects/ee-gdss2016/assets/Projectarea';
var DONOR_BLOCK_ASSET   = 'projects/ee-gdss2016/assets/Donorarea';

var PROJECT_START_YEAR = 2026;

var GRID_SIZE_M = 1000;
var CELL_HALF_M = GRID_SIZE_M / 2;

var CRS = 'EPSG:32644';
var PROJ_GRID = ee.Projection(CRS).atScale(GRID_SIZE_M);

// For project start year 2026, this becomes 2016 to 2025.
var LOSS_START_YEAR = PROJECT_START_YEAR - 10;
var LOSS_END_YEAR = PROJECT_START_YEAR - 1;

// Hansen lossyear code:
// 1 = 2001, 2 = 2002, ..., 25 = 2025.
var LOSS_YEAR_START_CODE = LOSS_START_YEAR - 2000;
var LOSS_YEAR_END_CODE = LOSS_END_YEAR - 2000;

// Screening thresholds.
var MAX_RECENT_LOSS_PCT_STRICT = 5;
var MAX_RECENT_LOSS_PCT_PRIORITY = 10;

// Main selected-cell rule.
var MIN_PRIORITY_ELIGIBLE_AREA_HA = 5;
var MIN_PRIORITY_CANDIDATE_FRACTION = 0.05;


/***************************************************************
1. OUTPUT ASSETS
****************************************************************/

var PROJECT_OUTPUT_ASSET =
  ASSET_ROOT + '01_Project_EligibleScored_1000m_';

var DONOR_OUTPUT_ASSET =
  ASSET_ROOT + '01_Donor_EligibleScored_1000m_';


/***************************************************************
2. LOAD BOUNDARIES
****************************************************************/

var projectBlock = ee.FeatureCollection(PROJECT_BLOCK_ASSET);
var donorBlock = ee.FeatureCollection(DONOR_BLOCK_ASSET);

Map.centerObject(projectBlock, 8);
Map.addLayer(projectBlock, {color: 'red'}, 'Project block');
Map.addLayer(donorBlock, {color: 'blue'}, 'Donor block');

print('Running Final Script 1 conservative ARR priority version');
print('Project block asset', PROJECT_BLOCK_ASSET);
print('Donor block asset', DONOR_BLOCK_ASSET);
print('Project block area ha', projectBlock.geometry().area(1).divide(10000));
print('Donor block area ha', donorBlock.geometry().area(1).divide(10000));
print('Recent forest-loss check period', LOSS_START_YEAR + '-' + LOSS_END_YEAR);
print('Minimum priority eligible area ha', MIN_PRIORITY_ELIGIBLE_AREA_HA);
print('Minimum priority candidate fraction', MIN_PRIORITY_CANDIDATE_FRACTION);
print('Project output', PROJECT_OUTPUT_ASSET);
print('Donor output', DONOR_OUTPUT_ASSET);


/***************************************************************
3. LOAD DATASETS
****************************************************************/

var worldCover = ee.ImageCollection('ESA/WorldCover/v200')
  .first()
  .select('Map')
  .unmask(0);

var hansen = ee.Image('UMD/hansen/global_forest_change_2025_v1_13');

// Critical fix:
// Hansen no-loss and no-tree pixels must be treated as 0, not masked.
var treecover2000 = hansen.select('treecover2000').unmask(0);
var lossYear = hansen.select('lossyear').unmask(0);


/***************************************************************
4. BASE LAND-COVER AND FOREST-LOSS MASKS
****************************************************************/

var wcTree = worldCover.eq(10).unmask(0).rename('wc_tree');
var wcShrub = worldCover.eq(20).unmask(0).rename('wc_shrub');
var wcGrass = worldCover.eq(30).unmask(0).rename('wc_grass');
var wcCrop  = worldCover.eq(40).unmask(0).rename('wc_crop');
var wcBuilt = worldCover.eq(50).unmask(0).rename('wc_built');
var wcBare  = worldCover.eq(60).unmask(0).rename('wc_bare');
var wcWater = worldCover.eq(80).unmask(0).rename('wc_water');
var wcWetland = worldCover.eq(90).unmask(0).rename('wc_wetland');
var wcMangrove = worldCover.eq(95).unmask(0).rename('wc_mangrove');

// Broad candidate area for diagnostics only.
// Includes cropland.
var wcCandidateBroad = wcShrub
  .or(wcGrass)
  .or(wcCrop)
  .or(wcBare)
  .unmask(0)
  .rename('wc_candidate_broad');

// Conservative ARR priority candidate area.
// Excludes cropland from main eligible area.
var wcPriorityARR = wcShrub
  .or(wcGrass)
  .or(wcBare)
  .unmask(0)
  .rename('wc_priority_arr');

var hardExcludedPixel = wcBuilt
  .or(wcWater)
  .or(wcWetland)
  .or(wcMangrove)
  .unmask(0)
  .rename('hard_excluded_pixel');

var forest2000_10pct = treecover2000
  .gte(10)
  .unmask(0)
  .rename('forest2000_10pct');

var forest2000_30pct = treecover2000
  .gte(30)
  .unmask(0)
  .rename('forest2000_30pct');

// Recent Hansen tree-cover loss in the 10 years before project start.
// For project start year 2026, this checks 2016 to 2025.
var recentForestLoss = lossYear
  .gte(LOSS_YEAR_START_CODE)
  .and(lossYear.lte(LOSS_YEAR_END_CODE))
  .unmask(0)
  .rename('recent_forest_loss');


/***************************************************************
5. PIXEL-LEVEL ELIGIBILITY AREA MASKS
****************************************************************/

var strictEligiblePixel = wcPriorityARR
  .and(forest2000_10pct.not())
  .and(recentForestLoss.not())
  .and(hardExcludedPixel.not())
  .unmask(0)
  .rename('strict_eligible_pixel');

var priorityEligiblePixel = wcPriorityARR
  .and(recentForestLoss.not())
  .and(hardExcludedPixel.not())
  .unmask(0)
  .rename('priority_eligible_pixel');

var broadEligiblePixel = wcCandidateBroad
  .and(recentForestLoss.not())
  .and(hardExcludedPixel.not())
  .unmask(0)
  .rename('broad_eligible_pixel');


/***************************************************************
6. POINT GRID CREATION
Includes edge-intersecting cells.
****************************************************************/

function makePointGrid(fc, prefix) {
  var geom = fc.geometry()
    .simplify(250)
    .transform(PROJ_GRID, 1);

  var sampleRegion = geom
    .buffer(CELL_HALF_M, 1)
    .bounds(1, PROJ_GRID);

  var coords = ee.Image.pixelCoordinates(PROJ_GRID)
    .reproject({
      crs: PROJ_GRID,
      scale: GRID_SIZE_M
    });

  var points = coords.sample({
    region: sampleRegion,
    projection: PROJ_GRID,
    scale: GRID_SIZE_M,
    geometries: true,
    tileScale: 16
  });

  return points.map(function(f) {
    var x = ee.Number(f.get('x'));
    var y = ee.Number(f.get('y'));

    return ee.Feature(f.geometry()).set({
      cell_id: ee.String(prefix)
        .cat('_')
        .cat(x.format('%.0f'))
        .cat('_')
        .cat(y.format('%.0f')),
      grid_source: prefix,
      grid_size_m: GRID_SIZE_M,
      nominal_area_ha: ee.Number(GRID_SIZE_M)
        .multiply(GRID_SIZE_M)
        .divide(10000),
      x_index: x,
      y_index: y
    });
  });
}


/***************************************************************
7. GEOMETRY AND SAFE NUMBER HELPERS
****************************************************************/

function makeCellSquare(feature) {
  return feature
    .geometry()
    .transform(PROJ_GRID, 1)
    .buffer(CELL_HALF_M, 1)
    .bounds(1, PROJ_GRID);
}

function safeDictNumber(dict, key, fallback) {
  dict = ee.Dictionary(dict);
  key = ee.String(key);

  return ee.Number(
    ee.Algorithms.If(
      dict.contains(key),
      dict.get(key),
      fallback
    )
  );
}

function safeFraction(numeratorHa, denominatorHa) {
  return ee.Number(
    ee.Algorithms.If(
      ee.Number(denominatorHa).gt(0),
      ee.Number(numeratorHa).divide(denominatorHa),
      0
    )
  );
}

function flagToNumber(condition) {
  return ee.Number(
    ee.Algorithms.If(
      condition,
      1,
      0
    )
  );
}


/***************************************************************
8. AREA IMAGE STACK HELPERS
****************************************************************/

function makeAreaStack10m(boundaryFc) {
  var boundaryMask = ee.Image.constant(1)
    .clip(boundaryFc.geometry())
    .selfMask();

  var boundaryArea = ee.Image.pixelArea()
    .updateMask(boundaryMask)
    .rename('boundary_area_m2');

  var broadCandidateArea = ee.Image.pixelArea()
    .updateMask(wcCandidateBroad)
    .updateMask(boundaryMask)
    .rename('broad_candidate_area_m2');

  var priorityCandidateArea = ee.Image.pixelArea()
    .updateMask(wcPriorityARR)
    .updateMask(boundaryMask)
    .rename('priority_candidate_area_m2');

  var croplandCandidateArea = ee.Image.pixelArea()
    .updateMask(wcCrop)
    .updateMask(boundaryMask)
    .rename('cropland_candidate_area_m2');

  var hardExcludedArea = ee.Image.pixelArea()
    .updateMask(hardExcludedPixel)
    .updateMask(boundaryMask)
    .rename('hard_excluded_area_m2');

  var treeArea = ee.Image.pixelArea()
    .updateMask(wcTree)
    .updateMask(boundaryMask)
    .rename('tree_area_m2');

  var builtArea = ee.Image.pixelArea()
    .updateMask(wcBuilt)
    .updateMask(boundaryMask)
    .rename('built_area_m2');

  var waterArea = ee.Image.pixelArea()
    .updateMask(wcWater)
    .updateMask(boundaryMask)
    .rename('water_area_m2');

  var wetlandArea = ee.Image.pixelArea()
    .updateMask(wcWetland)
    .updateMask(boundaryMask)
    .rename('wetland_area_m2');

  var mangroveArea = ee.Image.pixelArea()
    .updateMask(wcMangrove)
    .updateMask(boundaryMask)
    .rename('mangrove_area_m2');

  var strictEligibleArea = ee.Image.pixelArea()
    .updateMask(strictEligiblePixel)
    .updateMask(boundaryMask)
    .rename('strict_eligible_area_m2');

  var priorityEligibleArea = ee.Image.pixelArea()
    .updateMask(priorityEligiblePixel)
    .updateMask(boundaryMask)
    .rename('priority_eligible_area_m2');

  var broadEligibleArea = ee.Image.pixelArea()
    .updateMask(broadEligiblePixel)
    .updateMask(boundaryMask)
    .rename('broad_eligible_area_m2');

  return ee.Image.cat([
    boundaryArea,
    broadCandidateArea,
    priorityCandidateArea,
    croplandCandidateArea,
    hardExcludedArea,
    treeArea,
    builtArea,
    waterArea,
    wetlandArea,
    mangroveArea,
    strictEligibleArea,
    priorityEligibleArea,
    broadEligibleArea
  ]);
}

function makeAreaStack30m(boundaryFc) {
  var boundaryMask = ee.Image.constant(1)
    .clip(boundaryFc.geometry())
    .selfMask();

  var forest10Area = ee.Image.pixelArea()
    .updateMask(forest2000_10pct)
    .updateMask(boundaryMask)
    .rename('forest10_area_m2');

  var forest30Area = ee.Image.pixelArea()
    .updateMask(forest2000_30pct)
    .updateMask(boundaryMask)
    .rename('forest30_area_m2');

  var recentLossArea = ee.Image.pixelArea()
    .updateMask(recentForestLoss)
    .updateMask(boundaryMask)
    .rename('recent_loss_area_m2');

  return ee.Image.cat([
    forest10Area,
    forest30Area,
    recentLossArea
  ]);
}


/***************************************************************
9. SCORE ONE GRID COLLECTION
****************************************************************/

function scoreGrid(points, processingType, boundaryFc) {
  var areaStack10m = makeAreaStack10m(boundaryFc);
  var areaStack30m = makeAreaStack30m(boundaryFc);

  var scored = points.map(function(f) {
    var cellGeom = makeCellSquare(f);

    var stats10m = areaStack10m.reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: cellGeom,
      scale: 10,
      maxPixels: 1e8,
      bestEffort: true,
      tileScale: 8
    });

    var stats30m = areaStack30m.reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: cellGeom,
      scale: 30,
      maxPixels: 1e8,
      bestEffort: true,
      tileScale: 8
    });

    var trueCellAreaHa = safeDictNumber(stats10m, 'boundary_area_m2', 0).divide(10000);

    var broadCandidateAreaHa = safeDictNumber(stats10m, 'broad_candidate_area_m2', 0).divide(10000);
    var priorityCandidateAreaHa = safeDictNumber(stats10m, 'priority_candidate_area_m2', 0).divide(10000);
    var croplandCandidateAreaHa = safeDictNumber(stats10m, 'cropland_candidate_area_m2', 0).divide(10000);
    var hardExcludedAreaHa = safeDictNumber(stats10m, 'hard_excluded_area_m2', 0).divide(10000);

    var treeAreaHa = safeDictNumber(stats10m, 'tree_area_m2', 0).divide(10000);
    var builtAreaHa = safeDictNumber(stats10m, 'built_area_m2', 0).divide(10000);
    var waterAreaHa = safeDictNumber(stats10m, 'water_area_m2', 0).divide(10000);
    var wetlandAreaHa = safeDictNumber(stats10m, 'wetland_area_m2', 0).divide(10000);
    var mangroveAreaHa = safeDictNumber(stats10m, 'mangrove_area_m2', 0).divide(10000);

    var strictEligibleAreaHa = safeDictNumber(stats10m, 'strict_eligible_area_m2', 0).divide(10000);
    var priorityEligibleAreaHa = safeDictNumber(stats10m, 'priority_eligible_area_m2', 0).divide(10000);
    var broadEligibleAreaHa = safeDictNumber(stats10m, 'broad_eligible_area_m2', 0).divide(10000);

    var forest10AreaHa = safeDictNumber(stats30m, 'forest10_area_m2', 0).divide(10000);
    var forest30AreaHa = safeDictNumber(stats30m, 'forest30_area_m2', 0).divide(10000);
    var recentForestLossAreaHa = safeDictNumber(stats30m, 'recent_loss_area_m2', 0).divide(10000);

    var broadCandidateFraction = safeFraction(broadCandidateAreaHa, trueCellAreaHa);
    var priorityCandidateFraction = safeFraction(priorityCandidateAreaHa, trueCellAreaHa);
    var croplandCandidateFraction = safeFraction(croplandCandidateAreaHa, trueCellAreaHa);
    var hardExcludedFraction = safeFraction(hardExcludedAreaHa, trueCellAreaHa);

    var treeFraction = safeFraction(treeAreaHa, trueCellAreaHa);
    var builtFraction = safeFraction(builtAreaHa, trueCellAreaHa);
    var waterFraction = safeFraction(waterAreaHa, trueCellAreaHa);
    var wetlandFraction = safeFraction(wetlandAreaHa, trueCellAreaHa);
    var mangroveFraction = safeFraction(mangroveAreaHa, trueCellAreaHa);

    var forest2000_10pct_fraction = safeFraction(forest10AreaHa, trueCellAreaHa);
    var forest2000_30pct_fraction = safeFraction(forest30AreaHa, trueCellAreaHa);
    var recentForestLossFraction = safeFraction(recentForestLossAreaHa, trueCellAreaHa);

    var recentForestLossPct = recentForestLossFraction.multiply(100);

    var passMinPriorityAreaBool = priorityEligibleAreaHa.gte(
      MIN_PRIORITY_ELIGIBLE_AREA_HA
    );

    var passMinPriorityFractionBool = priorityCandidateFraction.gte(
      MIN_PRIORITY_CANDIDATE_FRACTION
    );

    var passForest30Bool = forest2000_30pct_fraction.lte(0.30);

    var passRecentLossBool = recentForestLossPct.lte(
      MAX_RECENT_LOSS_PCT_PRIORITY
    );

    var strictPassBool = priorityEligibleAreaHa.gte(MIN_PRIORITY_ELIGIBLE_AREA_HA)
      .and(priorityCandidateFraction.gte(MIN_PRIORITY_CANDIDATE_FRACTION))
      .and(forest2000_10pct_fraction.lte(0.20))
      .and(recentForestLossPct.lte(MAX_RECENT_LOSS_PCT_STRICT));

    var priorityPassBool = passMinPriorityAreaBool
      .and(passMinPriorityFractionBool)
      .and(passForest30Bool)
      .and(passRecentLossBool);

    var passMinPriorityArea = flagToNumber(passMinPriorityAreaBool);
    var passMinPriorityFraction = flagToNumber(passMinPriorityFractionBool);
    var passForest30 = flagToNumber(passForest30Bool);
    var passRecentLoss = flagToNumber(passRecentLossBool);

    var strictPass = flagToNumber(strictPassBool);
    var priorityPass = flagToNumber(priorityPassBool);

    var feasibilityScore = priorityCandidateFraction.multiply(3)
      .subtract(forest2000_10pct_fraction.multiply(1.5))
      .subtract(recentForestLossFraction.multiply(2))
      .subtract(hardExcludedFraction.multiply(2))
      .subtract(waterFraction.multiply(1.5))
      .subtract(wetlandFraction.multiply(1.5))
      .subtract(mangroveFraction.multiply(1.5))
      .subtract(croplandCandidateFraction.multiply(0.5));

    return f.set({
      processing_type: processingType,
      project_start_year: PROJECT_START_YEAR,

      recent_loss_check_done_in_script: 'Script 1',
      recent_loss_period: LOSS_START_YEAR + '-' + LOSS_END_YEAR,
      recent_lossyear_min_code: LOSS_YEAR_START_CODE,
      recent_lossyear_max_code: LOSS_YEAR_END_CODE,

      min_priority_eligible_area_ha: MIN_PRIORITY_ELIGIBLE_AREA_HA,
      min_priority_candidate_fraction: MIN_PRIORITY_CANDIDATE_FRACTION,

      true_cell_area_inside_boundary_ha: trueCellAreaHa,

      broad_candidate_area_ha: broadCandidateAreaHa,
      priority_candidate_area_ha: priorityCandidateAreaHa,
      cropland_candidate_area_ha: croplandCandidateAreaHa,
      hard_excluded_area_ha: hardExcludedAreaHa,

      wc_broad_candidate_fraction: broadCandidateFraction,
      wc_priority_candidate_fraction: priorityCandidateFraction,
      wc_cropland_candidate_fraction: croplandCandidateFraction,
      hard_excluded_fraction: hardExcludedFraction,

      pass_min_priority_area: passMinPriorityArea,
      pass_min_priority_fraction: passMinPriorityFraction,
      pass_forest2000_30pct: passForest30,
      pass_recent_loss_threshold: passRecentLoss,

      wc_tree_fraction: treeFraction,
      wc_built_fraction: builtFraction,
      wc_water_fraction: waterFraction,
      wc_wetland_fraction: wetlandFraction,
      wc_mangrove_fraction: mangroveFraction,

      forest2000_10pct_fraction: forest2000_10pct_fraction,
      forest2000_30pct_fraction: forest2000_30pct_fraction,

      recent_forest_loss_fraction: recentForestLossFraction,
      recent_forest_loss_pct: recentForestLossPct,
      recent_forest_loss_area_ha: recentForestLossAreaHa,

      strict_eligible_area_ha: strictEligibleAreaHa,
      priority_eligible_area_ha: priorityEligibleAreaHa,
      broad_eligible_area_ha: broadEligibleAreaHa,

      // Main area carried into Script 2 and final area reporting.
      eligible_area_ha: priorityEligibleAreaHa,

      strict_eligibility_pass: strictPass,
      priority_eligibility_pass: priorityPass,

      // Main numeric flag used in Script 2.
      selected_for_matching: priorityPass,

      feasibility_score_non_carbon: feasibilityScore,

      area_method: 'Main eligible_area_ha uses conservative ARR priority land: shrubland + grassland + bare/sparse vegetation, excluding cropland, recent loss, built-up, water, wetland, and mangrove. Selection uses numeric 1/0 flags for GEE filtering.',

      note_script1: 'Eligibility, 10-year recent forest-loss screening, and conservative ARR priority area calculated before matching. No carbon calculated.'
    });
  });

  return scored.filter(
    ee.Filter.gt('true_cell_area_inside_boundary_ha', 0)
  );
}


/***************************************************************
10. RUN PROJECT AND DONOR PROCESSING
****************************************************************/

var projectPoints = makePointGrid(projectBlock, 'P');
var donorPoints = makePointGrid(donorBlock, 'D');

var projectScored = scoreGrid(
  projectPoints,
  'project',
  projectBlock
);

var donorScored = scoreGrid(
  donorPoints,
  'donor',
  donorBlock
);


/***************************************************************
11. DIAGNOSTIC CHECKS
If memory issues occur, comment out this section only.
****************************************************************/

var projectSelected = projectScored.filter(
  ee.Filter.gt('selected_for_matching', 0)
);

print('Project selected cells sample', projectSelected.limit(10));

print(
  'Project total true boundary area in sampled cells ha',
  projectScored.aggregate_sum('true_cell_area_inside_boundary_ha')
);

print(
  'Project total broad candidate area ha',
  projectScored.aggregate_sum('broad_candidate_area_ha')
);

print(
  'Project total priority candidate area ha',
  projectScored.aggregate_sum('priority_candidate_area_ha')
);

print(
  'Project total cropland candidate area ha',
  projectScored.aggregate_sum('cropland_candidate_area_ha')
);

print(
  'Project total priority eligible area ha',
  projectScored.aggregate_sum('priority_eligible_area_ha')
);

print(
  'Project selected priority eligible area ha',
  projectSelected.aggregate_sum('eligible_area_ha')
);

print(
  'Project selected cell count',
  projectSelected.size()
);

print(
  'Project cells passing min priority area only',
  projectScored
    .filter(ee.Filter.gt('pass_min_priority_area', 0))
    .size()
);

print(
  'Project cells passing min priority fraction only',
  projectScored
    .filter(ee.Filter.gt('pass_min_priority_fraction', 0))
    .size()
);

print(
  'Project cells passing min area AND min fraction',
  projectScored
    .filter(ee.Filter.gt('pass_min_priority_area', 0))
    .filter(ee.Filter.gt('pass_min_priority_fraction', 0))
    .size()
);

print(
  'Project cells passing min area + min fraction + forest30',
  projectScored
    .filter(ee.Filter.gt('pass_min_priority_area', 0))
    .filter(ee.Filter.gt('pass_min_priority_fraction', 0))
    .filter(ee.Filter.gt('pass_forest2000_30pct', 0))
    .size()
);

print(
  'Project cells passing min area + min fraction + forest30 + recent loss',
  projectScored
    .filter(ee.Filter.gt('pass_min_priority_area', 0))
    .filter(ee.Filter.gt('pass_min_priority_fraction', 0))
    .filter(ee.Filter.gt('pass_forest2000_30pct', 0))
    .filter(ee.Filter.gt('pass_recent_loss_threshold', 0))
    .size()
);

print(
  'Project cells with priority eligible area greater than 0',
  projectScored.filter(
    ee.Filter.gt('priority_eligible_area_ha', 0)
  ).size()
);

print(
  'Project cells failing forest2000_30pct threshold',
  projectScored.filter(
    ee.Filter.gt('forest2000_30pct_fraction', 0.30)
  ).size()
);

print(
  'Project cells with recent forest loss greater than 10 percent',
  projectScored.filter(
    ee.Filter.gt('recent_forest_loss_pct', 10)
  ).size()
);


/***************************************************************
12. EXPORT OUTPUTS
****************************************************************/

Export.table.toAsset({
  collection: projectScored,
  description: '01_Project_EligibleScored_1000m_',
  assetId: PROJECT_OUTPUT_ASSET
});

Export.table.toDrive({
  collection: projectScored,
  description: '01_Project_EligibleScored_1000m_' + '_CSV',
  fileNamePrefix: '01_Project_EligibleScored_1000m_',
  fileFormat: 'CSV'
});

Export.table.toAsset({
  collection: donorScored,
  description: '01_Donor_EligibleScored_1000m_',
  assetId: DONOR_OUTPUT_ASSET
});

Export.table.toDrive({
  collection: donorScored,
  description: '01_Donor_EligibleScored_1000m_'+ '_CSV',
  fileNamePrefix: '01_Donor_EligibleScored_1000m_',
  fileFormat: 'CSV'
});
