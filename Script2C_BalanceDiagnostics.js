/***************************************************************
FINAL SCRIPT 2C:
BALANCE CHECK ONLY FOR MATCHED PROJECT-CONTROL PAIRS

Input:
- 02B1_Matched_Project_Control_Pairs_1000m

Output:
- 02C_Matched_Project_Control_Balance_1000m

Purpose:
- Loads matched pairs from Script 2B1.
- Calculates sample-level balance diagnostics.
- Adds standardised difference of means for covariates.
- Adds pair-level absolute differences.

****************************************************************/


/***************************************************************
0. CONFIGURATION
****************************************************************/

var ASSET_ROOT = 'projects/ee-gdss2016/assets/Pre-feasibility_assessment/';

var MATCHED_INPUT_ASSET =
  ASSET_ROOT + '02B1_Matched_Project_Control_Pairs_1000m_';

var BALANCE_OUTPUT_ASSET =
  ASSET_ROOT + '02C_Matched_Project_Control_Balance_1000m_';

// Balance diagnostic threshold.
// Lower standardised difference is better.
// <= 0.25 is used here as a first-pass matching diagnostic.
var STANDARDISED_DIFF_THRESHOLD = 0.25;


/***************************************************************
1. LOAD MATCHED PAIRS
****************************************************************/

var matched = ee.FeatureCollection(MATCHED_INPUT_ASSET);

print('Running Final Script 2C balance check only');
print('Matched input asset', MATCHED_INPUT_ASSET);
print('Balance output asset', BALANCE_OUTPUT_ASSET);

print('Matched input count', matched.size());
print('Matched input sample', matched.limit(10));


/***************************************************************
2. COVARIATE PAIRS TO CHECK
****************************************************************/

var BALANCE_COVARIATES = [
  {
    name: 'NDVI',
    project: 'project_pre_NDVI_mean',
    donor: 'donor_pre_NDVI_mean'
  },
  {
    name: 'EVI',
    project: 'project_pre_EVI_mean',
    donor: 'donor_pre_EVI_mean'
  },
  {
    name: 'NDMI',
    project: 'project_pre_NDMI_mean',
    donor: 'donor_pre_NDMI_mean'
  },
  {
    name: 'NBR',
    project: 'project_pre_NBR_mean',
    donor: 'donor_pre_NBR_mean'
  },
  {
    name: 'AGB_reference',
    project: 'project_pre_AGB_reference',
    donor: 'donor_pre_AGB_reference'
  },
  {
    name: 'elevation',
    project: 'project_elevation',
    donor: 'donor_elevation'
  },
  {
    name: 'slope',
    project: 'project_slope',
    donor: 'donor_slope'
  },
  {
    name: 'eligible_area_ha',
    project: 'project_eligible_area_ha',
    donor: 'donor_eligible_area_ha'
  },
  {
    name: 'priority_eligible_area_ha',
    project: 'project_priority_eligible_area_ha',
    donor: 'donor_priority_eligible_area_ha'
  },
  {
    name: 'broad_eligible_area_ha',
    project: 'project_broad_eligible_area_ha',
    donor: 'donor_broad_eligible_area_ha'
  },
  {
    name: 'cropland_candidate_area_ha',
    project: 'project_cropland_candidate_area_ha',
    donor: 'donor_cropland_candidate_area_ha'
  },
  {
    name: 'priority_candidate_fraction',
    project: 'project_wc_priority_candidate_fraction',
    donor: 'donor_wc_priority_candidate_fraction'
  },
  {
    name: 'cropland_candidate_fraction',
    project: 'project_wc_cropland_candidate_fraction',
    donor: 'donor_wc_cropland_candidate_fraction'
  },
  {
    name: 'hard_excluded_fraction',
    project: 'project_hard_excluded_fraction',
    donor: 'donor_hard_excluded_fraction'
  },
  {
    name: 'forest2000_10pct_fraction',
    project: 'project_forest2000_10pct_fraction',
    donor: 'donor_forest2000_10pct_fraction'
  },
  {
    name: 'forest2000_30pct_fraction',
    project: 'project_forest2000_30pct_fraction',
    donor: 'donor_forest2000_30pct_fraction'
  },
  {
    name: 'recent_forest_loss_fraction',
    project: 'project_recent_forest_loss_fraction',
    donor: 'donor_recent_forest_loss_fraction'
  },
  {
    name: 'recent_forest_loss_pct',
    project: 'project_recent_forest_loss_pct',
    donor: 'donor_recent_forest_loss_pct'
  }
];


/***************************************************************
3. FILTER VALID PAIRS
****************************************************************/

function filterValidPairs(fc) {
  var out = fc;

  BALANCE_COVARIATES.forEach(function(c) {
    out = out
      .filter(ee.Filter.notNull([c.project]))
      .filter(ee.Filter.notNull([c.donor]));
  });

  return out;
}

var validMatched = filterValidPairs(matched);

print('Valid matched pairs count', validMatched.size());
print('Valid matched sample', validMatched.limit(10));


/***************************************************************
4. NULL-SAFE HELPERS
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

function flagNumber(condition) {
  return ee.Number(
    ee.Algorithms.If(
      condition,
      1,
      0
    )
  );
}


/***************************************************************
5. CALCULATE SAMPLE-LEVEL BALANCE METRICS
****************************************************************/

var balanceMetrics = ee.Dictionary({});

BALANCE_COVARIATES.forEach(function(c) {
  var pMean = safeNumber(validMatched.aggregate_mean(c.project), 0);
  var dMean = safeNumber(validMatched.aggregate_mean(c.donor), 0);

  var pSd = safeSd(validMatched.aggregate_total_sd(c.project));
  var dSd = safeSd(validMatched.aggregate_total_sd(c.donor));

  var pooledSd = pSd.pow(2)
    .add(dSd.pow(2))
    .divide(2)
    .sqrt();

  pooledSd = safeSd(pooledSd);

  var absMeanDiff = pMean.subtract(dMean).abs();

  var standardisedDiff = absMeanDiff.divide(pooledSd);

  var balancePass = flagNumber(
    standardisedDiff.lte(STANDARDISED_DIFF_THRESHOLD)
  );

  balanceMetrics = balanceMetrics
    .set('project_mean_' + c.name, pMean)
    .set('donor_mean_' + c.name, dMean)
    .set('project_sd_' + c.name, pSd)
    .set('donor_sd_' + c.name, dSd)
    .set('abs_mean_diff_' + c.name, absMeanDiff)
    .set('standardised_diff_' + c.name, standardisedDiff)
    .set('balance_pass_' + c.name, balancePass);
});


/***************************************************************
6. OVERALL BALANCE FLAG
****************************************************************/

var allBalancePassNumber = ee.Number(1);

BALANCE_COVARIATES.forEach(function(c) {
  var pass = ee.Number(
    balanceMetrics.get('balance_pass_' + c.name)
  );

  allBalancePassNumber = allBalancePassNumber.multiply(pass);
});

var overallBalancePass = allBalancePassNumber;


/***************************************************************
7. ADD PAIR-LEVEL ABSOLUTE DIFFERENCES
****************************************************************/

var balancedPairs = validMatched.map(function(f) {
  var out = f;

  BALANCE_COVARIATES.forEach(function(c) {
    var pVal = safeNumber(f.get(c.project), 0);
    var dVal = safeNumber(f.get(c.donor), 0);

    var absDiff = pVal.subtract(dVal).abs();

    out = out.set(
      'pair_abs_diff_' + c.name,
      absDiff
    );
  });

  return out
    .set(balanceMetrics)
    .set({
      standardised_diff_threshold: STANDARDISED_DIFF_THRESHOLD,
      overall_balance_pass: overallBalancePass,
      matched_sample_size_used: validMatched.size(),
      balance_check_script: 'Script 2C',
      note_script2c: 'Balance diagnostics for matched project-control pairs. No carbon calculated.'
    });
});


/***************************************************************
8. DIAGNOSTICS
****************************************************************/

print('Balanced pairs count', balancedPairs.size());
print('Balanced pairs sample', balancedPairs.limit(10));

print(
  'Overall balance pass numeric flag',
  overallBalancePass
);

print(
  'Standardised diff NDVI',
  balanceMetrics.get('standardised_diff_NDVI')
);

print(
  'Standardised diff EVI',
  balanceMetrics.get('standardised_diff_EVI')
);

print(
  'Standardised diff NDMI',
  balanceMetrics.get('standardised_diff_NDMI')
);

print(
  'Standardised diff NBR',
  balanceMetrics.get('standardised_diff_NBR')
);

print(
  'Standardised diff AGB reference',
  balanceMetrics.get('standardised_diff_AGB_reference')
);

print(
  'Standardised diff elevation',
  balanceMetrics.get('standardised_diff_elevation')
);

print(
  'Standardised diff slope',
  balanceMetrics.get('standardised_diff_slope')
);

print(
  'Standardised diff eligible area ha',
  balanceMetrics.get('standardised_diff_eligible_area_ha')
);

print(
  'Standardised diff priority candidate fraction',
  balanceMetrics.get('standardised_diff_priority_candidate_fraction')
);

print(
  'Standardised diff recent forest loss fraction',
  balanceMetrics.get('standardised_diff_recent_forest_loss_fraction')
);


/***************************************************************
9. EXPORT OUTPUT
****************************************************************/

Export.table.toAsset({
  collection: balancedPairs,
  description: '02C_Matched_Project_Control_Balance_1000m_',
  assetId: BALANCE_OUTPUT_ASSET
});

Export.table.toDrive({
  collection: balancedPairs,
  description: '02C_Matched_Project_Control_Balance_1000m_' + '_CSV',
  fileNamePrefix: '02C_Matched_Project_Control_Balance_1000m_',
  fileFormat: 'CSV'
});
