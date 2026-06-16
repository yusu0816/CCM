var MOD09A1 = ee.ImageCollection("MODIS/061/MOD09A1");
var MOD11A2 = ee.ImageCollection("MODIS/061/MOD11A2");
var MOD13Q1 = ee.ImageCollection("MODIS/061/MOD13Q1");
var roi = table;

var startYear = 2000;
var endYear = 2025;
var startMonth = 4;
var endMonth = 10;
var outScale = 1000;
var outCRS = 'EPSG:4326';

function maskMOD09A1(image) {
  var qa = image.select('QA');
  var modland = qa.bitwiseAnd(3).lte(1);
  var b1 = qa.rightShift(2).bitwiseAnd(15).lte(1);
  var b2 = qa.rightShift(6).bitwiseAnd(15).lte(1);
  var b3 = qa.rightShift(10).bitwiseAnd(15).lte(1);
  var b4 = qa.rightShift(14).bitwiseAnd(15).lte(1);
  var b5 = qa.rightShift(18).bitwiseAnd(15).lte(1);
  var b6 = qa.rightShift(22).bitwiseAnd(15).lte(1);
  var b7 = qa.rightShift(26).bitwiseAnd(15).lte(1);
  var atm = qa.rightShift(30).bitwiseAnd(1).eq(1);
  var mask = modland.and(b1).and(b2).and(b3).and(b4).and(b5).and(b6).and(b7).and(atm);
  return image.updateMask(mask);
}

var years = ee.List.sequence(startYear, endYear);

years.evaluate(function(yearList) {
  yearList.forEach(function(y) {
    var year = ee.Number(y);
    var startDate = ee.Date.fromYMD(year, startMonth, 1);
    var endDate = ee.Date.fromYMD(year, endMonth, 1).advance(1, 'month');

    var srIMG = MOD09A1.filterDate(startDate, endDate)
      .filterBounds(roi)
      .map(maskMOD09A1)
      .mean()
      .multiply(0.0001)
      .clip(roi);

    var mndwi = srIMG.normalizedDifference(['sur_refl_b04', 'sur_refl_b06']);
    var waterMask = mndwi.lte(0);

    var red = srIMG.select('sur_refl_b01');
    var nir = srIMG.select('sur_refl_b02');
    var blue = srIMG.select('sur_refl_b03');
    var green = srIMG.select('sur_refl_b04');
    var swir1 = srIMG.select('sur_refl_b06');
    var swir2 = srIMG.select('sur_refl_b07');

    var viIMG = MOD13Q1.filterDate(startDate, endDate)
      .filterBounds(roi)
      .mean()
      .clip(roi);

    var NDVI = viIMG.select('NDVI')
      .multiply(0.0001)
      .rename('NDVI')
      .updateMask(waterMask)
      .clip(roi);

    var EVI = viIMG.select('EVI')
      .multiply(0.0001)
      .rename('EVI')
      .updateMask(waterMask)
      .clip(roi);

    var kNDVI = NDVI.pow(2).tanh()
      .rename('kNDVI')
      .updateMask(waterMask)
      .clip(roi);

    var WET = srIMG.expression(
      '0.1147*B1 + 0.2489*B2 + 0.2408*B3 + 0.3132*B4 - 0.3122*B5 - 0.6416*B6 - 0.5087*B7', {
        'B1': red,
        'B2': nir,
        'B3': blue,
        'B4': green,
        'B5': srIMG.select('sur_refl_b05'),
        'B6': swir1,
        'B7': swir2
      }).rename('WET')
      .updateMask(waterMask)
      .clip(roi);

    var LST = MOD11A2.filterDate(startDate, endDate)
      .filterBounds(roi)
      .mean()
      .select('LST_Day_1km')
      .multiply(0.02)
      .subtract(273.15)
      .rename('LST')
      .updateMask(waterMask)
      .clip(roi);

    var si = swir1.add(red).subtract(nir.add(blue))
      .divide(swir1.add(red).add(nir).add(blue));

    var ibi_a = swir1.multiply(2).divide(swir1.add(nir));
    var ibi_b = nir.divide(nir.add(red)).add(green.divide(green.add(swir1)));
    var ibi = ibi_a.subtract(ibi_b).divide(ibi_a.add(ibi_b));

    var NDBSI = si.add(ibi).divide(2)
      .rename('NDBSI')
      .updateMask(waterMask)
      .clip(roi);

    var yearStr = y.toString();

    var exportTasks = [
      {img: NDVI, name: 'NDVI'},
      {img: EVI, name: 'EVI'},
      {img: kNDVI, name: 'kNDVI'},
      {img: WET, name: 'WET'},
      {img: NDBSI, name: 'NDBSI'},
      {img: LST, name: 'LST'}
    ];

    exportTasks.forEach(function(task) {
      Export.image.toDrive({
        image: task.img,
        description: task.name + '_' + yearStr,
        fileNamePrefix: task.name + '_' + yearStr,
        folder: task.name,
        scale: outScale,
        crs: outCRS,
        region: roi,
        maxPixels: 1e13
      });
    });
  });
});
