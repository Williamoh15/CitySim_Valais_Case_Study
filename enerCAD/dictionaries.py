# year : WallType, RoofType, FloorType, Ninf, glazing_Uvalue, glazing_Gvalue, glazing_ratio
THRESHOLDS = {
    1945: (100, 200, 300, 1.4, 2.3, 0.47, 0.25),
    1960: (101, 201, 301, 1.3, 2.3, 0.47, 0.25),
    1970: (102, 202, 302, 1.2, 2.3, 0.47, 0.25),
    1980: (103, 203, 303, 1.1, 2.3, 0.47, 0.25),
    1990: (104, 204, 304, 1, 2.3, 0.47, 0.25),
    2000: (105, 205, 305, 0.8, 2.3, 0.47, 0.25),
    2010: (106, 206, 306, 0.7, 1.5, 0.48, 0.35),
    2015: (107, 207, 307, 0.7, 1.3, 0.53, 0.4),
    2023: (108, 207, 307, 0.35, 1.2, 0.53, 0.4),
}

# RegBL GBAUP : construction year threshold
PERIODS = {
    "8010": 1945,
    "8011": 1945,
    "8012": 1945,
    "8013": 1960,
    "8014": 1970,
    "8015": 1980,
    "8016": 1985,
    "8017": 1990,
    "8018": 1995,
    "8019": 2000,
    "8020": 2005,
    "8021": 2010,
    "8022": 2015,
    "8023": 2023,
}

# RegBL GKLAS : building type
TYPE = {
    "1110": 1, #Residential
    "1121": 1,
    "1122": 1,
    "1130": 1,
    "1275": 1,
    "1220": 2, #Office
    "1242": 3, #Garage
    "1230": 4, #Commercial
    "1231": 5, #Restaurant
    "1211": 6, #Hotel
    "1212": 6,
    "1264": 7, #Hospital
    "1261": 8, #Education
    "1262": 8,
    "1263": 8,
    "1251": 9, #Industrie
    "1241": 10, #Other
    "1252": 10, 
    "1271": 10,
    "1272": 10,
    "1273": 10,
    "1274": 10,
    "1276": 10,
    "1277": 10,
    "1278": 10,
    "1265": 11, #Sports installations
}

# building type : SIA Surface [m2]/person
SURFACE = {
    "1": 50, 
    "2": 20,
    "4": 10,
    "5": 5,
    "6": 40,
    "7": 30,
    "8": 10,
    "9": 20,
    "11": 20,
}

# building type : minimum ambient temperature [°C]
TEMPERATURE = {
    "1": 20, 
    "2": 20,
    "4": 20,
    "5": 20,
    "6": 20,
    "7": 22,
    "8": 20,
    "9": 18,
    "11": 18,
}