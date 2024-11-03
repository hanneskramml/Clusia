

class RPi(object):
    MEASUREMENT_INTERVAL = 1  # every x minute(s)
    EXCEL_FILENAME = 'clusia_15072021.xlsx'
    ENV_FILE = '/home/pi/clusia/data/env.csv'
    LEAF_FILE = '/home/pi/clusia/data/leaf.csv'
    SOIL_FILE = '/home/pi/clusia/data/soil.csv'
    CLOUD_URL = 'https://ucloud.univie.ac.at/index.php/s/e01SkTTwaDvvJdz'
    CLOUD_PASSWORD = ''


class ArduinoDARK(object):
    NAME = 'ArduinoDARK'
    PORT = '/dev/ttyACM0'
    CONDITION = 'DARK'
    SENSORS = {
        # Environmental parameters
        "PAR2": [ "ENV", "Light", "" ],
        "DHTt": [ "ENV", "Temp", "" ],
        "DHTh": [ "ENV", "Humidity", "" ],
        # Infrared sensors
        "IR5l": [ "LEAF", "LeafTemp", "FD2" ],
        "IR5a": [ "LEAF", "AmbientTemp", "FD2" ],
        "IR6l": [ "LEAF", "LeafTemp", "MD1" ],
        "IR6a": [ "LEAF", "AmbientTemp", "MD1" ],
        "IR7l": [ "LEAF", "LeafTemp", "RD3" ],
        "IR7a": [ "LEAF", "AmbientTemp", "RD3" ],
        # Soil temperature
        "ST1B": ["SOIL", "SoilTemp", "MD1" ],
        "ST2B": ["SOIL", "SoilTemp", "FD2" ],
        "ST3B": ["SOIL", "SoilTemp", "RD3" ],
        "ST4B": ["SOIL", "SoilTemp", "MD6" ],
        # Soil moisture
        "MO10": ["SOIL", "SoilMoisture", "MD1" ],
        "MO11": ["SOIL", "SoilMoisture", "RD3" ],
        "MO12": ["SOIL", "SoilMoisture", "FD2" ],
        "MO13": ["SOIL", "SoilMoisture", "MD6" ],
        "MO14": ["SOIL", "SoilMoisture", "FD1" ],
        "MO15": ["SOIL", "SoilMoisture", "RD6" ],
        "MO16": ["SOIL", "SoilMoisture", "MD5" ],
        "MO17": ["SOIL", "SoilMoisture", "FD6" ],
        "MO18": ["SOIL", "SoilMoisture", "RD10" ]
    }


class ArduinoHL(object):
    NAME = 'ArduinoHL'
    PORT = '/dev/ttyACM1'
    CONDITION = 'HighLight'
    SENSORS = {
        # Environmental parameters
        "CO2x": [ "ENV", "CO2", "" ],
        "PAR1": [ "ENV", "Light", "" ],
        "DHTt": [ "ENV", "Temp", "" ],
        "DHTh": [ "ENV", "Humidity", "" ],
        # Infrared sensors
        "IR1l": [ "LEAF", "LeafTemp", "RH2" ],
        "IR1a": [ "LEAF", "AmbientTemp", "RH2" ],
        "IR2l": [ "LEAF", "LeafTemp", "MH1" ],
        "IR2a": [ "LEAF", "AmbientTemp", "MH1" ],
        "IR3l": [ "LEAF", "LeafTemp", "FH5" ],
        "IR3a": [ "LEAF", "AmbientTemp", "FH5" ],
        "IR4l": ["LEAF", "LeafTemp", "MH5"],
        "IR4a": ["LEAF", "AmbientTemp", "MH5"],
        # Soil temperature
        "ST1A": [ "SOIL", "SoilTemp", "RH2" ],
        "ST2A": [ "SOIL", "SoilTemp", "FH5" ],
        "ST3A": [ "SOIL", "SoilTemp", "MH1" ],
        "ST4A": [ "SOIL", "SoilTemp", "MH5" ],
        "ST5A": [ "SOIL", "SoilTemp", "FH2" ],
        "ST5B": [ "SOIL", "SoilTemp", "RH10" ],
        # Soil moisture
        "MO01": [ "SOIL", "SoilMoisture", "RH2" ],
        "MO02": [ "SOIL", "SoilMoisture", "FH5" ],
        "MO03": [ "SOIL", "SoilMoisture", "MH1" ],
        "MO04": [ "SOIL", "SoilMoisture", "MH5" ],
        "MO05": [ "SOIL", "SoilMoisture", "FH2" ],
        "MO06": [ "SOIL", "SoilMoisture", "RH10" ],
        "MO07": [ "SOIL", "SoilMoisture", "RH5" ],
        "MO08": [ "SOIL", "SoilMoisture", "MH6" ],
        "MO09": [ "SOIL", "SoilMoisture", "FH4" ]
    }
