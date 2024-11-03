
#include "sensors.h"
#define RPi_BAUDRATE 115200
#define RPi_TIMEOUT 4

//Conditions
//#define HL
#define DARK

//Configuration for condition: HighLight
#ifdef HL

  #define CO2_RX_PIN 14             //CO2 sensor RX=>TX (MH-Z19C)
  #define CO2_TX_PIN 15             //CO2 sensor TX=>RX (MH-Z19C)
  
  #define DHT22_PIN 5               //Air temperature and humidity sensor (DHT22)
  #define SOILTEMP_BUS 16           //Soil temperature sensors (DS18B20)
  
  #define IR1_I2CADDR 0x5A          //Infrared sensor (Adafruit MLX90614)
  #define IR2_I2CADDR 0x5B
  #define IR3_I2CADDR 0x5C
  #define IR4_I2CADDR 0x5D
  
  #define SOILTEMP_NUM 6
  DeviceAddress SOILTEMP_ADDR_1A = { 0x28, 0x35, 0x97, 0xE2, 0x5C, 0x20, 0x01, 0xEF };
  DeviceAddress SOILTEMP_ADDR_2A = { 0x28, 0x65, 0x62, 0xB1, 0x5C, 0x20, 0x01, 0xE5 };
  DeviceAddress SOILTEMP_ADDR_3A = { 0x28, 0x69, 0x2D, 0xF2, 0x5C, 0x20, 0x01, 0x76 };
  DeviceAddress SOILTEMP_ADDR_4A = { 0x28, 0x9C, 0x31, 0xBE, 0x5C, 0x20, 0x01, 0x40 };
  DeviceAddress SOILTEMP_ADDR_5A = { 0x28, 0xA3, 0x80, 0xE1, 0x5C, 0x20, 0x01, 0x49 };
  DeviceAddress SOILTEMP_ADDR_5B = { 0x28, 0xFE, 0x51, 0xE8, 0x5C, 0x20, 0x01, 0x2D };
  
  #define MO1_PIN A0                //Soil moisture sensor
  #define MO2_PIN A1
  #define MO3_PIN A2
  #define MO4_PIN A3
  #define MO5_PIN A4
  #define MO6_PIN A5
  #define MO7_PIN A6
  #define MO8_PIN A7
  #define MO9_PIN A8

  CO2 *co2;
  Light *light;
  Dht22 *dht;
  
  IRThermo *ir1, *ir2, *ir3, *ir4;
  SoilTemps *sts;
  SoilMoisture *mo01, *mo02, *mo03, *mo04, *mo05, *mo06, *mo07, *mo08, *mo09;

#endif


//Configuration for condition: DARK
#ifdef DARK

  #define DHT22_PIN 5               //Air temperature and humidity sensor (DHT22)
  #define SOILTEMP_BUS 16           //Soil temperature sensors (DS18B20)
  
  #define IR5_I2CADDR 0x5A          //Infrared sensor (Adafruit MLX90614)
  #define IR6_I2CADDR 0x5B
  #define IR7_I2CADDR 0x5C

  #define SOILTEMP_NUM 4
  DeviceAddress SOILTEMP_ADDR_1B = { 0x28, 0xA4, 0xBA, 0x6D, 0x5C, 0x20, 0x01, 0x3E };
  DeviceAddress SOILTEMP_ADDR_2B = { 0x28, 0xA6, 0x9E, 0xD4, 0x5C, 0x20, 0x01, 0xB5 };
  DeviceAddress SOILTEMP_ADDR_3B = { 0x28, 0xDC, 0x0E, 0xE4, 0x5C, 0x20, 0x01, 0x1D };
  DeviceAddress SOILTEMP_ADDR_4B = { 0x28, 0xED, 0xAF, 0xE1, 0x5C, 0x20, 0x01, 0xB8 };
  
  #define MO10_PIN A0
  #define MO11_PIN A1
  #define MO12_PIN A2
  #define MO13_PIN A3
  #define MO14_PIN A4
  #define MO15_PIN A5
  #define MO16_PIN A6
  #define MO17_PIN A7
  #define MO18_PIN A8

  Light *light;
  Dht22 *dht;
  
  IRThermo *ir1, *ir2, *ir3;
  SoilTemps *sts;
  SoilMoisture *mo10, *mo11, *mo12, *mo13, *mo14, *mo15, *mo16, *mo17, *mo18;

#endif


void setup() {
  
  Serial.begin(RPi_BAUDRATE); 
  Serial.setTimeout(RPi_TIMEOUT);
  while (!Serial) {
    delay(1); // wait for serial port to connect...
  }
  Serial.println("Serial connection established! (Arduino=>RPi)");
  Serial.println("Setting up sensors...");

  #ifdef HL

    co2 = new CO2("CO2x", CO2_RX_PIN, CO2_TX_PIN);
    light = new Light("PAR1");
    dht = new Dht22("DHT", DHT22_PIN);
  
    ir1 = new IRThermo("IR1", IR1_I2CADDR);
    ir2 = new IRThermo("IR2", IR2_I2CADDR);
    ir3 = new IRThermo("IR3", IR3_I2CADDR);
    ir4 = new IRThermo("IR4", IR4_I2CADDR);
  
    sts = new SoilTemps("ST", SOILTEMP_NUM, SOILTEMP_BUS);
    sts->setAddr(0, "1A", SOILTEMP_ADDR_1A);
    sts->setAddr(1, "2A", SOILTEMP_ADDR_2A);
    sts->setAddr(2, "3A", SOILTEMP_ADDR_3A);
    sts->setAddr(3, "4A", SOILTEMP_ADDR_4A);
    sts->setAddr(4, "5A", SOILTEMP_ADDR_5A);
    sts->setAddr(5, "5B", SOILTEMP_ADDR_5B);
  
    mo01 = new SoilMoisture("MO01", MO1_PIN);
    mo02 = new SoilMoisture("MO02", MO2_PIN);
    mo03 = new SoilMoisture("MO03", MO3_PIN);
    mo04 = new SoilMoisture("MO04", MO4_PIN);
    mo05 = new SoilMoisture("MO05", MO5_PIN);
    mo06 = new SoilMoisture("MO06", MO6_PIN);
    mo07 = new SoilMoisture("MO07", MO7_PIN);
    mo08 = new SoilMoisture("MO08", MO8_PIN);
    mo09 = new SoilMoisture("MO09", MO9_PIN);

  #endif
  #ifdef DARK

    light = new Light("PAR2");
    dht = new Dht22("DHT", DHT22_PIN);
  
    ir1 = new IRThermo("IR5", IR5_I2CADDR);
    ir2 = new IRThermo("IR6", IR6_I2CADDR);
    ir3 = new IRThermo("IR7", IR7_I2CADDR);
  
    sts = new SoilTemps("ST", SOILTEMP_NUM, SOILTEMP_BUS);
    sts->setAddr(0, "1B", SOILTEMP_ADDR_1B);
    sts->setAddr(1, "2B", SOILTEMP_ADDR_2B);
    sts->setAddr(2, "3B", SOILTEMP_ADDR_3B);
    sts->setAddr(3, "4B", SOILTEMP_ADDR_4B);
  
    mo10 = new SoilMoisture("MO10", MO10_PIN);
    mo11 = new SoilMoisture("MO11", MO11_PIN);
    mo12 = new SoilMoisture("MO12", MO12_PIN);
    mo13 = new SoilMoisture("MO13", MO13_PIN);
    mo14 = new SoilMoisture("MO14", MO14_PIN);
    mo15 = new SoilMoisture("MO15", MO15_PIN);
    mo16 = new SoilMoisture("MO16", MO16_PIN);
    mo17 = new SoilMoisture("MO17", MO17_PIN);
    mo18 = new SoilMoisture("MO18", MO18_PIN);

  #endif

  Serial.println("Initialization completed!\n");
  
}

void loop() {

  Serial.println("Waiting for next measurement...");
  while (Serial.available() == 0);
  Serial.print(Serial.readString());

  Serial.println("\nMeasuring environmental parameters...");
  co2->measure();
  light->measure();
  dht->measure();

  Serial.println("\nMeasuring leaf temperatures...");
  ir1->measure();
  ir2->measure();
  ir3->measure();
  ir4->measure();


  Serial.println("\nMeasuring soil temperatures...");
  sts->measure();

  Serial.println("\nMeasuring soil moistures...");
  mo01->measure();
  mo02->measure();
  mo03->measure();
  mo04->measure();
  mo05->measure();
  mo06->measure();
  mo07->measure();
  mo08->measure();
  mo09->measure();

  Serial.println("\nMeasurement completed!");

}
