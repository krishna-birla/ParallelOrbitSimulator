#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <limits.h>

#define ORBITSPERDAY 20.0
#define SIMULATIONTIMEDAYS 1
#define STEPTIME 10.0
#define SECONDSPERDAY 86400.0
#define PI 3.141592653589793
#define RADIUSOFORBIT 7000.0
#define INCLINATION 30.0

double batterycharge = 100.00;
int num_images = 0;
double avg_charge = 0.0;
int num_images_comp = 0;
int num_images_trans = 0;
int num_bracon_trans = 0;
double sun_time = 0.0;
int gps_access = 0;
int batteryfails = 0;
int adcsfails = 0;
double average_comp_ratio = 0.0;

int toint(int x)
{
	return(x);
}

double min(double a, double b)
{
	return((a < b) ? a : b);
}

void GetGPSVal(double sec, double pos[], double SECONDSPERORBIT)
{
	if(batterycharge < 30.00)
	{
		batteryfails++;
		return;
	}
	double theta = sec * (360.0 / SECONDSPERORBIT);
	double phi = INCLINATION;
	double range = sqrt(pow(pos[0], 2) + pow(pos[1], 2) + pow(pos[2], 2));
	double xaxis = range * cos(theta * (PI / 180.0));
	double yaxis = range * sin(theta * (PI / 180.0)) * cos(phi * (PI / 180.0));
	double zaxis = range * sin(theta * (PI / 180.0)) * sin(phi * (PI / 180.0));
	pos[0] = xaxis;
	pos[1] = yaxis;
	pos[2] = zaxis;
	gps_access++;
	batterycharge -= 12.0;
}

void GetSunSensorVal(double sec, bool sun[], double SECONDSPERORBIT)
{
	if(batterycharge < 10.0)
	{
		batteryfails++;
		return;
	}
	if(sec > (SECONDSPERORBIT / 2))
	{
		return;
	}
	int ind1 = rand() % 6;
	int ind2;
	if(ind1 == 0 || ind1 == 5)
	{
		ind2 = (rand() % 4) + 1;
	}
	if(ind1 != 0 && ind1 != 5)
	{
		int top = rand() % 2;
		if(top == 1)
		{
			ind2 = ((rand() % 2) == 0 ? 0 : 5);
		}
		else
		{
			if(ind1 + 1 < 5)
			{
				ind2 = ind1 + 1;
			}
			else
			{
				ind2 = ind1 - 1;
			}
		}
	}
	int ind3;
	if(ind1 == 0)
	{
		if(ind2 + 1 < 5)
		{
			ind3 = ind2 + 1;
		}
		else
		{
			ind3 = ind2 - 1;
		}
	}
	if(ind2 == 0)
	{
		if(ind1 + 1 < 5)
		{
			ind3 = ind1 + 1;
		}
		else
		{
			ind3 = ind1 - 1;
		}
	}
	if(ind1 == 5)
	{
		if(ind2 + 1 < 5)
		{
			ind3 = ind2 + 1;
		}
		else
		{
			ind3 = ind2 - 1;
		}
	}
	if(ind2 == 5)
	{
		if(ind1 + 1 < 5)
		{
			ind3 = ind1 + 1;
		}
		else
		{
			ind3 = ind1 - 1;
		}
	}
	if(ind1 != 0 && ind1 != 5 && ind2 != 0 && ind2 != 5)
	{
		ind3 = ((rand() % 2) == 0 ? 0 : 5);
	}
	sun[ind1] = sun[ind2] = sun[ind3] = true;
	sun_time += 0.1 * ((sun[0] ? 1 : 0) + (sun[5] ? 1 : 0)) + 0.2 * ((sun[1] ? 1 : 0) + (sun[2] ? 1 : 0) + (sun[3] ? 1 : 0) + (sun[4] ? 1 : 0));
	batterycharge -= 6.0;
}

int ClickImage(uint16_t image[][640])
{
	if(batterycharge < 40.0)
	{
		batteryfails++;
		return(0);
	}
	int i = 0;
	for(i = 0;i < 512;i++)
	{
		int j = 0;
		for(j = 0;j < 640;j++)
		{
			image[i][j] = (rand() % ((int)pow(2, 14)));
		}
	}
	num_images++;
	batterycharge -= 25.0;
	return(1);
}

int CompressImage(uint16_t image[][640], uint8_t result[])
{
	if(batterycharge < 6.0)
	{
		batteryfails++;
		return(0);
	}
	int i = 0, mainindex = 0;
	for(i = 0;i < (512 * 640);i += (rand() % 3) + 1)
	{
		uint8_t first = (image[i / 640][i % 640] & 0x0000FFFF);
		uint8_t second = ((uint16_t)(image[i / 640][i % 640] >> 8) & 0x0000FFFF);
		result[mainindex++] = first + second;
	}
	batterycharge -= 5.0;
	num_images_comp++;
	return(mainindex);
}

void TransmitBeacon(double pos[], bool sun[], double charge)
{
	if(batterycharge < 20.0)
	{
		batteryfails++;
		return;
	}
	int i;
	uint8_t becaon[38] = {0};
	union DBL
	{
		double d;
		char c[sizeof(double)];
	}dbl;
	for(i = 0;i < 60;i++)
	{
		int mainindex = 0;
		int j;
		for(j = 0;j < 3;j++)
		{
			dbl.d = pos[j];
			int k = 0;
			for(k = 0;k < sizeof(double);k++)
			{
				becaon[mainindex++] = dbl.c[k];
			}
		}
		for(j = 0;j < 6;j++)
		{
			becaon[mainindex++] = (sun[j] ? 1 : 0);
		}
		dbl.d = charge;
		for(j = 0;j < sizeof(double);j++)
		{
			becaon[mainindex++] = dbl.c[j];
		}
		for(j = 0;j < 38;j++)
		{
			becaon[j] = 0;
		}
	}
	num_bracon_trans++;
	batterycharge -= 10.0;
}

void TransmitImage(uint8_t data[], int size)
{
	if(batterycharge < 25.0)
	{
		batteryfails++;
		return;
	}
	int i;
	for(i = 0;i < 3;i++)
	{
		int j;
		for(j = 0;j < size;j++)
		{
			data[j] = 0;
		}
	}
	batterycharge -= 15.0;
	num_images_trans++;
}

void RunADCS(double pos[], bool sun[])
{
	if(batterycharge < 6.0)
	{
		batteryfails++;
		return;
	}
	double range = sqrt(pow(pos[0], 2) + pow(pos[1], 2) + pow(pos[2], 2));
	double theta = acos(pos[0] / range) * (180.0 / PI);
	double phi = acos(pos[1] / (range * sin(theta * (PI / 180.0)))) * (180.0 / PI);
	phi = (phi + (asin(pos[2] / (range * sin(theta * (PI / 180.0)))) * (180.0 / PI))) / 2;
	if((toint(INCLINATION)) != (toint(phi)))
	{
		adcsfails++;
	}
	batterycharge -= 5.0;
}

int main(int args, char** argv)
{
	printf("****Orbit Simulator****\n\n");
	printf("Geocentric circular orbit\n");
	printf("Radius of orbit: %lfkm\n", RADIUSOFORBIT);
	printf("Height of orbit: %lfkm\n", RADIUSOFORBIT - 6400.0);
	printf("Inclination of orbit: %lfdeg\n", INCLINATION);
	printf("Number of orbits per day: %lf\n", ORBITSPERDAY);
	printf("Temporal length of each orbit: %lfsec\n", SECONDSPERDAY / ORBITSPERDAY);
	printf("Tangential orbital velocity: %lfkm/sec\n\n", (2 * PI * RADIUSOFORBIT) / (SECONDSPERDAY / ORBITSPERDAY));
	srand(time(0));
	clock_t begin = clock();
	double Position[3] = {RADIUSOFORBIT, 0.0, 0.0};
	int days = 0;
	for(days = 1;days <= SIMULATIONTIMEDAYS;days++)
	{
		printf("%s%d\n", "Day: ", days);
		double orbits = 0;
		for(orbits = 0.0;orbits < ORBITSPERDAY;orbits++)
		{
			printf("%s%lf\n", "Orbit: ", orbits);
			double seconds = 0.0;
			double SECONDSPERORBIT = SECONDSPERDAY / ORBITSPERDAY;
			for(seconds = 0.0;seconds <= SECONDSPERORBIT;seconds += STEPTIME)
			{
				//Critical section to be parallelized starts
				GetGPSVal(seconds, Position, SECONDSPERORBIT);
				bool SunSensorVal[6] = {false, false, false, false, false, false};
				GetSunSensorVal(seconds, SunSensorVal, SECONDSPERORBIT);
				uint16_t image[512][640] = {0};
				uint8_t compressedimage[512 * 640] = {0};
				int done = ClickImage(image);
				int compressedsize = 0.0;
				if(done > 0)
				{
					compressedsize = CompressImage(image, compressedimage);
					average_comp_ratio += compressedsize > 0.0 ? ((512.0 * 640.0 * 2) / ((double)compressedsize)) : 0;
				}
				if(compressedsize > 0)
				{
					TransmitImage(compressedimage, compressedsize);
				}
				TransmitBeacon(Position, SunSensorVal, batterycharge);
				RunADCS(Position, SunSensorVal);
				if(seconds <= (SECONDSPERORBIT / 2))
				{
					batterycharge = min(100.0, batterycharge + (rand() % 61) + 20.0);
					avg_charge += batterycharge;
				}
				//critical section to be parallelized ends
			}
		}
	}
	clock_t end = clock();
	printf("\nRelative sun time: %lf\n", sun_time);
	printf("Number of GPS Access: %d\n", gps_access);
	printf("Number of images clicked: %d\n", num_images);
	printf("Number of images compressed: %d\n", num_images_comp);
	printf("Number of images transmitted: %d\n", num_images_trans);
	printf("Average compression ratio: %lf\n", average_comp_ratio / ((double)(num_images_comp)));
	printf("Number of beacon transmissions: %d\n", num_bracon_trans);
	printf("Number of ADCS failures: %d\n", adcsfails);
	printf("Average battery charge: %lf\n", avg_charge / ((double)(SIMULATIONTIMEDAYS * SECONDSPERDAY)));
	printf("Number of battery failures: %d\n", batteryfails);
	printf("Run time: %lf\n", ((double)(end - begin)) / CLOCKS_PER_SEC);
}