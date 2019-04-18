#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <limits.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#define ORBITSPERDAY 20.0
#define SIMULATIONTIMEDAYS 1
#define STEPTIME 50.0
#define SECONDSPERDAY 86400.0
#define PI 3.141592653589793
#define RADIUSOFORBIT 7000.0
#define INCLINATION 30.0

double batterycharge = 100.00;
int num_images = 0;
double avg_charge = 0.0;
int num_images_comp = 0;
int num_images_trans = 0;
int num_beacon_trans = 0;
double sun_time = 0.0;
int gps_access = 0;
int batteryfails = 0;
int adcsfails = 0;
double average_comp_ratio = 0.0;

__global__ void masterkernel(int* imgsizecuda,double* pos, int* suncuda, uint16_t* imgcuda, uint8_t* compimgcuda, double* batterychargecuda,\
                                int* batteryfailscuda, double* avg_chargecuda, int* adcsfailscuda, int* gps_accesscuda,\
                                int* num_imagescuda, int* num_images_transcuda, int* num_images_compcuda,\
                                double* average_comp_ratiocuda, int* num_beacon_transcuda,\
                                double* sun_timecuda, double SECONDSPERORBIT, double inc, double sec)
{
    int bid = blockIdx.x;
    int tid = threadIdx.x;
    switch(bid)
    {
        case 0:{
            switch(tid)
            {
                case 0:{
                    if(*batterychargecuda < 30.00)
                    {
                        (*batteryfailscuda)++;
                        return;
                    }
                    double theta = sec * (360.0 / SECONDSPERORBIT);
                    double phi = inc;
                    double range = sqrt(pow(pos[0], 2) + pow(pos[1], 2) + pow(pos[2], 2));
                    double xaxis = range * cos(theta * (PI / 180.0));
                    double yaxis = range * sin(theta * (PI / 180.0)) * cos(phi * (PI / 180.0));
                    double zaxis = range * sin(theta * (PI / 180.0)) * sin(phi * (PI / 180.0));
                    pos[0] = xaxis;
                    pos[1] = yaxis;
                    pos[2] = zaxis;
                    (gps_accesscuda)++;
                    *batterychargecuda -= 12.0;
                    break;
                }
                case 1:{
                    if(*batterychargecuda < 10.0)
                    {
                        (*batteryfailscuda)++;
                        return;
                    }
                    if(sec > (SECONDSPERORBIT / 2))
                    {
                        return;
                    }
                    int ind1 = 0;
                    int ind2 = 3;
                    int ind3 = 4;
                    suncuda[ind1] = suncuda[ind2] = suncuda[ind3] = 1;
                    *sun_timecuda += 0.1 * ((suncuda[0] ? 1 : 0) + (suncuda[5] ? 1 : 0)) + 0.2 * ((suncuda[1] ? 1 : 0) + (suncuda[2] ? 1 : 0) + (suncuda[3] ? 1 : 0) + (suncuda[4] ? 1 : 0));
                    *batterychargecuda -= 6.0;
                    break;
                }
                case 2:{
                    if(*batterychargecuda < 40.0)
                    {
                        (*batteryfailscuda)++;
                        return;
                    }
                    int i = 0;
                    for(i = 0;i < 512;i++)
                    {
                        int j = 0;
                        for(j = 0;j < 640;j++)
                        {
                            imgcuda[i * 640 + j] = 5524;
                        }
                    }
                    (*num_imagescuda)++;
                    *batterychargecuda -= 25.0;
                    break;
                }
                case 3:{
                    if(*batterychargecuda < 6.0)
                    {
                        (*batteryfailscuda)++;
                        return;
                    }
                    int i = 0, mainindex = 0;
                    for(i = 0;i < (512 * 640);i += 2)
                    {
                        uint8_t first = (imgcuda[i] & 0x0000FFFF);
                        uint8_t second = (((uint16_t)(imgcuda[i] >> 8)) & 0x0000FFFF);
                        compimgcuda[mainindex++] = first + second;
                    }
                    *batterychargecuda -= 5.0;
                    (*num_images_compcuda)++;
                    *imgsizecuda = mainindex;
                    *average_comp_ratiocuda += (*imgsizecuda) > 0.0 ? ((512.0 * 640.0 * 2) / ((double)(*imgsizecuda))) : 0;
                    break;
                }
            }
            break;
        }
        case 1:{
            switch(tid)
            {
                case 0:{
                    if(*batterychargecuda < 25.0)
                    {
                        (*batteryfailscuda)++;
                        return;
                    }
                    int i;
                    for(i = 0;i < 3;i++)
                    {
                        int j;
                        for(j = 0;j < *imgsizecuda;j++)
                        {
                            compimgcuda[j] = 0;
                        }
                    }
                    *batterychargecuda -= 15.0;
                    (*num_images_compcuda)++;
                    break;
                }
                case 1:{
                    if(*batterychargecuda < 20.0)
                    {
                        (*batteryfailscuda)++;
                        return;
                    }
                    int i;
                    uint8_t becaon[38];
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
                            becaon[mainindex++] = (suncuda[j] ? 1 : 0);
                        }
                        dbl.d = *batterychargecuda;
                        for(j = 0;j < sizeof(double);j++)
                        {
                            becaon[mainindex++] = dbl.c[j];
                        }
                        for(j = 0;j < 38;j++)
                        {
                            becaon[j] = 0;
                        }
                    }
                    (*num_beacon_transcuda)++;
                    *batterychargecuda -= 10.0;
                    break;
                }
                case 2:{
                    if(*batterychargecuda < 6.0)
                    {
                        (*batteryfailscuda)++;
                        return;
                    }
                    double range = sqrt(pow(pos[0], 2) + pow(pos[1], 2) + pow(pos[2], 2));
                    double theta = acos(pos[0] / range) * (180.0 / PI);
                    double phi = acos(pos[1] / (range * sin(theta * (PI / 180.0)))) * (180.0 / PI);
                    phi = (phi + (asin(pos[2] / (range * sin(theta * (PI / 180.0)))) * (180.0 / PI))) / 2;
                    if(((int)(inc)) != ((int)(phi)))
                    {
                        (*adcsfailscuda)++;
                    }
                    *batterychargecuda -= 5.0;
                    break;
                }
                case 3:{
                    if(sec <= (SECONDSPERORBIT / 2))
                    {
                        *batterychargecuda = min(100.0, *batterychargecuda + 60.0);
                        *avg_chargecuda += *batterychargecuda;
                    }
                    break;
                }
            }
            break;
        }
    }
}

int main()
{
    printf("****Orbit Simulator****\n\n");
    printf("Geocentric circular orbit\n");
    printf("Radius of orbit: %lfkm\n", RADIUSOFORBIT);
    printf("Height of orbit: %lfkm\n", RADIUSOFORBIT - 6400.0);
    printf("Inclination of orbit: %lfdeg\n", INCLINATION);
    printf("Number of orbits per day: %lf\n", ORBITSPERDAY);
    printf("Temporal length of each orbit: %lfsec\n", SECONDSPERDAY / ORBITSPERDAY);
    printf("Tangential orbital velocity: %lfkm/sec\n\n", (2 * PI * RADIUSOFORBIT) / (SECONDSPERDAY / ORBITSPERDAY));
    double Position[3] = {RADIUSOFORBIT, 0.0, 0.0};
    uint16_t image[512][640] = {0};
    uint8_t compressedimage[512 * 640] = {0};
    int SunSensorVal[6] = {0};
    int days = 0;

    int *suncuda, *imgsizecuda, *num_imagescuda, *num_images_compcuda, *num_images_transcuda, *num_beacon_transcuda, *batteryfailscuda, *adcsfailscuda, *gps_accesscuda;
    uint8_t* compimgcuda;
    uint16_t* imgcuda;
    int compressedsize = 0.0;
    double *avg_chargecuda, *sun_timecuda, *batterychargecuda, *poscuda, *average_comp_ratiocuda;

    cudaMalloc(&poscuda, sizeof(double) * 3);
    cudaMalloc(&suncuda, sizeof(int) * 6);
    cudaMalloc(&imgcuda, sizeof(uint16_t) * 512 * 640);
    cudaMalloc(&compimgcuda, sizeof(uint8_t) * 512 * 640);
    cudaMalloc(&num_images_compcuda, sizeof(int));
    cudaMalloc(&num_imagescuda, sizeof(int));
    cudaMalloc(&num_images_transcuda, sizeof(int));
    cudaMalloc(&gps_accesscuda, sizeof(int));
    cudaMalloc(&num_beacon_transcuda, sizeof(int));
    cudaMalloc(&batteryfailscuda, sizeof(int));
    cudaMalloc(&imgsizecuda, sizeof(int));
    cudaMalloc(&adcsfailscuda, sizeof(int));
    cudaMalloc(&avg_chargecuda, sizeof(double));
    cudaMalloc(&average_comp_ratiocuda, sizeof(double));
    cudaMalloc(&batterychargecuda, sizeof(double));
    cudaMalloc(&sun_timecuda, sizeof(double));

    cudaMemcpy(poscuda, Position, sizeof(double) * 3, cudaMemcpyHostToDevice);
    cudaMemcpy(avg_chargecuda, &avg_charge, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(average_comp_ratiocuda, &average_comp_ratio, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(sun_timecuda, &sun_time, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(suncuda, SunSensorVal, sizeof(int) * 6, cudaMemcpyHostToDevice);
    cudaMemcpy(batterychargecuda, &batterycharge, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(batteryfailscuda, &batteryfails, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(imgsizecuda, &compressedsize, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(adcsfailscuda, &adcsfails, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(imgcuda, image, sizeof(uint16_t) * 512 * 640, cudaMemcpyHostToDevice);
    cudaMemcpy(compimgcuda, compressedimage, sizeof(uint8_t) * 512 * 640, cudaMemcpyHostToDevice);
    cudaMemcpy(gps_accesscuda, &gps_access, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(num_imagescuda, &num_images, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(num_beacon_transcuda, &num_beacon_trans, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(num_images_transcuda, &num_images_trans, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(num_images_compcuda, &num_images_comp, sizeof(int), cudaMemcpyHostToDevice);

    clock_t begin = clock();

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
                masterkernel<<<2, 4>>>(imgsizecuda, poscuda, suncuda, imgcuda,\
                                        compimgcuda, batterychargecuda, batteryfailscuda, avg_chargecuda,\
                                        adcsfailscuda, gps_accesscuda, num_imagescuda, num_images_transcuda,\
                                        num_images_compcuda, average_comp_ratiocuda, num_beacon_transcuda,\
                                        sun_timecuda, SECONDSPERORBIT, INCLINATION, seconds);
            }
        }
    }

    clock_t end = clock();

    cudaMemcpy(Position, poscuda, sizeof(double) * 3, cudaMemcpyDeviceToHost);
    cudaMemcpy(&avg_charge, avg_chargecuda, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&average_comp_ratio, average_comp_ratiocuda, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&sun_time, sun_timecuda, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(SunSensorVal, suncuda, sizeof(int) * 6, cudaMemcpyDeviceToHost);
    cudaMemcpy(&batterycharge, batterychargecuda, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&batteryfails, batteryfailscuda, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&adcsfails, adcsfailscuda, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&compressedsize, imgsizecuda, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(image, imgcuda, sizeof(uint16_t) * 512 * 640, cudaMemcpyDeviceToHost);
    cudaMemcpy(compressedimage, compimgcuda, sizeof(uint8_t) * 512 * 640, cudaMemcpyDeviceToHost);
    cudaMemcpy(&gps_access, gps_accesscuda, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&num_images, num_imagescuda, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&num_beacon_trans, num_beacon_transcuda, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&num_images_trans, num_images_transcuda, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&num_images_comp, num_images_compcuda, sizeof(int), cudaMemcpyDeviceToHost);

    printf("\nRelative sun time: %lf\n", sun_time);
    printf("Number of GPS Access: %d\n", gps_access);
    printf("Number of images clicked: %d\n", num_images);
    printf("Number of images compressed: %d\n", num_images_comp);
    printf("Number of images transmitted: %d\n", num_images_trans);
    printf("Average compression ratio: %lf\n", average_comp_ratio / ((double)(num_images_comp)));
    printf("Number of beacon transmissions: %d\n", num_beacon_trans);
    printf("Number of ADCS failures: %d\n", adcsfails);
    printf("Average battery charge: %lf\n", avg_charge / ((double)(SIMULATIONTIMEDAYS * SECONDSPERDAY)));
    printf("Number of battery failures: %d\n", batteryfails);
    printf("Run time: %lf\n", ((double)(end - begin)) / CLOCKS_PER_SEC);
}
